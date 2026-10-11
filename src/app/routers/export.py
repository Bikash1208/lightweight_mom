"""Phase 13 Part 3, page 6 ("Export"): full-array STEP + layout CSV/JSON + HFSS-ready README,
chaining off a completed reflectarray (or HMSA) design job's cell assignments. Thin glue over
`export.array_export.export_array` (itself `phase12c.grounded_supercell.build_grounded_supercell`
+ `cad_export.export_supercell_step`, unmodified) -- grounded-array-only, the same limitation
`export_array`'s own module docstring states (square N x N grid).
"""

import json
import math
import threading
import time
import zipfile
from pathlib import Path

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

from ...errors import config_invalid
from ...export.aedt_script import generate_aedt_script
from ...export.array_export import export_array, export_multi_port_array
from ...export.board_readme import build_board_readme
from ...export.exact_board import (
    DEFAULT_MONO_H_MM,
    DEFAULT_R_PIN_MM,
    exact_board_from_ports_result,
    write_exact_step,
    write_feeds_csv,
    write_layout_csv,
    write_layout_dxf,
    write_layout_json,
)
from ...materials import Metal, Substrate
from ..jobs import JobStatus
from ..jsonutil import to_jsonable
from ..storage import job_output_dir

router = APIRouter(prefix="/api/export", tags=["export"])

# PHASE13D item 3: STEP export has no incremental progress hook (one blocking cadquery call) --
# these are Phase 12D's own MEASURED runtimes (docs/USER_GUIDE.md's Verify page costing, the same
# CAD boolean-union cost family) log-log interpolated/extrapolated to estimate a percentage while
# the blocking call runs, so the UI shows *something* better than a silent freeze. Labelled
# "estimated" wherever shown -- never presented as a measured number.
_STEP_BENCHMARK_N = (6, 8, 12)
_STEP_BENCHMARK_S = (7.0, 18.5, 124.0)


def estimate_step_export_seconds(n_cells_per_side):
    """Log-log linear interpolation/extrapolation across `_STEP_BENCHMARK_N/_S`."""
    xs = [math.log(v) for v in _STEP_BENCHMARK_N]
    ys = [math.log(v) for v in _STEP_BENCHMARK_S]
    x = math.log(max(n_cells_per_side, 1))
    if x <= xs[0]:
        i = 0
    elif x >= xs[-1]:
        i = len(xs) - 2
    else:
        i = next(k for k in range(len(xs) - 1) if xs[k] <= x <= xs[k + 1])
    t = (x - xs[i]) / (xs[i + 1] - xs[i])
    return float(math.exp(ys[i] + t * (ys[i + 1] - ys[i])))


class _StepProgressEstimator:
    """Ticks `job.progress`/`job.progress_message` toward (but never reaching) 1.0 over the
    estimated duration while the blocking STEP export call runs in the main worker thread."""

    def __init__(self, job, estimated_s, n_cells_per_side):
        self._job = job
        self._estimated_s = max(estimated_s, 1.0)
        self._n = n_cells_per_side
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def _run(self):
        t0 = time.perf_counter()
        while not self._stop.wait(0.5):
            frac = min(0.95, (time.perf_counter() - t0) / self._estimated_s)
            self._job.progress = frac
            self._job.progress_message = (
                f"Building STEP (CAD) for {self._n}x{self._n} cells -- "
                f"~{self._estimated_s:.0f}s estimated, {frac * 100:.0f}%"
            )

    def __enter__(self):
        self._job.progress_message = (
            f"Building STEP (CAD) for {self._n}x{self._n} cells -- ~{self._estimated_s:.0f}s estimated"
        )
        self._thread.start()
        return self

    def __exit__(self, *exc_info):
        self._stop.set()
        self._thread.join(timeout=2.0)
        self._job.progress = 1.0
        self._job.progress_message = "STEP export complete."


class ExportRequest(BaseModel):
    source_job_id: str
    period_mm: float
    eps_r: float
    tan_delta: float
    thickness_mm: float
    include_step: bool = False  # PHASE13D item 3: STEP is opt-in -- can be 100+ MB and take minutes


def _lw_grid_from_assignments(assignments, nx, ny):
    grid = [[None] * nx for _ in range(ny)]
    for a in assignments:
        if a.chosen_l_mm is None or a.chosen_w_mm is None:
            raise config_invalid(
                message=f"cell ({a.i},{a.j}) has no library match -- cannot export an incomplete array"
            )
        grid[a.i][a.j] = (a.chosen_l_mm * 1e-3, a.chosen_w_mm * 1e-3)
    return grid


def _download(job_id, path):
    if path is None:
        return None
    return f"/api/files/{job_id}/{Path(path).name}"


def _check_lw_grid_complete(lw_grid_m):
    for i, row in enumerate(lw_grid_m):
        for j, cell in enumerate(row):
            if cell is None:
                raise config_invalid(message=f"cell ({i},{j}) has no geometry -- cannot export an incomplete array")


def _run_export_job(job, source_result, params):
    substrate = Substrate(
        name="user", eps_r=params["eps_r"], tan_delta=params["tan_delta"], thickness=params["thickness_mm"] * 1e-3
    )
    metal = Metal()
    period_m = params["period_mm"] * 1e-3
    out_dir = job_output_dir(job.id)
    # PHASE15 Part C: mark each feed location in the STEP/README, if the source design recorded
    # one (every pre-Phase-15 synthesize result has no "feeds_mm" key -- .get(...) keeps this
    # export byte-for-byte unaffected for those).
    feeds_mm = source_result.get("export_feeds_mm", source_result.get("feeds_mm"))

    if source_result.get("assignments") is not None:
        assignments = source_result["assignments"]
        ny, nx = source_result["x_m"].shape
        lw_grid_m = _lw_grid_from_assignments(assignments, nx, ny)
        if params["include_step"]:
            n = len(lw_grid_m)
            with _StepProgressEstimator(job, estimate_step_export_seconds(n), n):
                result = export_array(
                    lw_grid_m, period_m, period_m, substrate, metal, out_dir, include_step=True, feeds_mm=feeds_mm
                )
        else:
            result = export_array(
                lw_grid_m, period_m, period_m, substrate, metal, out_dir, include_step=False, feeds_mm=feeds_mm
            )
        return {
            "n_cells_per_side": result.n_cells_per_side,
            "total_cells": result.total_cells,
            "file_size_bytes": result.file_size_bytes,
            "export_time_s": result.export_time_s,
            "step_download_url": _download(job.id, result.step_path),
            "layout_csv_download_url": _download(job.id, result.layout_csv_path),
            "layout_json_download_url": _download(job.id, result.layout_json_path),
            "layout_dxf_download_url": _download(job.id, result.layout_dxf_path),
            "readme_download_url": _download(job.id, result.readme_path),
        }

    # PHASE15 Part C/E: a multi-port board (possibly non-square, `"ports"` present) or a 1-bit
    # Coded design (possibly mixing plain-rectangle and annular-slot cells) -- both hand back a
    # ready-made `lw_grid_m` directly (metres, mixed-cell-safe) rather than `assignments`
    # (reflectarray-shaped, rectangle-only); `export_multi_port_array` handles either.
    lw_grid_m = source_result.get("lw_grid_m")
    if lw_grid_m is None:
        raise config_invalid(message="source job produced neither 'assignments' nor 'lw_grid_m' -- cannot export")
    _check_lw_grid_complete(lw_grid_m)
    ports_info = None
    if "ports" in source_result:
        ports_info = [(p["port_id"], p["feed_x_mm"], p["feed_y_mm"]) for p in source_result["ports"]]
        if feeds_mm and len(feeds_mm) == len(ports_info):
            ports_info = [(pid, fx, fy) for (pid, _x, _y), (fx, fy) in zip(ports_info, feeds_mm, strict=True)]
    n = max(len(lw_grid_m), len(lw_grid_m[0]) if lw_grid_m else 0)
    if params["include_step"]:
        with _StepProgressEstimator(job, estimate_step_export_seconds(n), n):
            result = export_multi_port_array(
                lw_grid_m,
                period_m,
                period_m,
                substrate,
                metal,
                out_dir,
                include_step=True,
                feeds_mm=feeds_mm,
                ports_info=ports_info,
            )
    else:
        result = export_multi_port_array(
            lw_grid_m,
            period_m,
            period_m,
            substrate,
            metal,
            out_dir,
            include_step=False,
            feeds_mm=feeds_mm,
            ports_info=ports_info,
        )
    return {
        "n_cells_x": result.n_cells_x,
        "n_cells_y": result.n_cells_y,
        "total_cells": result.total_cells,
        "file_size_bytes": result.file_size_bytes,
        "export_time_s": result.export_time_s,
        "step_download_url": _download(job.id, result.step_path),
        "layout_csv_download_url": _download(job.id, result.layout_csv_path),
        "layout_json_download_url": _download(job.id, result.layout_json_path),
        "layout_dxf_download_url": _download(job.id, result.layout_dxf_path),
        "feeds_csv_download_url": _download(job.id, result.feeds_csv_path) if result.feeds_csv_path else None,
        "readme_download_url": _download(job.id, result.readme_path),
    }


@router.post("/array")
def export_reflectarray(body: ExportRequest, request: Request):
    source = request.app.state.jobs.get(body.source_job_id)
    if source is None:
        raise HTTPException(404, f"job {body.source_job_id!r} not found")
    if source.status != JobStatus.DONE:
        raise HTTPException(409, f"source job {body.source_job_id!r} is {source.status.value}, not done")
    result = source.result or {}
    if result.get("assignments") is None and result.get("lw_grid_m") is None:
        raise HTTPException(422, f"source job {body.source_job_id!r} did not produce cell assignments")

    job = request.app.state.jobs.submit("array_export", _run_export_job, source.result, body.model_dump())
    return job.to_dict()


class BoardExportRequest(BaseModel):
    """Phase 21 D: exact export of a multi-port board (STEP from exact vertex coordinates, layout CSV,
    AEDT/HFSS script, README) from a completed `/synthesize_ports` job."""

    source_job_id: str
    period_mm: float
    eps_r: float
    tan_delta: float
    thickness_mm: float
    frequency_hz: float
    include_step: bool = True
    include_pins: bool = False  # add PEC monopole pins (cylinders) to the STEP
    r_pin_mm: float = float(DEFAULT_R_PIN_MM)
    mono_h_mm: float = float(DEFAULT_MONO_H_MM)
    name: str = "HMSA_board"


def _write_board_png(board, path):
    """Patch layout with every feed marked and labelled (best effort: a plotting failure never fails the export)."""
    try:
        from ..plots import patch_layout_png

        period = float(board.period_mm)
        cells = [(float((p.x0 + p.x1) / 2), float((p.y0 + p.y1) / 2), float(p.x1 - p.x0), float(p.y1 - p.y0)) for p in board.patches]
        nx = max(1, round(float(board.width_mm) / period))
        ny = max(1, round(float(board.height_mm) / period))
        return patch_layout_png(
            cells, period, nx, ny, path, title=f"{len(board.feeds_mm)}-port board, {len(board.patches)} patches",
            feed_xy_mm=[(float(x), float(y)) for x, y in board.feeds_mm],
        )
    except Exception:
        return None


def _run_board_export_job(job, source_result, params):
    out_dir = job_output_dir(job.id)
    board = exact_board_from_ports_result(source_result, params["period_mm"], params["thickness_mm"])
    files = []
    step_path = None
    if params["include_step"]:
        job.progress_message = "Writing the exact STEP..."
        step_path = write_exact_step(
            board,
            out_dir / "board.step",
            include_pins=params["include_pins"],
            r_pin_mm=params["r_pin_mm"],
            mono_h_mm=params["mono_h_mm"],
        )
        files.append("board.step")
    csv_path = write_layout_csv(board, out_dir / "layout.csv")
    feeds_path = write_feeds_csv(board, out_dir / "feeds.csv")
    json_path = write_layout_json(board, out_dir / "layout.json")
    dxf_path = write_layout_dxf(board, out_dir / "layout.dxf")
    files += ["layout.csv", "feeds.csv", "layout.json", "layout.dxf"]
    png_path = _write_board_png(board, out_dir / "layout.png")
    if png_path:
        files.append("layout.png")
    script_path = out_dir / "hfss_build_board.py"
    script_path.write_text(
        generate_aedt_script(board, params["name"], params["eps_r"], params["tan_delta"], params["frequency_hz"])
    )
    files.append("hfss_build_board.py")
    readme_path = out_dir / "README.md"
    files.append("README.md")
    readme_path.write_text(
        build_board_readme(
            board, source_result, params["name"], params["eps_r"], params["tan_delta"], params["frequency_hz"], files,
            with_pins=bool(params["include_step"] and params["include_pins"]),
        )
    )
    return {
        "n_patches": len(board.patches),
        "board_width_mm": float(board.width_mm),
        "board_height_mm": float(board.height_mm),
        "max_rounding_mm": board.max_rounding_mm,
        "step_download_url": _download(job.id, step_path),
        "layout_csv_download_url": _download(job.id, csv_path),
        "feeds_csv_download_url": _download(job.id, feeds_path),
        "layout_json_download_url": _download(job.id, json_path),
        "layout_dxf_download_url": _download(job.id, dxf_path),
        "layout_png_download_url": _download(job.id, png_path) if png_path else None,
        "has_monopole_pins": bool(params["include_step"] and params["include_pins"]),
        "hfss_script_download_url": _download(job.id, script_path),
        "readme_download_url": _download(job.id, readme_path),
    }


@router.post("/board")
def export_board(body: BoardExportRequest, request: Request):
    source = request.app.state.jobs.get(body.source_job_id)
    if source is None:
        raise HTTPException(404, f"job {body.source_job_id!r} not found")
    if source.status != JobStatus.DONE:
        raise HTTPException(409, f"source job {body.source_job_id!r} is {source.status.value}, not done")
    ports = (source.result or {}).get("ports")
    if not ports or "patches_mm" not in ports[0]:
        raise HTTPException(422, f"source job {body.source_job_id!r} is not a multi-port design")
    job = request.app.state.jobs.submit("board_export", _run_board_export_job, source.result, body.model_dump())
    return job.to_dict()


def _download_url_to_path(download_url):
    """`download_url`s this app hands out are always `/api/files/{job_id}/{filename}` -- resolve
    straight back to the file `files.py`'s route would serve, without an HTTP round-trip."""
    if not download_url:
        return None
    job_id, filename = download_url.rsplit("/", 2)[-2:]
    return job_output_dir(job_id) / filename


class ExportPackageRequest(BaseModel):
    export_job_id: str
    pattern_job_id: "str | None" = None  # optional: bundles the pattern PNGs too
    synthesize_job_id: "str | None" = None  # optional: bundles the layout PNG (PHASE13D item 4 --
    # the layout preview is a plot the HMSA/reflectarray synthesis job produces, not the export
    # job; the export result never had a `layout_png_url` key to read, so the zip silently never
    # got layout.png until this was wired through)


def _run_package_job(job, export_result, pattern_result, synthesize_result):
    """13C/13D Step 3: "one 'Download package' (zip) ... the project record (inputs + results, to
    reload later)". The project record here captures every RESULT value already computed
    (export summary, pattern summary/plots) -- not the full original request bodies (beam list,
    board substrate, ...), which the app does not currently persist anywhere the export job could
    read back; a documented simplification, not a silent one (docs/history/PHASE13_WEBAPP.md)."""
    out_dir = job_output_dir(job.id)
    zip_path = out_dir / "package.zip"
    record = {
        "export": to_jsonable(export_result),
        "pattern": to_jsonable(pattern_result),
        "synthesis": to_jsonable(synthesize_result),
    }
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for key, arcname in (
            ("step_download_url", "array.step"),
            ("layout_csv_download_url", "layout.csv"),
            ("layout_json_download_url", "layout.json"),
            ("layout_dxf_download_url", "layout.dxf"),
            ("readme_download_url", "HFSS_README.md"),
            ("hfss_script_download_url", "hfss_build_board.py"),
            ("feeds_csv_download_url", "feeds.csv"),
            ("layout_png_download_url", "layout_board.png"),
        ):
            src = _download_url_to_path(export_result.get(key))
            if src is not None and src.is_file():
                zf.write(src, arcname)
        if pattern_result is not None:
            uv_src = _download_url_to_path(pattern_result.get("pattern_uv_png_url"))
            if uv_src is not None and uv_src.is_file():
                zf.write(uv_src, "pattern_uv.png")
            for i, cut_url in enumerate(pattern_result.get("pattern_cut_png_urls") or ()):
                cut_src = _download_url_to_path(cut_url)
                if cut_src is not None and cut_src.is_file():
                    zf.write(cut_src, f"pattern_cut_{i}.png")
        if synthesize_result is not None:
            layout_src = _download_url_to_path(synthesize_result.get("layout_png_url"))
            if layout_src is not None and layout_src.is_file():
                zf.write(layout_src, "layout.png")
            patch_src = _download_url_to_path(synthesize_result.get("patch_layout_png_url"))
            if patch_src is not None and patch_src.is_file():
                zf.write(patch_src, "patch_layout.png")
        zf.writestr("project_record.json", json.dumps(record, indent=2))
    return {"download_url": f"/api/files/{job.id}/{zip_path.name}"}


@router.post("/package")
def export_package(body: ExportPackageRequest, request: Request):
    export_job = request.app.state.jobs.get(body.export_job_id)
    if export_job is None:
        raise HTTPException(404, f"job {body.export_job_id!r} not found")
    if export_job.status != JobStatus.DONE:
        raise HTTPException(409, f"export job {body.export_job_id!r} is {export_job.status.value}, not done")

    pattern_result = None
    if body.pattern_job_id is not None:
        pattern_job = request.app.state.jobs.get(body.pattern_job_id)
        if pattern_job is None:
            raise HTTPException(404, f"job {body.pattern_job_id!r} not found")
        if pattern_job.status != JobStatus.DONE:
            raise HTTPException(409, f"pattern job {body.pattern_job_id!r} is {pattern_job.status.value}, not done")
        pattern_result = pattern_job.result

    synthesize_result = None
    if body.synthesize_job_id is not None:
        synthesize_job = request.app.state.jobs.get(body.synthesize_job_id)
        if synthesize_job is None:
            raise HTTPException(404, f"job {body.synthesize_job_id!r} not found")
        if synthesize_job.status != JobStatus.DONE:
            raise HTTPException(
                409, f"synthesize job {body.synthesize_job_id!r} is {synthesize_job.status.value}, not done"
            )
        synthesize_result = synthesize_job.result

    job = request.app.state.jobs.submit(
        "export_package", _run_package_job, export_job.result, pattern_result, synthesize_result
    )
    return job.to_dict()
