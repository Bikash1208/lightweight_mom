"""Phase 6 orchestration: for one `ValidationCase`, write the validation-case folder described in
the Phase 6 brief. Reuses the existing pipeline unchanged:

    config (Phase 1) -> runner.run_sweep (Phase 5, in memory, one cell x this case's frequency
    grid) -> Database (Phase 5) -> mom_results.json / mom_results.csv (Database.to_json/to_csv,
    byte for byte the Phase 5 format) -> plots (matplotlib, from Database.series /
    phase_series_deg) -> README.md (a template filled from the case, its config and its results).

cad_export.py builds the STEP file, independently, from the same `SimulationConfig.unit_cell`.
Nothing here touches Z, V, the Green's function or `compute_s_parameters`.

Folder layout (one per case, see VALIDATION_EXPORT.md):

    <case_id>/
      unit_cell.step
      geometry.json
      solver_config.json
      mom_results.csv
      mom_results.json
      plots/{s11,s21,co_polarized,cross_polarized}.png
      README.md
"""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import __version__, units
from .cad_export import DEFAULT_METAL_EXPORT_THICKNESS_M, export_unit_cell_step
from .database import CONVENTIONS, Database
from .runner import run_sweep
from .sweep import SweepDefinition
from .validation_cases import ValidationCase

_COMPONENTS_XX_YY = ("S11_xx", "S21_xx", "S11_yy", "S21_yy")
_COMPONENTS_CROSS = ("S11_yx", "S21_yx", "S11_xy", "S21_xy")
_DB_FLOOR = -300.0  # dB shown for a magnitude at or below floating-point noise


@dataclass(frozen=True, eq=False)
class CaseExportResult:
    case: ValidationCase
    out_dir: Path
    step_path: Path
    step_result: object  # cad_export.StepExportResult
    database: Database
    geometry_json_path: Path
    solver_config_json_path: Path
    mom_results_json_path: Path
    mom_results_csv_path: Path
    plot_paths: dict  # name -> Path
    readme_path: Path


# ------------------------------------------------------------------------------- solving


def _solve_case(case: ValidationCase):
    """The case's Database: run_sweep with no geometry axes, so exactly one cell, over the
    case's own frequency grid. Phase 5 code, unmodified; in memory (out_dir=None)."""
    definition = SweepDefinition(frequency=case.frequency, axes=())
    return run_sweep(case.config, out_dir=None, definition=definition).database


# --------------------------------------------------------------------------------- JSON


def geometry_dict(case: ValidationCase, metal_export_thickness_m):
    uc = case.config.unit_cell
    geo, lat, sub = uc.geometry, uc.lattice, uc.substrate
    return {
        "case_id": case.case_id,
        "category": case.category,
        "description": case.description,
        "code_version": __version__,
        "lattice": {"period_x_m": lat.period_x, "period_y_m": lat.period_y},
        "substrate": {
            "name": sub.name,
            "eps_r": sub.eps_r,
            "tan_delta": sub.tan_delta,
            "thickness_m": sub.thickness,
        },
        "metal": {"model": uc.metal.model},
        "geometry": {
            "type": geo.kind,
            "polarity": geo.polarity,
            "length_m": geo.length,
            "width_m": geo.width,
            "center_x_m": geo.center_x,
            "center_y_m": geo.center_y,
            "rotation_deg": units.to_deg(geo.rotation),
            "metal_area_m2": uc.metal_area,
            "metal_fill_fraction": uc.metal_fill_fraction,
        },
        "cad_export": {
            "metal_export_thickness_m": metal_export_thickness_m,
            "note": (
                "CAD solid convenience only (a zero-thickness solid cannot be exported); the MoM "
                "model is a zero-thickness PEC surface current at z=0 and does not use this value. "
                "See cad_export.py and VALIDATION_EXPORT.md Sec. 4."
            ),
        },
    }


def solver_config_dict(case: ValidationCase, database: Database):
    cfg, rp = case.config, case.config.reference_planes
    mesh = cfg.mesh
    return {
        "mesh": {"cells_x": mesh.cells_x, "cells_y": mesh.cells_y, "floquet_order": mesh.floquet_order},
        "reference_planes": {
            "top_air_thickness_m": rp.top_air_thickness,
            "bottom_air_thickness_m": rp.bottom_air_thickness,
            "z_top_m": rp.z_top,
            "z_bottom_m": rp.z_bottom,
        },
        "frequency_grid": database.metadata["frequency_grid"],
        "excitation": {
            "description": (
                "Unit-amplitude (1 V/m) normal-incidence plane wave, phase referenced to the "
                "pattern plane z=0; both x- and y-polarized incidence are solved for every "
                "frequency (both give a column of the 4x4 S-matrix, see S_PARAMETERS.md)."
            ),
            "incidence": "normal (theta=0)",
            "polarizations_solved": ["x", "y"],
        },
        "conventions": CONVENTIONS,
    }


# -------------------------------------------------------------------------------- plots


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")  # headless: no display, no GUI backend dependency
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for the Phase 6 plots (src/validation_export.py); install it "
            "with 'pip install -r requirements-cad.txt'."
        ) from e
    return plt


def _mag_phase(database, name):
    freqs, s = database.series(0, name)
    _, phase_deg = database.phase_series_deg(0, name)
    return freqs / 1e9, np.abs(s), phase_deg


def _db(mag):
    return 20.0 * np.log10(np.maximum(mag, 10 ** (_DB_FLOOR / 20.0)))


def _plot_family(plt, database, path, family, title):
    """One figure, two subplots (magnitude, phase) for <family>_xx and <family>_yy."""
    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(6.4, 6.0), sharex=True)
    for comp, style in ((f"{family}_xx", "-"), (f"{family}_yy", "--")):
        f_ghz, mag, phase_deg = _mag_phase(database, comp)
        ax_mag.plot(f_ghz, mag, style, marker="o", markersize=3, label=comp)
        ax_phase.plot(f_ghz, phase_deg, style, marker="o", markersize=3, label=comp)
    ax_mag.set_ylabel(f"|{family}| (linear)")
    ax_mag.set_title(title)
    ax_mag.grid(True, alpha=0.3)
    ax_mag.legend()
    ax_phase.set_ylabel(f"phase({family}) [deg]")
    ax_phase.set_xlabel("frequency [GHz]")
    ax_phase.set_ylim(-180, 180)
    ax_phase.grid(True, alpha=0.3)
    ax_phase.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _plot_overview(plt, database, path, components, title):
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for comp in components:
        f_ghz, mag, _ = _mag_phase(database, comp)
        ax.plot(f_ghz, _db(mag), marker="o", markersize=3, label=comp)
    ax.set_xlabel("frequency [GHz]")
    ax.set_ylabel("magnitude [dB]")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def generate_plots(database, plots_dir):
    """Write the four PNGs into ``plots_dir`` (created if needed); return {name: path}."""
    plt = _require_matplotlib()
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "s11": plots_dir / "s11.png",
        "s21": plots_dir / "s21.png",
        "co_polarized": plots_dir / "co_polarized.png",
        "cross_polarized": plots_dir / "cross_polarized.png",
    }
    _plot_family(plt, database, paths["s11"], "S11", "S11 (reflection): magnitude and phase vs frequency")
    _plot_family(plt, database, paths["s21"], "S21", "S21 (transmission): magnitude and phase vs frequency")
    _plot_overview(
        plt, database, paths["co_polarized"], _COMPONENTS_XX_YY, "Co-polarized response (xx, yy) vs frequency"
    )
    _plot_overview(
        plt, database, paths["cross_polarized"], _COMPONENTS_CROSS, "Cross-polarized response (yx, xy) vs frequency"
    )
    return paths


# ------------------------------------------------------------------------------- README


def render_readme(case: ValidationCase, step_result, database: Database, plot_paths):
    uc = case.config.unit_cell
    geo, lat, sub = uc.geometry, uc.lattice, uc.substrate
    rp = case.config.reference_planes
    mesh = case.config.mesh
    grid = database.metadata["frequency_grid"]
    t_metal = step_result.metal_export_thickness_m
    n_unknowns = database.records[0].diagnostics["n_unknowns"]
    plot_lines = "\n".join(f"- `plots/{p.name}`" for p in plot_paths.values())

    return f"""# {case.case_id}

**Category:** {case.category}. {case.description}

Generated by `src/validation_export.py` (Phase 6) from the same `SimulationConfig` used to run
the MoM solver (`src/validation_cases.py`); nothing in this folder was hand-edited. Time
convention throughout this project, including every number below: **`exp(+j*omega*t)`** (a
phasor `F` represents `Re{{F * exp(+j*omega*t)}}`; a delay is a *negative* phase). Full
conventions: `CONVENTIONS.md`; S-parameter/port definitions: `S_PARAMETERS.md`.

## Frequency range
{units.to_ghz(grid["start_hz"]):.4g} GHz to {units.to_ghz(grid["stop_hz"]):.4g} GHz, step
{units.to_ghz(grid["step_hz"]):.4g} GHz ({grid["count"]} points): `{", ".join(f"{units.to_ghz(f):.4g}" for f in grid["frequencies_hz"])}` GHz.

## Unit-cell period
`period_x` = {units.to_mm(lat.period_x):.6g} mm, `period_y` = {units.to_mm(lat.period_y):.6g} mm
(square lattice, rectangular grid).

## Substrate
`{sub.name}`: `eps_r` = {sub.eps_r:g}, `tan_delta` = {sub.tan_delta:g}, `thickness` =
{units.to_mm(sub.thickness):.6g} mm. Non-magnetic (`mu_r` = 1), constant (non-dispersive) over
frequency (`materials.py`). Complex permittivity `eps_c = eps_r*(1 - j*tan_delta)`.

## Metal geometry
`geometry.type` = `{geo.kind}` ({geo.polarity} polarity{"" if geo.polarity == "metal" else ": the metal is the complement of this feature (a solid plate with the feature cut out)"}),
`length` = {units.to_mm(geo.length):.6g} mm, `width` = {units.to_mm(geo.width):.6g} mm,
`center` = ({units.to_mm(geo.center_x):.6g}, {units.to_mm(geo.center_y):.6g}) mm,
`rotation` = {units.to_deg(geo.rotation):.6g} deg (counter-clockwise about the feature centre).
Metal fill fraction of the cell: {uc.metal_fill_fraction:.4f}. Metal model: `{uc.metal.model}`
(idealized PEC, zero thickness in the MoM model — see "CAD export" below).

## Mesh parameters
Rooftop mesh: `{mesh.cells_x}` x `{mesh.cells_y}` cells ({units.to_mm(lat.period_x / mesh.cells_x):.4g} x
{units.to_mm(lat.period_y / mesh.cells_y):.4g} mm each); `{n_unknowns}` rooftop unknowns for this
geometry (`basis_functions.py`).

## Floquet truncation
`floquet_order` = `{mesh.floquet_order}` (Floquet harmonics `|m|, |n| <= {mesh.floquet_order}` in
the impedance-matrix assembly; rank rule `2*floquet_order + 1 >= max(cells_x, cells_y)`,
`MOM_FORMULATION.md`).

## Excitation convention
Unit-amplitude (`|E_inc| = 1 V/m`) normal-incidence (`theta = 0`) plane wave, phase referenced to
the pattern plane `z = 0`; **both** x- and y-polarized incidence are solved (giving the two
columns of the 4x4 S-matrix), regardless of any single `polarization` value in the case config
(`excitation.py`, `sparameters.py`).

## Reference-plane locations
Port 1 (`top`) at `z_top` = {units.to_mm(rp.z_top):.6g} mm (`{units.to_mm(rp.top_air_thickness):.6g}` mm
above the pattern plane); port 2 (`bottom`) at `z_bottom` = {units.to_mm(rp.z_bottom):.6g} mm
(`{units.to_mm(rp.bottom_air_thickness):.6g}` mm below the substrate's bottom face). Both planes are
in air; only the specular (0,0) Floquet order is a port mode (`S_PARAMETERS.md`).

## exp(+jwt) convention
This project uses **`exp(+j*omega*t)`** throughout: `exp(-j*k*z)` travels toward `+z`,
`exp(+j*k*z)` travels toward `-z`, and a lossy medium has `Im(eps) <= 0`, `Im(k) <= 0`
(`CONVENTIONS.md`). If you compare against a tool whose default is `exp(-j*omega*t)`, conjugate
every S-parameter (magnitude unchanged, phase and the sign of any loss term flipped) before
comparing.

## Files in this folder
- `unit_cell.step` — substrate + metal solids (see "CAD export" below).
- `geometry.json` — every geometry/material parameter above, machine-readable, exactly as given
  to the solver (`config_from_dict`), plus the `cad_export.metal_export_thickness_m` used only
  for the STEP solid.
- `solver_config.json` — mesh, Floquet order, reference planes, frequency grid, excitation
  convention and the full sign/unit convention reference (`database.CONVENTIONS`).
- `mom_results.json` / `mom_results.csv` — the Phase 5 unit-cell database format (`DATABASE.md`)
  for this one cell: complex S11/S12/S21/S22 (linear x/y **and** derived circular R/L bases),
  magnitude and phase (principal value, degrees), plus per-record diagnostics
  (`incoming_residual`, `reciprocity_error`, `passivity_min_eig`, `cells_per_wavelength`,
  `specular_only`, `checks_passed`, `n_unknowns`).
{plot_lines}

## CAD export
`unit_cell.step` contains two named solids, built directly from the same geometry parameters
listed above (`cad_export.py`; millimetres, the STEP default unit — the solver works internally
in metres and this file is the one place besides `units.py` that converts):

- **`Substrate`**: a box of `period_x` x `period_y` x `substrate.thickness`, `z` from
  `-{units.to_mm(sub.thickness):.6g}` mm to `0` mm (matching the solver's stack: substrate between
  the pattern plane `z=0` and the substrate's bottom face at `z=-thickness`).
- **`Metal`**: a thin solid slab, `z` from `0` to `{units.to_mm(t_metal):.6g}` mm
  (`{units.to_mm(t_metal) * 1000:.0f}` um, a nominal 1 oz copper thickness), flush with the
  substrate's top face. **This thickness has no effect on any solved quantity in this project** —
  the MoM model is a strictly zero-thickness PEC surface current at `z=0`; the thickness exists
  only because a STEP solid cannot have zero volume. Do not read S-parameter sensitivity to metal
  thickness into this export; that has not been modelled or tested.

## Reproducing this model in HFSS (manual — not run or checked by this project)
1. Set the modeler length unit to **millimetres** (the STEP file's native unit) before importing,
   then `Import` -> `unit_cell.step`. Confirm two objects, `Substrate` and `Metal`.
2. Assign `Substrate` a dielectric material with relative permittivity `{sub.eps_r:g}` and loss
   tangent `{sub.tan_delta:g}` (create a new material if `{sub.name}` is not in the library).
3. Assign the `Metal` object a **Perfect E** boundary (Assign Boundary -> Perfect E on its
   surface). The small `{units.to_mm(t_metal) * 1000:.0f}` um export thickness is not part of the
   MoM model (see "CAD export" above) — a Perfect E boundary makes the volume electromagnetically
   irrelevant, matching the zero-thickness assumption; do not assign a finite-conductivity copper
   material unless you are deliberately testing thickness/loss sensitivity outside this project.
4. Draw the periodic unit cell's air region as a box of footprint `{units.to_mm(lat.period_x):.6g}`
   x `{units.to_mm(lat.period_y):.6g}` mm, from `z = {units.to_mm(rp.z_bottom):.6g}` mm to
   `z = {units.to_mm(rp.z_top):.6g}` mm (matching `z_top`/`z_bottom` above), centred on the
   `Substrate`/`Metal` solids (both already centred at `x=y=0`).
5. Assign **periodic (lattice-pair / master-slave) boundaries** on the four side faces of that
   box, pairing `x = -{units.to_mm(lat.period_x) / 2:.6g}` mm with `x = +{units.to_mm(lat.period_x) / 2:.6g}`
   mm and `y = -{units.to_mm(lat.period_y) / 2:.6g}` mm with `y = +{units.to_mm(lat.period_y) / 2:.6g}` mm,
   zero phase delay (normal incidence).
6. Assign **Floquet ports** on the top face (`z = {units.to_mm(rp.z_top):.6g}` mm) and bottom face
   (`z = {units.to_mm(rp.z_bottom):.6g}` mm), each with two modes (TE00/TM00, i.e. the linear x and
   y polarizations of the (0,0) Floquet order) so the resulting S-matrix is the same 4x4
   [[S11,S12],[S21,S22]] over (port1 x, port1 y, port2 x, port2 y) as `mom_results.json`.
7. Set an interpolating frequency sweep from {units.to_ghz(grid["start_hz"]):.4g} GHz to
   {units.to_ghz(grid["stop_hz"]):.4g} GHz (at least the {grid["count"]} points listed above).
8. Solve, then compare magnitude directly; compare phase after checking your HFSS version's time
   convention against "exp(+jwt) convention" above (conjugate if it differs) and, if your Floquet
   port reference planes are not exactly `z_top`/`z_bottom` above, re-reference this project's
   values with `SParameters.at_planes` (or HFSS's own port de-embedding) before comparing phase.

This project does **not** perform this comparison automatically and does **not** tune its MoM
results to match HFSS; any difference you find is data for you to interpret; see
`VALIDATION_EXPORT.md` for what has and has not been checked.
"""


# --------------------------------------------------------------------------- orchestration


def export_case(case: ValidationCase, out_dir, metal_export_thickness_m=DEFAULT_METAL_EXPORT_THICKNESS_M):
    """Write the full validation-case folder for ``case`` into ``out_dir`` (created if needed)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = out_dir / "plots"

    database = _solve_case(case)

    step_path = out_dir / "unit_cell.step"
    step_result = export_unit_cell_step(case.config.unit_cell, step_path, metal_export_thickness_m)

    geometry_json_path = out_dir / "geometry.json"
    geometry_json_path.write_text(json.dumps(geometry_dict(case, metal_export_thickness_m), indent=1) + "\n", encoding="utf-8", newline="\n")

    solver_config_json_path = out_dir / "solver_config.json"
    solver_config_json_path.write_text(json.dumps(solver_config_dict(case, database), indent=1) + "\n", encoding="utf-8", newline="\n")

    mom_results_json_path = out_dir / "mom_results.json"
    database.to_json(mom_results_json_path)
    mom_results_csv_path = out_dir / "mom_results.csv"
    database.to_csv(mom_results_csv_path)

    plot_paths = generate_plots(database, plots_dir)

    readme_path = out_dir / "README.md"
    readme_path.write_text(render_readme(case, step_result, database, plot_paths), encoding="utf-8", newline="\n")

    return CaseExportResult(
        case,
        out_dir,
        step_path,
        step_result,
        database,
        geometry_json_path,
        solver_config_json_path,
        mom_results_json_path,
        mom_results_csv_path,
        plot_paths,
        readme_path,
    )


def export_all_cases(base_dir, cases=None, metal_export_thickness_m=DEFAULT_METAL_EXPORT_THICKNESS_M, progress=None):
    """``export_case`` for every case (default: `validation_cases.CASES`) into
    ``base_dir/<case_id>/``. Returns {case_id: CaseExportResult}."""
    from .validation_cases import CASES as _ALL_CASES

    cases = cases if cases is not None else _ALL_CASES
    base_dir = Path(base_dir)
    results = {}
    for i, case in enumerate(cases):
        results[case.case_id] = export_case(case, base_dir / case.case_id, metal_export_thickness_m)
        if progress is not None:
            progress({"done": i + 1, "total": len(cases), "case_id": case.case_id})
    return results
