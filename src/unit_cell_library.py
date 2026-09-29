"""Phase 7: `unit_cell_library.json`/`.csv` (req. 5, 20), the manifest (req. 19),
`LIBRARY_REPORT.md` (req. 31), and the read-only query/inversion API (req. 24-26).

Every S-parameter value written or served by this module is read from an already-solved
`database.Database` (`library_builder.py`); nothing here re-solves, interpolates or extrapolates.
"""

import csv
import hashlib
import json
import platform
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from . import __version__ as SOLVER_VERSION
from .database import CONVENTIONS
from .library_schema import (
    LIBRARY_CSV_COLUMNS,
    PHASE_COLLISION_BIN_DEG,
    PHASE_USABLE_MIN_MAG_DB,
    RESONANCE_MIN_PROMINENCE,
    LibraryConfig,
)
from .validation_export import _db

LIBRARY_FORMAT_NAME = "lightweight_mom.unit_cell_library"
LIBRARY_FORMAT_VERSION = 1

_S_FAMILIES = ("S11", "S21", "S12", "S22")
POLARIZATION_CONVENTION = f"{CONVENTIONS['linear_basis']}; {CONVENTIONS['s_matrix']}"
TIME_CONVENTION = CONVENTIONS["time_dependence"]


# ------------------------------------------------------------------------------ JSON / CSV export


def _s_entry(record, name):
    z = record.component(name)
    mag = abs(z)
    return {
        "real": float(z.real),
        "imag": float(z.imag),
        "mag": float(mag),
        "mag_db": float(_db(np.array([mag]))[0]),
        "phase_deg": float(record.phase_deg(name)),
    }


def _cell_doc(cell, database):
    l_mm = round(cell.parameters["length_m"] * 1e3, 9)
    w_mm = round(cell.parameters["width_m"] * 1e3, 9)
    geometry_status = "VALID" if cell.status == "ok" else "INVALID"
    if cell.status != "ok":
        return {
            "cell_id": cell.cell_id,
            "L_mm": l_mm,
            "W_mm": w_mm,
            "geometry_status": geometry_status,
            "solver_status": "FAILED",
            "invalid_reason": cell.reason,
            "frequency_GHz": [],
            "S11": [],
            "S21": [],
            "S12": [],
            "S22": [],
        }
    records = sorted(database.records_of(cell.index), key=lambda r: r.frequency_index)
    ctx = database.metadata["context"]
    doc = {
        "cell_id": cell.cell_id,
        "L_mm": l_mm,
        "W_mm": w_mm,
        "period_x_mm": ctx["lattice"]["period_x_m"] * 1e3,
        "period_y_mm": ctx["lattice"]["period_y_m"] * 1e3,
        "substrate_thickness_mm": ctx["substrate"]["thickness_m"] * 1e3,
        "epsilon_r": ctx["substrate"]["eps_r"],
        "tan_delta": ctx["substrate"]["tan_delta"],
        "polarization_convention": POLARIZATION_CONVENTION,
        "time_convention": TIME_CONVENTION,
        "geometry_status": geometry_status,
        "solver_status": "OK",
        "frequency_GHz": [r.frequency / 1e9 for r in records],
    }
    for fam in _S_FAMILIES:
        doc[fam] = [_s_entry(r, f"{fam}_xx") for r in records]
    return doc


def write_library_json(database, path):
    """`unit_cell_library.json` (req. 18): one entry per geometry (req. 5's fields, nested by
    frequency rather than flattened -- semantically the same data as the flat CSV)."""
    doc = {
        "format": LIBRARY_FORMAT_NAME,
        "version": LIBRARY_FORMAT_VERSION,
        "conventions": CONVENTIONS,
        "cells": [_cell_doc(c, database) for c in database.cells],
    }
    path = Path(path)
    path.write_text(json.dumps(doc, indent=1, allow_nan=False) + "\n")
    return path


def write_library_csv(database, path):
    """`unit_cell_library.csv` (req. 18, 20): exactly `library_schema.LIBRARY_CSV_COLUMNS`, one
    row per (cell, frequency), successful cells only (a failed cell has no S-parameters to
    report; it is recorded in `unit_cell_library.json` and the manifest's failure count, req. 4,
    never as a row of invented numbers here)."""
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(LIBRARY_CSV_COLUMNS)
        for cell in database.cells:
            if cell.status != "ok":
                continue
            l_mm = round(cell.parameters["length_m"] * 1e3, 9)
            w_mm = round(cell.parameters["width_m"] * 1e3, 9)
            for r in sorted(database.records_of(cell.index), key=lambda r: r.frequency_index):
                row = [cell.cell_id, repr(l_mm), repr(w_mm), repr(r.frequency / 1e9)]
                for fam in _S_FAMILIES:
                    e = _s_entry(r, f"{fam}_xx")
                    row += [repr(e["real"]), repr(e["imag"]), repr(e["mag"]), repr(e["mag_db"]), repr(e["phase_deg"])]
                writer.writerow(row)
    return path


def write_phase_lookup_csv(table, path):
    """`phase_lookup_12p5GHz.csv` (req. 18): req. 10's table, one row per successful geometry."""
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["L_mm", "W_mm", "phase_deg", "magnitude_db", "S11_real", "S11_imag"])
        for row in table:
            writer.writerow(
                [
                    repr(row["L_mm"]),
                    repr(row["W_mm"]),
                    repr(row["phase_deg"]),
                    repr(row["magnitude_db"]),
                    repr(row["S11_real"]),
                    repr(row["S11_imag"]),
                ]
            )
    return path


# --------------------------------------------------------------------------------------- manifest


def _configuration_hash(config: LibraryConfig):
    payload = json.dumps(
        {
            "period_m": config.period_m,
            "substrate": config.substrate,
            "f_min_hz": config.f_min_hz,
            "f_max_hz": config.f_max_hz,
            "f_step_hz": config.f_step_hz,
            "f_design_hz": config.f_design_hz,
            "l_values_mm": config.l_values_mm,
            "w_values_mm": config.w_values_mm,
            "mesh": config.mesh,
            "reference_planes": config.reference_planes,
        },
        sort_keys=True,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def build_manifest(config, build, phase_source, deterministic_timestamp=None):
    """Req. 19: every field req. 19 lists, plus `configuration_hash` (sha256 of every physical
    input) so "exactly which solver configuration generated the library" (req. 19) is checkable
    without re-reading every config field by hand. ``deterministic_timestamp``: tests pass a
    fixed value so manifest-generation tests do not depend on wall-clock time; the real build
    (examples/build_unit_cell_library.py) leaves it `None` and gets the real UTC time -- the one
    intentionally non-reproducible field, same as every other phase's file headers."""
    timestamp = deterministic_timestamp or datetime.now(UTC).isoformat()
    return {
        "library_version": 1,
        "creation_timestamp": timestamp,
        "solver_version": SOLVER_VERSION,
        "configuration_hash": _configuration_hash(config),
        "number_of_cells": config.n_geometries,
        "number_of_successful_cells": build.n_ok_cells,
        "number_of_failed_cells": build.n_invalid_cells,
        "frequency_count": int(round((config.f_max_hz - config.f_min_hz) / config.f_step_hz)) + 1,
        "frequency_min_GHz": config.f_min_hz / 1e9,
        "frequency_max_GHz": config.f_max_hz / 1e9,
        "design_frequency_GHz": config.f_design_hz / 1e9,
        "period_x_mm": config.period_m * 1e3,
        "period_y_mm": config.period_m * 1e3,
        "substrate_thickness_mm": config.substrate["thickness_m"] * 1e3,
        "phase_source": phase_source.component,
        "phase_source_status": phase_source.status,
        "time_convention": TIME_CONVENTION,
        "polarization_convention": POLARIZATION_CONVENTION,
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
    }


def write_manifest(manifest, path):
    path = Path(path)
    path.write_text(json.dumps(manifest, indent=1, allow_nan=False) + "\n")
    return path


# --------------------------------------------------------------------------------- query API


@dataclass(frozen=True)
class PhaseCandidate:
    cell_id: str
    l_mm: float
    w_mm: float
    phase_deg: float
    phase_error_deg: float
    magnitude_db: float


class UnitCellLibrary:
    """Read-only query interface over an exported `unit_cell_library.json` (req. 24). Every
    lookup is an *exact* match against the library's own stored geometries/frequencies -- no
    interpolation, no extrapolation (req. 26)."""

    def __init__(self, doc):
        self._doc = doc
        self._cells = {(c["L_mm"], c["W_mm"]): c for c in doc["cells"]}

    @classmethod
    def from_json(cls, path):
        return cls(json.loads(Path(path).read_text()))

    def get_cell(self, L_mm, W_mm):
        """The raw cell document for (L_mm, W_mm), or `None` if that exact geometry is not in
        the library."""
        return self._cells.get((round(L_mm, 9), round(W_mm, 9)))

    def response(self, L_mm, W_mm, frequency_GHz, family="S11", rtol=1e-6):
        """The stored S-parameter entry ({real, imag, mag, mag_db, phase_deg}) of ``family``
        (default S11) for an exact (L_mm, W_mm) at an exact frequency already in the library.
        Raises KeyError/ValueError rather than interpolating if either is not present."""
        cell = self.get_cell(L_mm, W_mm)
        if cell is None:
            raise KeyError(f"no cell at L={L_mm} mm, W={W_mm} mm in this library")
        if cell["solver_status"] != "OK":
            raise ValueError(
                f"cell {cell['cell_id']} has solver_status={cell['solver_status']!r}, no response available"
            )
        freqs = cell["frequency_GHz"]
        matches = [i for i, f in enumerate(freqs) if abs(f - frequency_GHz) <= rtol * max(abs(frequency_GHz), 1.0)]
        if not matches:
            raise ValueError(
                f"frequency {frequency_GHz} GHz is not one of this cell's {len(freqs)} stored "
                f"frequencies ({freqs[0]:g}..{freqs[-1]:g} GHz); no interpolation is performed"
            )
        return cell[family][matches[0]]

    def find_phase_candidates(self, target_phase_deg, frequency_GHz, tolerance_deg=None, family="S11"):
        """Every cell whose ``family`` phase at ``frequency_GHz`` (an exact stored frequency,
        req. 26) is within ``tolerance_deg`` (default `library_schema.DEFAULT_PHASE_TOLERANCE_DEG`)
        of ``target_phase_deg`` -- **all** matches, never a single chosen winner (req. 24)."""
        from .library_schema import DEFAULT_PHASE_TOLERANCE_DEG
        from .phase import wrapped_difference_deg

        tolerance_deg = DEFAULT_PHASE_TOLERANCE_DEG if tolerance_deg is None else tolerance_deg
        candidates = []
        for cell in self._doc["cells"]:
            if cell["solver_status"] != "OK":
                continue
            freqs = cell["frequency_GHz"]
            matches = [i for i, f in enumerate(freqs) if abs(f - frequency_GHz) <= 1e-6 * max(abs(frequency_GHz), 1.0)]
            if not matches:
                continue
            entry = cell[family][matches[0]]
            error = float(wrapped_difference_deg(entry["phase_deg"], target_phase_deg))
            if abs(error) <= tolerance_deg:
                candidates.append(
                    PhaseCandidate(
                        cell["cell_id"],
                        cell["L_mm"],
                        cell["W_mm"],
                        entry["phase_deg"],
                        error,
                        entry["mag_db"],
                    )
                )
        return sorted(candidates, key=lambda c: abs(c.phase_error_deg))

    def invert_phase(self, target_phase_deg, frequency_GHz, family="S11"):
        """Req. 25: the single closest-phase candidate at ``frequency_GHz`` (an exact library
        frequency), or the string `"NO_VALID_CANDIDATE"` if none is within
        `library_schema.DEFAULT_PHASE_TOLERANCE_DEG` -- never an extrapolated/interpolated
        geometry (req. 26)."""
        candidates = self.find_phase_candidates(target_phase_deg, frequency_GHz, family=family)
        if not candidates:
            return "NO_VALID_CANDIDATE"
        best = candidates[0]
        return {
            "target_phase_deg": target_phase_deg,
            "selected_phase_deg": best.phase_deg,
            "phase_error_deg": best.phase_error_deg,
            "L_mm": best.l_mm,
            "W_mm": best.w_mm,
            "response_magnitude_db": best.magnitude_db,
            "cell_id": best.cell_id,
        }


def phase_lookup_filename(f_design_hz):
    ghz = f_design_hz / 1e9
    return f"phase_lookup_{ghz:g}GHz".replace(".", "p") + ".csv"


# ------------------------------------------------------------------------------------- report


def render_library_report_md(
    config, build, coverage, collisions, classifications, geometry_exports, phase_source, manifest, consistency
):
    n_resonant = sum(1 for c in classifications if c.resonance is not None and c.resonance.is_resonant)
    n_non_resonant = sum(1 for c in classifications if c.resonance is not None and not c.resonance.is_resonant)
    n_phase_usable = sum(1 for c in classifications if c.phase_status == "PHASE_USABLE")
    n_phase_unusable = sum(1 for c in classifications if c.phase_status == "PHASE_UNUSABLE")
    n_grating_risk = sum(1 for c in classifications if c.grating is not None and c.grating["is_grating_risk"])
    failed = [c for c in classifications if c.solver_status != "OK" or c.geometry_status != "VALID"]
    onset_hz = next((c.grating["grating_onset_hz"] for c in classifications if c.grating is not None), None)

    lines = [
        "# LIBRARY_REPORT.md",
        "",
        "Production rectangular-PEC-patch unit-cell library (Phase 7), generated from the "
        "existing, unmodified lightweight periodic MoM solver (Phases 1-6A). No EM formulation, "
        "Green's function, MoM matrix, basis functions or S-parameter code was changed to build "
        "this library.",
        "",
        "## 1. Physical configuration",
        "",
        f"- Period: {config.period_m * 1e3:g} x {config.period_m * 1e3:g} mm (square lattice).",
        f"- Substrate: `{config.substrate['name']}` (authoritative source: `config.yaml`, this "
        f"project's own root config), eps_r = {config.substrate['eps_r']:g}, tan_delta = "
        f"{config.substrate['tan_delta']:g}, thickness = {config.substrate['thickness_m'] * 1e3:g} mm.",
        "- Metal: PEC (`materials.py`, Phase 1's only supported model).",
        f"- Patch sweep: L, W in {{{config.l_values_mm[0]:g}, ..., {config.l_values_mm[-1]:g}}} mm, "
        f"step {config.l_values_mm[1] - config.l_values_mm[0]:g} mm "
        f"({len(config.l_values_mm)} x {len(config.w_values_mm)} = {config.n_geometries} geometries).",
        f"- Mesh: {config.mesh['cells_x']}x{config.mesh['cells_y']} cells, Floquet order "
        f'{config.mesh["floquet_order"]} (accuracy trade-off documented in Sec. 16, "Limitations").',
        "",
        "## 2. Number of geometries",
        "",
        f"{config.n_geometries} enumerated ({len(config.l_values_mm)}x{len(config.w_values_mm)}); "
        f"{build.n_ok_cells} solved OK, "
        f"{build.n_invalid_cells} invalid/skipped (a priori geometry check, before any solve). "
        f"No geometry was silently dropped or substituted (req. 4/28).",
        "",
        "## 3. Frequency range",
        "",
        f"{config.f_min_hz / 1e9:g} GHz to {config.f_max_hz / 1e9:g} GHz, step {config.f_step_hz / 1e6:g} MHz, "
        f"{manifest['frequency_count']} points inclusive; design frequency "
        f"{config.f_design_hz / 1e9:g} GHz is an exact grid point (verified, never interpolated).",
        "",
        "## 4. Solver convention",
        "",
        f"Time convention: `{TIME_CONVENTION}` (CONVENTIONS.md). S-parameters: "
        f"{CONVENTIONS['s_matrix']}. Phase: {CONVENTIONS['phase']}. No phase unwrapping is "
        "performed anywhere in this library (req. 7, 9).",
        "",
        "## 5. Polarization convention",
        "",
        POLARIZATION_CONVENTION,
        "",
        "## 6. Database structure",
        "",
        "`database/unit_cell_library.json` (one entry per geometry, nested by frequency), "
        "`database/unit_cell_library.csv` (flat, one row per (cell, frequency), "
        f"{len(LIBRARY_CSV_COLUMNS)} columns), `database/{phase_lookup_filename(config.f_design_hz)}` "
        "(req. 10's design-frequency lookup table). Both JSON and CSV re-export the same "
        "underlying `database.Database` (Phase 5, unmodified) this build produced; internal "
        "consistency between them: "
        f"{'**consistent**' if consistency.ok else f'**{len(consistency.mismatches)} mismatch(es)**'} "
        "(`reports/LIBRARY_VALIDATION.md`).",
        "",
        "**HFSS comparison compatibility (req. 23):** this library uses exactly the conventions "
        "Phase 6A's `HFSS_COMPARISON.md` already documents per case, not a second incompatible "
        "one: `S11`/`S21`/`S12`/`S22` here are the same `S_ij[received, incident]`, index 0=x, "
        "1=y, `exp(+j*omega*t)` quantities as `S_PARAMETERS.md`; `MoM S11xx <-> HFSS S(Port1:x, "
        "Port1:x)`, `MoM S21xx <-> HFSS S(Port2:x, Port1:x)` (Floquet port modes x/y = the "
        "(0,0)-order linear polarizations, `HFSS_COMPARISON.md` Sec. 2). If your HFSS export "
        "uses `exp(-j*omega*t)`, conjugate every HFSS S-parameter before comparing.",
        "",
        "## 7. Geometry generation",
        "",
        f"{len(geometry_exports)} STEP files written to `geometry/` (`UC_L<L>_W<W>.step`), one "
        "per successfully-solved geometry, built directly from the same `UnitCell` the solver "
        "used (`cad_export.py`, Phase 6, unmodified). Every export passed its independent bounds "
        "check (req. 17: 0 < L < period_x, 0 < W < period_y, and the patch bounding box strictly "
        "inside the cell). STEP cannot represent a periodic boundary condition itself; the "
        "unit-cell footprint (substrate/metal solids) stands in for it, per `cad_export.py`.",
        "",
        "## 8. Phase extraction method",
        "",
        f"Component: `{phase_source.component}` (status: `{phase_source.status}`). {phase_source.note}",
        "",
        "## 9. Phase coverage",
        "",
        f"{coverage.n_points} successful cells at {config.f_design_hz / 1e9:g} GHz: phase range "
        f"[{coverage.phase_min_deg:.2f}, {coverage.phase_max_deg:.2f}] deg (naive, non-circular "
        f"span {coverage.naive_span_deg:.2f} deg); circularly-aware coverage = "
        f"**{coverage.coverage_deg:.2f} deg of 360 deg**. See `plots/phase_coverage_map.png` and "
        "`plots/phase_vs_geometry.png`.",
        "",
        "## 10. Largest phase gap",
        "",
        f"**{coverage.largest_gap_deg:.2f} deg**, between {coverage.largest_gap_before_deg:.2f} "
        f"and {coverage.largest_gap_after_deg:.2f} deg (wrapping through +-180 deg if the gap "
        "straddles it). This is the single largest arc of phase this library cannot represent "
        f"with any geometry in the sweep; {len(collisions)} phase region(s) "
        f"(`library_schema.PHASE_COLLISION_BIN_DEG` = {PHASE_COLLISION_BIN_DEG:g} deg bins) have more than one "
        "candidate geometry (no winner chosen automatically, req. 11).",
        "",
        "## 11. Resonance analysis",
        "",
        f"{n_resonant} cell(s) classified RESONANT, {n_non_resonant} NON_RESONANT, out of "
        f"{n_resonant + n_non_resonant} solved cells (criterion: local |S11xx(f)| maximum at a "
        f"strictly interior frequency sample with prominence >= {RESONANCE_MIN_PROMINENCE:g}, "
        "`library_classifier.detect_resonance`). See `plots/resonance_map.png`.",
        "",
        "## 12. Classification rules",
        "",
        f"- Phase usability: PHASE_USABLE if |S11xx| at design frequency >= "
        f"{PHASE_USABLE_MIN_MAG_DB:g} dB (else PHASE_UNUSABLE); {n_phase_usable} usable, "
        f"{n_phase_unusable} unusable (`library_schema.PHASE_USABLE_MIN_MAG_DB`).",
        f"- Grating risk: GRATING_RISK if a non-specular Floquet order propagates at any "
        f"frequency in the sweep (`discretization.is_specular_only`, unmodified); onset "
        f"frequency for this period = {onset_hz / 1e9 if onset_hz else float('nan'):.2f} GHz. "
        f"{n_grating_risk} of {len(classifications)} cells flagged.",
        "- Every threshold above is a named constant in `library_schema.py`, documented there "
        "and exercised at its boundary in `tests/test_library_classifier.py`.",
        "",
        "## 13. Failed cells",
        "",
        (
            f"{len(failed)} cell(s): "
            + "; ".join(f"{c.cell_id} (L={c.l_mm:g}, W={c.w_mm:g})" for c in failed[:20])
            + (f" (+{len(failed) - 20} more)" if len(failed) > 20 else "")
        )
        if failed
        else "None.",
        "",
        "## 14. STEP export status",
        "",
        f"{len(geometry_exports)} of {build.n_ok_cells} successfully-solved geometries exported "
        "as STEP; every export's substrate/metal solids independently re-validated as "
        "topologically valid solids on re-import (`cad_export.solid_validity`, "
        "`BRepCheck_Analyzer`) in `tests/test_geometry_export.py`.",
        "",
        "## 15. HFSS validation subset",
        "",
        "See `reports/HFSS_VALIDATION_SUBSET.csv` and `reports/LIBRARY_VALIDATION.md` -- a "
        "deterministic recommendation only. **No HFSS simulation has been run or imported for "
        "this library; no HFSS agreement is claimed anywhere in this report** (req. 22, 31.18).",
        "",
        "## 16. Limitations",
        "",
        '- No ground plane in the underlying solver (README "Assumptions"); S21 is not '
        "negligible for these cells, unlike a typical grounded reflectarray element -- see Sec. "
        "8's phase-source discussion.",
        '- Mesh (32x32 cells, Floquet order 32) is below this project\'s own "2-4x cells" '
        "practical accuracy guideline (`MOM_FORMULATION.md`) purely for build-time cost; every "
        "cell's a priori `metal_area_error`/`area_ok` (Phase 5, unmodified) is stored in "
        "`database/unit_cell_library.json`'s source `database.json` "
        "(and reachable via `database.Database.cell(...).validation`) rather than hidden.",
        "- STEP export cannot represent a periodic lattice-pair boundary condition; only the "
        "physical solids are exported (Sec. 7).",
        "- No interpolation or extrapolation anywhere in this library or its query API (req. 26); "
        "`find_phase_candidates`/`invert_phase` only ever return geometries actually simulated.",
        "",
        "## 17. Reproducibility",
        "",
        f"`manifest.json`'s `configuration_hash` ({manifest['configuration_hash'][:16]}...) "
        "identifies every physical input to this build. `library_builder.py` reuses "
        "`runner.run_sweep`'s checkpointed, resumable, deterministic solve (Phase 5, unmodified); "
        "`tests/test_library_builder.py` verifies numerically identical output across two builds "
        "of the same small configuration. `creation_timestamp` is the one field that is not "
        "reproducible run to run, by design (same convention as every earlier phase's file "
        "headers).",
        "",
        "## 18. HFSS validation status",
        "",
        "**HFSS validation is pending.** This report and this library are MoM-only results. No "
        "claim of agreement with HFSS (or any other full-wave tool) is made anywhere in this "
        "project until actual HFSS results are imported and quantitatively compared -- see "
        "`reports/LIBRARY_VALIDATION.md` and `VALIDATION_EXPORT.md`.",
    ]
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------ orchestration


@dataclass(frozen=True)
class FullLibraryResult:
    base_dir: Path
    build: object  # library_builder.LibraryBuildResult
    manifest: dict
    manifest_path: Path
    library_json_path: Path
    library_csv_path: Path
    phase_lookup_path: Path
    step_paths: dict
    plot_paths: dict
    hfss_subset_path: Path
    library_report_path: Path
    library_validation_path: Path
    coverage: object
    consistency: object


def build_full_library(config: LibraryConfig, base_dir, *, resume=True, progress=None, deterministic_timestamp=None):
    """The complete Phase 7 pipeline (req. 18's directory layout), end to end: solve (Phase 5,
    unmodified) -> phase extraction -> classification -> STEP export -> plots -> database/
    manifest/report files. Every step reads only what an earlier step already computed; nothing
    is solved twice and nothing is invented."""
    from .geometry_export import export_library_geometries
    from .library_builder import build_library_database
    from .library_classifier import classify_cell
    from .library_plots import generate_all_plots
    from .library_validation import (
        check_json_csv_consistency,
        render_library_validation_md,
        select_hfss_validation_subset,
        write_hfss_validation_subset_csv,
    )
    from .phase_extraction import (
        extract_design_frequency,
        phase_collision_groups,
        phase_coverage,
        phase_vs_geometry_table,
        resolve_phase_source,
    )

    base_dir = Path(base_dir)
    db_dir, geom_dir, plots_dir, reports_dir = (base_dir / d for d in ("database", "geometry", "plots", "reports"))
    for d in (db_dir, geom_dir, plots_dir, reports_dir):
        d.mkdir(parents=True, exist_ok=True)

    build = build_library_database(config, out_dir=base_dir / "_solve", resume=resume, progress=progress)
    database = build.database

    rows = extract_design_frequency(database, config)
    table = phase_vs_geometry_table(rows)
    coverage = phase_coverage(table)
    collisions = phase_collision_groups(table)
    phase_source = resolve_phase_source()

    by_lw = {(r.l_mm, r.w_mm): r for r in rows}
    classifications = [
        classify_cell(
            cell,
            database,
            by_lw[(round(cell.parameters["length_m"] * 1e3, 9), round(cell.parameters["width_m"] * 1e3, 9))],
        )
        for cell in database.cells
    ]

    library_json_path = write_library_json(database, db_dir / "unit_cell_library.json")
    library_csv_path = write_library_csv(database, db_dir / "unit_cell_library.csv")
    phase_lookup_path = write_phase_lookup_csv(table, db_dir / phase_lookup_filename(config.f_design_hz))

    step_paths = export_library_geometries(database, config, geom_dir)
    plot_paths = generate_all_plots(database, config, table, coverage, classifications, plots_dir)

    subset = select_hfss_validation_subset(config, table, classifications)
    hfss_subset_path = write_hfss_validation_subset_csv(subset, reports_dir / "HFSS_VALIDATION_SUBSET.csv")

    manifest = build_manifest(config, build, phase_source, deterministic_timestamp)
    manifest_path = write_manifest(manifest, base_dir / "manifest.json")

    consistency = check_json_csv_consistency(library_json_path, library_csv_path)

    library_validation_path = reports_dir / "LIBRARY_VALIDATION.md"
    library_validation_path.write_text(render_library_validation_md(config, subset, consistency, build))

    library_report_path = reports_dir / "LIBRARY_REPORT.md"
    library_report_path.write_text(
        render_library_report_md(
            config,
            build,
            coverage,
            collisions,
            classifications,
            step_paths,
            phase_source,
            manifest,
            consistency,
        )
    )

    return FullLibraryResult(
        base_dir,
        build,
        manifest,
        manifest_path,
        library_json_path,
        library_csv_path,
        phase_lookup_path,
        step_paths,
        plot_paths,
        hfss_subset_path,
        library_report_path,
        library_validation_path,
        coverage,
        consistency,
    )
