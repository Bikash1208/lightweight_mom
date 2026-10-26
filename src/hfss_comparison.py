"""Phase 6A: exact MoM reference plots + CSV for one already-exported validation case, formatted
for a direct, point-by-point overlay against exported HFSS data.

**Reads only files Phase 6 already wrote** (`validation_export.py`, unmodified by this module):
`mom_results.json` (`Database.from_json`, Phase 5, unmodified) for every S-parameter value, and
`geometry.json`/`solver_config.json` for the geometry/mesh/reference-plane facts quoted in
`HFSS_COMPARISON.md`. Nothing here re-solves `Z I = V`, re-derives an S-parameter, or invents a
number that is not already in one of those three files -- this module is a read-and-format step
on top of Phase 6's own output, not a new computation.

dB definition, exactly `20*log10(|S|)`, with one documented exception: a magnitude of exactly 0
(machine-precision cross-pol of an unrotated cell) has no finite dB value, so it is floored at
`validation_export._DB_FLOOR` dB (imported and reused directly rather than re-implemented, so
both Phase 6 dashboards and this Phase 6A export apply the same floor to the same inputs). No
normalization, no arbitrary phase shift, no unwrapping and no interpolation are applied anywhere
in this module.
"""

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import units
from .database import CONVENTIONS, Database
from .validation_export import _DB_FLOOR, _db  # same dB definition/floor as the Phase 6 dashboards

# Exact column order requested for the overlay CSV.
HFSS_CSV_COLUMNS = (
    "frequency_GHz",
    "S11xx_mag_dB",
    "S11yy_mag_dB",
    "S21xx_mag_dB",
    "S21yy_mag_dB",
    "S11xx_phase_deg",
    "S11yy_phase_deg",
    "S21xx_phase_deg",
    "S21yy_phase_deg",
    "S11xy_mag_dB",
    "S11yx_mag_dB",
    "S21xy_mag_dB",
    "S21yx_mag_dB",
)

_CO_POL = ("S11_xx", "S11_yy", "S21_xx", "S21_yy")
_CROSS_POL = ("S11_xy", "S11_yx", "S21_xy", "S21_yx")
_CELL = 0  # every validation case has exactly one cell (validation_cases.py, axes=())


@dataclass(frozen=True)
class CaseArtifacts:
    """The three already-exported files this module reads, loaded and nothing else."""

    case_dir: Path
    database: Database
    geometry: dict
    solver_config: dict


def load_case_artifacts(case_dir) -> CaseArtifacts:
    """Load `mom_results.json`, `geometry.json`, `solver_config.json` from an existing Phase 6
    case folder (`validation_export.export_case`'s output). Raises FileNotFoundError with the
    missing path if the case has not been exported yet -- this module never runs the solver."""
    case_dir = Path(case_dir)
    mom_path = case_dir / "mom_results.json"
    geometry_path = case_dir / "geometry.json"
    solver_config_path = case_dir / "solver_config.json"
    for p in (mom_path, geometry_path, solver_config_path):
        if not p.exists():
            raise FileNotFoundError(
                f"{p} not found; run validation_export.export_case (Phase 6) for this case first "
                f"-- src/hfss_comparison.py only reads already-exported files, never solves."
            )
    database = Database.from_json(mom_path)
    geometry = json.loads(geometry_path.read_text(encoding="utf-8"))
    solver_config = json.loads(solver_config_path.read_text(encoding="utf-8"))
    return CaseArtifacts(case_dir, database, geometry, solver_config)


# ------------------------------------------------------------------------------- data extraction


def _series(database, name):
    """(freq_ghz, mag, phase_deg) for one S-parameter component, straight from the stored
    Database -- `database.series`/`phase_series_deg` are Phase 5, unmodified; `unwrap=False`
    (its default) returns the principal phase exactly as stored, per item 7 of the brief."""
    freqs_hz, s = database.series(_CELL, name)
    _, phase_deg = database.phase_series_deg(_CELL, name, unwrap=False)
    return freqs_hz / 1.0e9, np.abs(s), phase_deg


def hfss_comparison_table(database: Database):
    """{column_name: array} for every `HFSS_CSV_COLUMNS` entry, plus the raw magnitude/phase used
    to build the dB/deg columns (for the plots, so the table is built once and both the CSV and
    the plots read the same numbers)."""
    table = {}
    freq_ghz = None
    mags, phases = {}, {}
    for name in _CO_POL + _CROSS_POL:
        f_ghz, mag, phase_deg = _series(database, name)
        if freq_ghz is None:
            freq_ghz = f_ghz
        elif not np.array_equal(freq_ghz, f_ghz):
            raise ValueError(f"frequency grid mismatch for {name}: every component must share one grid")
        mags[name], phases[name] = mag, phase_deg

    table["frequency_GHz"] = freq_ghz
    for name in _CO_POL:
        key = name.replace("_", "")  # "S11_xx" -> "S11xx"
        table[f"{key}_mag_dB"] = _db(mags[name])
        table[f"{key}_phase_deg"] = phases[name]
    for name in _CROSS_POL:
        key = name.replace("_", "")
        table[f"{key}_mag_dB"] = _db(mags[name])
    return table, mags, phases


def write_hfss_comparison_csv(database: Database, path):
    """Write `hfss_comparison_reference.csv`: exactly `HFSS_CSV_COLUMNS`, one row per stored
    frequency, full precision (`repr`), no rounding -- the same values the plots show."""
    table, _, _ = hfss_comparison_table(database)
    n = len(table["frequency_GHz"])
    path = Path(path)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(HFSS_CSV_COLUMNS)
        for i in range(n):
            writer.writerow([repr(float(table[col][i])) for col in HFSS_CSV_COLUMNS])
    return path


# --------------------------------------------------------------------------------------- plots


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for the Phase 6A HFSS comparison plots "
            "(src/hfss_comparison.py); install it with 'pip install -r requirements-cad.txt'."
        ) from e
    return plt


def _single_panel(plt, path, freq_ghz, curves, ylabel, title):
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for label, y, style in curves:
        ax.plot(freq_ghz, y, style, marker="o", markersize=3, label=label)
    ax.set_xlabel("frequency [GHz]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def generate_hfss_plots(database: Database, plots_dir):
    """Write the 5 requested plots into ``plots_dir`` (created if needed); return {name: path}.
    Every curve is exactly the stored MoM data at the stored frequency points: no smoothing, no
    interpolation, no normalization, no phase shift."""
    plt = _require_matplotlib()
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    _, mags, phases = hfss_comparison_table(database)
    freq_ghz, _, _ = _series(database, "S11_xx")

    paths = {
        "s11_magnitude_db": plots_dir / "hfss_s11_magnitude_db.png",
        "s21_magnitude_db": plots_dir / "hfss_s21_magnitude_db.png",
        "s11_phase_deg": plots_dir / "hfss_s11_phase_deg.png",
        "s21_phase_deg": plots_dir / "hfss_s21_phase_deg.png",
        "cross_polarization_db": plots_dir / "hfss_cross_polarization_db.png",
    }

    _single_panel(
        plt,
        paths["s11_magnitude_db"],
        freq_ghz,
        [("S11xx", _db(mags["S11_xx"]), "-"), ("S11yy", _db(mags["S11_yy"]), "--")],
        "|S11| [dB]",
        "S11 magnitude: MoM reference (exp(+jwt), 20*log10|S|)",
    )
    _single_panel(
        plt,
        paths["s21_magnitude_db"],
        freq_ghz,
        [("S21xx", _db(mags["S21_xx"]), "-"), ("S21yy", _db(mags["S21_yy"]), "--")],
        "|S21| [dB]",
        "S21 magnitude: MoM reference (exp(+jwt), 20*log10|S|)",
    )
    _single_panel(
        plt,
        paths["s11_phase_deg"],
        freq_ghz,
        [("S11xx", phases["S11_xx"], "-"), ("S11yy", phases["S11_yy"], "--")],
        "phase(S11) [deg]",
        "S11 phase: MoM reference (exp(+jwt), principal value)",
    )
    _single_panel(
        plt,
        paths["s21_phase_deg"],
        freq_ghz,
        [("S21xx", phases["S21_xx"], "-"), ("S21yy", phases["S21_yy"], "--")],
        "phase(S21) [deg]",
        "S21 phase: MoM reference (exp(+jwt), principal value)",
    )
    _single_panel(
        plt,
        paths["cross_polarization_db"],
        freq_ghz,
        [
            ("S11xy", _db(mags["S11_xy"]), "-"),
            ("S11yx", _db(mags["S11_yx"]), "--"),
            ("S21xy", _db(mags["S21_xy"]), "-."),
            ("S21yx", _db(mags["S21_yx"]), ":"),
        ],
        "magnitude [dB]",
        "Cross-polarization: MoM reference (20*log10|S|)",
    )
    return paths


# ---------------------------------------------------------------------------------- HFSS_COMPARISON.md


def render_hfss_comparison_md(artifacts: CaseArtifacts, csv_path, plot_paths):
    geo, sc = artifacts.geometry, artifacts.solver_config
    grid = sc["frequency_grid"]
    rp = sc["reference_planes"]
    case_id = geo["case_id"]

    return f"""# HFSS comparison reference — {case_id}

Authoritative MoM reference for a manual, point-by-point overlay against HFSS-exported
S-parameters of the same STEP model (`unit_cell.step`, this folder). Generated by
`src/hfss_comparison.py` (Phase 6A) from **already-exported** `mom_results.json`,
`geometry.json` and `solver_config.json` (Phase 6) only -- no S-parameter here is recomputed,
re-derived or invented; see that module's docstring.

## 1. MoM quantity <-> HFSS quantity

| This plot/column | MoM quantity | What to export from HFSS |
|---|---|---|
| `S11xx_mag_dB`, `S11xx_phase_deg` | `S11_xx` = x received at port 1 / x incident at port 1 (`SParameters.component`, `S_PARAMETERS.md` Sec. 1) | `S(Port1:x, Port1:x)` (or your tool's equivalent port/mode naming) |
| `S11yy_*` | `S11_yy` | `S(Port1:y, Port1:y)` |
| `S21xx_*` | `S21_xx` = x received at port 2 / x incident at port 1 | `S(Port2:x, Port1:x)` |
| `S21yy_*` | `S21_yy` | `S(Port2:y, Port1:y)` |
| `S11xy_mag_dB` | `S11_xy` = y received at port 1 / x incident at port 1 (cross-pol) | `S(Port1:y, Port1:x)` |
| `S11yx_mag_dB` | `S11_yx` = x received at port 1 / y incident at port 1 | `S(Port1:x, Port1:y)` |
| `S21xy_mag_dB` | `S21_xy` = y received at port 2 / x incident at port 1 | `S(Port2:y, Port1:x)` |
| `S21yx_mag_dB` | `S21_yx` = x received at port 2 / y incident at port 1 | `S(Port2:x, Port1:y)` |

Every `S_ij[received, incident]` is indexed `0 = x, 1 = y` at both ports, same global axes at both
(no axis flip on the far side) -- `S_PARAMETERS.md` Sec. 1.

## 2. HFSS port/mode mapping for x/y polarization

Assign **Floquet ports** on the top face (`z = z_top`) and bottom face (`z = z_bottom`), each with
**two modes** corresponding to the `(0,0)` Floquet order's two linear polarizations (HFSS calls
these TE00/TM00 for a Floquet port with normal incidence): mode "x" = the mode whose E-field is
along the cell's x axis, mode "y" = along y. This gives the same 4x4 layout as `mom_results.json`
over `(port1 x, port1 y, port2 x, port2 y)`. Periodic (lattice-pair / master-slave) boundaries on
the four side faces, zero phase delay (normal incidence) -- see this case's own `README.md`,
section "Reproducing this model in HFSS", for the full step-by-step recipe (unit, material,
boundary and sweep settings for this specific case).

## 3. Curve <-> S-parameter (this export's 5 plots)

1. `hfss_s11_magnitude_db.png`: `S11xx`, `S11yy` magnitude, dB.
2. `hfss_s21_magnitude_db.png`: `S21xx`, `S21yy` magnitude, dB.
3. `hfss_s11_phase_deg.png`: `S11xx`, `S11yy` phase, degrees.
4. `hfss_s21_phase_deg.png`: `S21xx`, `S21yy` phase, degrees.
5. `hfss_cross_polarization_db.png`: `S11xy`, `S11yx`, `S21xy`, `S21yx` magnitude, dB (demonstrates
   negligible cross-polarization for this unrotated rectangular patch: `{case_id}` has
   `rotation_deg = {geo["geometry"]["rotation_deg"]:g}`, so cross-pol is expected at or near the
   solver's own numerical noise floor, not a modelled physical effect).

## 4. Time convention: `exp(+j*omega*t)`

This project uses **`exp(+j*omega*t)`** throughout (`{CONVENTIONS["time_dependence"]}`,
`database.CONVENTIONS`): a phasor `F` represents `Re{{F * exp(+j*omega*t)}}`; `exp(-j*k*z)`
travels toward `+z`; a propagation delay is a **negative** phase; a lossy medium has
`Im(eps) <= 0`. Full derivation: `CONVENTIONS.md`.

## 5. If HFSS uses a different time/sign convention

Most HFSS setups default to `exp(+j*omega*t)` as well (check Tools > Options or your project's
solution setup), but if your exported data instead uses `exp(-j*omega*t)`:
**conjugate every HFSS S-parameter before overlaying** (magnitude unchanged, phase sign
reversed). Do not conjugate or otherwise alter the MoM data in this export -- it is already in
this project's stated convention (Sec. 4).

## 6. Reference-plane locations

```
Port 1 (top)    z_top    = +{units.to_mm(rp["z_top_m"]):.6g} mm = {rp["z_top_m"]!r} m
Port 2 (bottom) z_bottom = {units.to_mm(rp["z_bottom_m"]):.6g} mm = {rp["z_bottom_m"]!r} m
```

(`z_bottom = -(substrate.thickness + bottom_air_thickness)` = `-({units.to_mm(geo["substrate"]["thickness_m"]):.6g} + {units.to_mm(rp["bottom_air_thickness_m"]):.6g})` mm, `S_PARAMETERS.md` Sec. 1.) Both planes
are in air; only the specular `(0,0)` Floquet order is a port mode. Set your HFSS Floquet ports at
these exact z-coordinates (or re-reference your HFSS result to them, e.g. via port de-embedding)
before comparing phase -- a different reference plane shifts phase by `k0 * delta_z` and this
export applies no such correction for you (Sec. 8).

## 7. No phase unwrapping

Every phase column/curve is the **principal value** on `(-180, 180]`, exactly as
`Database.phase_series_deg(cell, name, unwrap=False)` (its default) returns it -- no unwrapping is
applied in this export. If you unwrap your HFSS phase for readability, unwrap this data the same
way before comparing slopes (`src/phase.py`); do not compare an unwrapped HFSS curve against this
export's wrapped one without accounting for that difference.

## 8. No empirical correction or calibration

This export applies **no normalization, no arbitrary phase shift, no smoothing, no interpolation
and no empirical correction or calibration** of any kind to the MoM data. Every value in
`{csv_path.name}` and every point on the 5 plots is exactly what `mom_results.json` (Phase 6,
itself an unmodified `Database.to_json` write, DATABASE.md) already stored for `{grid["count"]}`
frequency points from {units.to_ghz(grid["start_hz"]):.4g} to {units.to_ghz(grid["stop_hz"]):.4g} GHz
in steps of {units.to_ghz(grid["step_hz"]):.4g} GHz:
`{", ".join(f"{units.to_ghz(f):.4g}" for f in grid["frequencies_hz"])}` GHz. Export your HFSS
sweep at these same frequency points, in the same units, using the same `20*log10(|S|)` dB
definition, and overlay directly -- any difference you see is data for you to interpret, not
something this project has tuned away (`VALIDATION_EXPORT.md` Sec. 7). One exception to the plain
`20*log10(|S|)` definition: a magnitude of exactly 0 has no finite dB value, so this export floors
it at `{_DB_FLOOR:g}` dB (`validation_export._DB_FLOOR`) rather than plotting `-inf`; this only
ever applies to a cross-pol term that is already at the solver's numerical noise floor.

## Files in this export

- `{csv_path.name}` — the 5 plots' data, one row per frequency, columns: `{", ".join(HFSS_CSV_COLUMNS)}`.
{chr(10).join(f"- `plots/{p.name}`" for p in plot_paths.values())}
"""


# --------------------------------------------------------------------------- orchestration


@dataclass(frozen=True)
class HfssComparisonResult:
    case_dir: Path
    csv_path: Path
    plot_paths: dict
    md_path: Path


def export_hfss_comparison(case_dir):
    """Full Phase 6A export for one already-exported Phase 6 case folder: reads its
    `mom_results.json`/`geometry.json`/`solver_config.json`, writes `hfss_comparison_reference.csv`
    and `HFSS_COMPARISON.md` into ``case_dir`` and the 5 plots into ``case_dir/plots``."""
    artifacts = load_case_artifacts(case_dir)
    case_dir = artifacts.case_dir

    csv_path = case_dir / "hfss_comparison_reference.csv"
    write_hfss_comparison_csv(artifacts.database, csv_path)

    plot_paths = generate_hfss_plots(artifacts.database, case_dir / "plots")

    md_path = case_dir / "HFSS_COMPARISON.md"
    md_path.write_text(render_hfss_comparison_md(artifacts, csv_path, plot_paths), encoding="utf-8", newline="\n")

    return HfssComparisonResult(case_dir, csv_path, plot_paths, md_path)
