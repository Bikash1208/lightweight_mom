"""Part 0 task 3 (PHASE12D_ADDENDUM.md): ingest G4 (``g4.s36p``, the zero-thickness Perfect-E
sheet run that overwrote the original 35 um run) and G6 (``g6.s36p``, 35 um patch) alongside G1,
comparing HFSS to both the fixed 0.2 mm MoM (Phase 12D's own baseline) and the Richardson-
converged MoM (``phase13.convergence``, this phase's new capability). Records the resonance-region
error budget (mesh vs patch-thickness vs unexplained residual) and produces a new, honest
single-cell HFSS gate status -- this module *supersedes* Phase 12D's HFSS status for G4/G6 without
editing any Phase 12D output in place (the project's standing rule; see PHASE12D_LOCAL_PERIODICITY.md
and PHASE12C_GROUNDED_LOCALITY.md for the same pattern one phase earlier).

Reuses, unmodified: ``phase12d.hfss_ingest.compare_single_cell`` (works for any case_id given its
own Touchstone path) and ``phase13.convergence.richardson_extrapolate_s11xx``.

The one piece of HFSS data that genuinely no longer exists on disk is the original 35 um-patch G4
run's raw Touchstone file (overwritten by the sheet run, per PHASE12D_ADDENDUM.md's own "Update"
note) -- only its already-processed comparison CSV (``G4_thick35um_comparison.csv``) survives.
That CSV is loaded as already-computed HFSS-vs-0.2mm-MoM data (not independently re-derivable),
and is used only for the "thickness" leg of the error budget, never presented as something this
phase re-measured from raw S-parameters.
"""

import csv
from dataclasses import dataclass

import numpy as np

from ..phase import wrapped_difference_deg
from ..phase12c import design_point as dp
from ..phase12d.hfss_ingest import EXPORT_ROOT, MOM_REFERENCE_CSV, compare_single_cell
from . import convergence as conv

HFSS_RESULTS_DIR = EXPORT_ROOT / "hfss_results"
G4_SHEET_TOUCHSTONE = HFSS_RESULTS_DIR / "g4.s36p"  # zero-thickness Perfect E sheet
G4_THICK_COMPARISON_CSV = HFSS_RESULTS_DIR / "G4_thick35um_comparison.csv"  # 35um run, raw file lost
G6_TOUCHSTONE = HFSS_RESULTS_DIR / "g6.s36p"  # 35um patch

# The addendum's own representative frequencies for the mesh-convergence study (cheap: 3 points,
# not the full 21/401-point grid -- a full Richardson family (7 mesh resolutions) at every
# frequency would be ~150x the cost for no added information at this stage).
RESONANCE_BUDGET_FREQUENCIES_GHZ = (11.5, 12.5, 13.5)

L_MM = {"G1": 2.0, "G2": 4.8, "G3": 5.6, "G4": 6.0, "G5": 6.8, "G6": 7.6}

RESONANCE_TOLERANCE_DEG = 5.0  # generous, resonance-region-specific tolerance (see report)


@dataclass(frozen=True)
class ConvergedPoint:
    frequency_ghz: float
    mesh_0p2mm_phase_deg: float  # the 40-cell (dx=0.2mm) point, Phase 12D/12C's production mesh
    extrapolated_phase_deg: float | None
    extrapolated_phase_error_deg: float | None
    fitted_order: float | None
    status: str


@dataclass(frozen=True)
class ErrorBudget:
    frequency_ghz: float
    mesh_deg: float  # |phi(0.2mm) - phi_inf|
    thickness_deg: float | None  # |phi_HFSS_thick - phi_HFSS_sheet|, None if no thick data
    residual_deg: float | None  # |phi_HFSS_sheet - phi_inf|, None if extrapolation unavailable
    total_hfss_vs_0p2mm_deg: float  # |phi_HFSS_thick(or sheet) - phi(0.2mm)|, the originally-observed gap


def _load_precomputed_comparison_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def converged_mom_curve(l_mm, frequencies_ghz=RESONANCE_BUDGET_FREQUENCIES_GHZ, cell_counts=conv.DEFAULT_CELL_COUNTS):
    """Richardson-extrapolated S11_xx phase at each of `frequencies_ghz`, for the isolated cell at
    `l_mm` -- reuses `convergence.richardson_extrapolate_s11xx` unmodified. Also returns the plain
    0.2 mm (40-cell) point from the same family, so both numbers come from ONE consistent set of
    solves (no separate, possibly-inconsistent lookup)."""
    unit_cell = dp.unit_cell_at(l_mm * 1e-3)
    points = []
    for f_ghz in frequencies_ghz:
        result = conv.richardson_extrapolate_s11xx(unit_cell, f_ghz * 1e9, dp.REFERENCE_PLANES, cell_counts=cell_counts)
        mesh_40 = next(p for p in result.points if p.cells == 40)
        points.append(
            ConvergedPoint(
                frequency_ghz=f_ghz,
                mesh_0p2mm_phase_deg=mesh_40.s11_phase_deg,
                extrapolated_phase_deg=result.phase.value_inf,
                extrapolated_phase_error_deg=result.phase.value_inf_error,
                fitted_order=result.phase.order,
                status=result.phase.status,
            )
        )
    return tuple(points)


def error_budget_for_g4(sheet_thick_rows=None):
    """The three-way error budget at G4 (L = 6.0 mm, resonance): mesh (0.2mm -> converged),
    patch thickness (sheet HFSS run vs the lost-raw-file 35 um run, via its surviving comparison
    CSV), and the unexplained residual (converged MoM vs HFSS sheet -- the apples-to-apples
    zero-thickness-vs-zero-thickness comparison)."""
    converged = {p.frequency_ghz: p for p in converged_mom_curve(L_MM["G4"])}
    sheet_rows = {
        round(float(r["f_GHz"]), 6): r
        for r in _load_precomputed_comparison_csv(HFSS_RESULTS_DIR / "G4_sheet_comparison.csv")
    }
    thick_rows = {round(float(r["f_GHz"]), 6): r for r in _load_precomputed_comparison_csv(G4_THICK_COMPARISON_CSV)}

    budgets = []
    for f_ghz in RESONANCE_BUDGET_FREQUENCIES_GHZ:
        cp = converged[f_ghz]
        sheet_row = sheet_rows.get(round(f_ghz, 6))
        thick_row = thick_rows.get(round(f_ghz, 6))
        if cp.status != conv.STATUS_CONVERGED or sheet_row is None:
            budgets.append(ErrorBudget(f_ghz, float("nan"), None, None, float("nan")))
            continue
        hfss_sheet_phase = cp.mesh_0p2mm_phase_deg + float(sheet_row["dphase_mode1_deg"])
        mesh_deg = abs(cp.extrapolated_phase_deg - cp.mesh_0p2mm_phase_deg)
        residual_deg = abs(
            wrapped_difference_deg(np.array([hfss_sheet_phase]), np.array([cp.extrapolated_phase_deg]))[0]
        )
        thickness_deg = None
        total_deg = abs(float(sheet_row["dphase_mode1_deg"]))
        if thick_row is not None:
            hfss_thick_phase = cp.mesh_0p2mm_phase_deg + float(thick_row["dphase_mode1_deg"])
            thickness_deg = abs(wrapped_difference_deg(np.array([hfss_thick_phase]), np.array([hfss_sheet_phase]))[0])
            total_deg = abs(float(thick_row["dphase_mode1_deg"]))
        budgets.append(ErrorBudget(f_ghz, mesh_deg, thickness_deg, residual_deg, total_deg))
    return tuple(budgets)


def ingest_g4_sheet():
    """G4 vs the zero-thickness Perfect-E sheet HFSS run (``g4.s36p``) -- the apples-to-apples
    case, reusing `compare_single_cell` unmodified with G4's own Touchstone path."""
    return compare_single_cell("G4", touchstone_path=G4_SHEET_TOUCHSTONE, mom_csv_path=MOM_REFERENCE_CSV)


def ingest_g6():
    """G6 vs its 35 um-patch HFSS run (``g6.s36p``) -- no separate sheet run exists for G6."""
    return compare_single_cell("G6", touchstone_path=G6_TOUCHSTONE, mom_csv_path=MOM_REFERENCE_CSV)
