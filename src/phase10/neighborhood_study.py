"""Sec. 5-17 orchestration: runs the complete 5-cell + diagonal + long-range test matrix for one
center cell, and repeats it for every center in ``environment_geometry.CENTER_CELLS_MM`` (Sec. 6:
"additional center cells must subsequently include near-small/mid-range/near-large/anisotropic
geometry" -- done here for all five, not just the primary center, since each full scan is cheap
(measured: a handful of seconds per center at this mesh resolution)).

Four contrast levels are tested per center: ``weak``/``moderate`` (Sec. 7-8, symmetric +/-0.2/
+/-0.6 mm perturbations, direction chosen to stay in-domain) and ``strong_low``/``strong_high``
(Sec. 9, the domain's own 1.0 mm and 4.4 mm extremes, tested SEPARATELY rather than only as a mix,
because how "strong" a given extreme actually is depends on the center: for a near-small center
like C1 (1.4 mm), the 1.0 mm extreme is a mild perturbation while 4.4 mm is a huge one, and
collapsing them into a single "strong" label would hide that asymmetry).
"""

from dataclasses import dataclass

from . import diagonal_study as ds
from . import environment_geometry as eg
from . import long_range as lr
from . import neighbor_sensitivity as ns
from . import neighborhood_decision as nd

CONTRAST_LABELS = ("weak", "moderate", "strong_low", "strong_high")


def _hits_discretization_floor(contributions):
    """True if any single-neighbor contribution came back EXACTLY zero -- verified (by direct
    finer-mesh comparison, see PHASE10_ENVIRONMENT_MODEL.md) to be a mesh-staircase quantization
    artifact for small perturbations of small-patch centers at this project's standard
    mesh_cells_per_subcell=8 (0.6 mm mesh pitch, coarser than a 0.2 mm "weak" perturbation), NOT
    evidence of zero physical coupling. Flagged rather than hidden or silently re-meshed."""
    return any(c.delta_s_abs == 0.0 for c in contributions.values())


def resolve_contrast_values(center_lw_mm):
    """{label: (L_mm, W_mm)} for this center, Sec. 7-9."""
    low, high = eg.strong_neighbor_low_high(center_lw_mm)
    return {
        "weak": eg.weak_neighbor(center_lw_mm).lw_mm,
        "moderate": eg.moderate_neighbor(center_lw_mm).lw_mm,
        "strong_low": (low, low),
        "strong_high": (high, high),
    }


@dataclass(frozen=True)
class CenterStudyResult:
    center_name: str
    center_lw_mm: tuple
    contrast_values_mm: dict
    five_cell_cases: dict            # label -> {case_name: FiveCellCase}
    contributions: dict              # label -> {dir: NeighborContribution}
    nonlinearity: dict                # label -> [NonlinearityResult]
    diagonal_cases: dict              # label -> {case_name: DiagonalCase}
    diagonal_contributions: dict      # label -> {dir: NeighborContribution}
    cardinal_vs_diagonal: dict        # label -> CardinalVsDiagonal
    nine_cell_true_px: dict           # label -> complex
    predictive_errors: dict           # label -> FiveCellPredictiveError
    correlated: object                # CorrelatedEnvironmentResult
    long_range_baseline_px: complex
    long_range_cases: dict            # direction -> LongRangeCase
    decision: object                  # NeighborhoodDecision
    discretization_floor_flags: dict  # label -> bool


def run_center_study(center_name, center_lw_mm, period_m, substrate, metal, frequency_hz,
                      mesh_cells_per_subcell=8):
    contrast_values = resolve_contrast_values(center_lw_mm)

    five_cell_cases, contributions, nonlinearity = {}, {}, {}
    diagonal_cases, diagonal_contribs, cvd = {}, {}, {}
    nine_cell_true, predictive_errors, floor_flags = {}, {}, {}

    for label, contrast_lw in contrast_values.items():
        cases = ns.run_five_cell_test_matrix(center_lw_mm, contrast_lw, period_m, substrate, metal,
                                              frequency_hz, mesh_cells_per_subcell)
        five_cell_cases[label] = cases
        contribs = ns.neighbor_contributions(cases)
        contributions[label] = contribs
        nonlinearity[label] = ns.nonlinearity_tests(cases)

        s0 = cases["5A_homogeneous"].px
        d_cases = ds.run_diagonal_single_tests(center_lw_mm, contrast_lw, period_m, substrate, metal,
                                                frequency_hz, mesh_cells_per_subcell)
        diagonal_cases[label] = d_cases
        d_contribs = ds.diagonal_contributions(d_cases, s0)
        diagonal_contribs[label] = d_contribs
        cvd[label] = ds.compare_cardinal_vs_diagonal(label, contribs, d_contribs)

        s9true, _ = ds.run_full_9cell_uniform(center_lw_mm, contrast_lw, period_m, substrate, metal,
                                               frequency_hz, mesh_cells_per_subcell)
        nine_cell_true[label] = s9true
        predictive_errors[label] = nd.five_cell_predictive_error(label, s0, cases["5D_all"].px, s9true)
        floor_flags[label] = _hits_discretization_floor(contribs) or _hits_discretization_floor(d_contribs)

    # Sec. 15: correlated environment test, using the two strong extremes as A/B (the most
    # demanding, unambiguous pair of distinct geometries available).
    low, high = eg.strong_neighbor_low_high(center_lw_mm)
    correlated = ds.run_correlated_environment_test(
        center_lw_mm, (low, low), (high, high), period_m, substrate, metal, frequency_hz,
        mesh_cells_per_subcell,
    )

    # Sec. 17: long-range (2-period) test, using the "strong_low" contrast as the representative,
    # most-demanding perturbation (a documented scoping choice, not tested at every contrast level
    # to keep the per-center cost bounded).
    lr_baseline = lr.solve_5x5_homogeneous(center_lw_mm, period_m, substrate, metal, frequency_hz,
                                            mesh_cells_per_subcell)
    lr_cases = lr.run_long_range_test(center_lw_mm, (low, low), lr_baseline, period_m, substrate,
                                       metal, frequency_hz, mesh_cells_per_subcell)

    decision = nd.decide_neighborhood_size(predictive_errors, lr_cases, lr_baseline)

    return CenterStudyResult(
        center_name, center_lw_mm, contrast_values, five_cell_cases, contributions, nonlinearity,
        diagonal_cases, diagonal_contribs, cvd, nine_cell_true, predictive_errors, correlated,
        lr_baseline, lr_cases, decision, floor_flags,
    )


def run_full_neighborhood_study(centers_mm, period_m, substrate, metal, frequency_hz,
                                 mesh_cells_per_subcell=8):
    """Sec. 6: repeats ``run_center_study`` for every center cell given (default: all five of
    ``environment_geometry.CENTER_CELLS_MM``)."""
    return {
        name: run_center_study(name, lw_mm, period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell)
        for name, lw_mm in centers_mm.items()
    }
