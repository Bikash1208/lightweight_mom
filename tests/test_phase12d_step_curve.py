"""Phase 12D Task A, at reduced scale (real solves, small mesh) so this stays fast while still
exercising the whole pipeline: reference curve, center classification, 1-D/2-D step cases, and the
design-relevant reflection-coefficient metric alongside Phase 9's dipole one."""

from src.phase12c.locality_study import classify_centers
from src.phase12d import step_curve as sc
from src.discretization import feature_edges_aligned
from src.phase12c import design_point as dp


def _small_family(mesh_cells_per_subcell):
    unit_mm = 2.0 * (dp.PERIOD_M * 1e3) / mesh_cells_per_subcell
    max_k = int((dp.PERIOD_M * 1e3 * 0.9) // unit_mm)   # stay safely inside the period
    return tuple(round(unit_mm * k, 10) for k in range(2, max_k + 1))


def test_l_family_is_aligned_at_its_own_mesh():
    for l in sc.L_FAMILY_MM:
        assert feature_edges_aligned(dp.PERIOD_M, sc.EMBEDDED_MESH_CELLS_PER_SUBCELL, l * 1e-3, 0.0)


def test_build_reference_curve_reduced_scale():
    family = _small_family(12)
    curve = sc.build_reference_curve(mesh_cells_per_subcell=12, l_values_mm=family)
    assert len(curve) == len(family)
    for p in curve:
        assert 0.0 < p.s11_mag <= 1.0 + 1e-9
        assert -180.0 < p.s11_phase_deg <= 180.0


def test_classify_centers_reused_from_phase12c_works_on_this_curve():
    family = _small_family(12)
    curve = sc.build_reference_curve(mesh_cells_per_subcell=12, l_values_mm=family)
    centers = classify_centers(curve)
    assert centers.off_resonance != centers.at_resonance


def test_build_step_cases_produces_both_kinds_and_real_errors():
    family = _small_family(12)
    curve = sc.build_reference_curve(mesh_cells_per_subcell=12, l_values_mm=family)
    centers = classify_centers(curve)
    cases = sc.build_step_cases(curve, [("off", centers.off_resonance)], max_k=1,
                                 mesh_cells_per_subcell=12)
    kinds = {c.kind for c in cases}
    assert kinds == {"1d_row", "2d_gradient"}
    for c in cases:
        assert c.achieved_step_deg >= 0.0
        assert c.reflection_phase_error_deg >= 0.0
        assert c.reflection_magnitude_error_db >= 0.0
        assert c.n_unknowns_embedded > c.n_unknowns_isolated


def test_zero_step_gives_zero_error():
    """A degenerate case (neighbour == centre) must show exactly zero deviation -- a sanity
    floor for the whole metric chain."""
    family = _small_family(12)
    curve = sc.build_reference_curve(mesh_cells_per_subcell=12, l_values_mm=family)
    center = curve[len(curve) // 2]
    grid = sc.build_1d_row_grid((center.l_mm, center.l_mm), (center.l_mm, center.l_mm), (center.l_mm, center.l_mm))
    case = sc._solve_case("degenerate", "1d_row", center, grid, 0.0, 0, mesh_cells_per_subcell=12)
    assert case.reflection_phase_error_deg < 1e-6
    assert case.reflection_magnitude_error_db < 1e-6
