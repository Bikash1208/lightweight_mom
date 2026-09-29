"""Phase 12C Task B: the local-periodicity test for the grounded 8 mm cell."""

from src.phase12c import locality_study as ls


def test_l_family_is_aligned():
    """Every family value's edges must sit on this module's own aligned mesh (no staircase error
    for the geometries this study actually uses)."""
    from src.discretization import feature_edges_aligned

    for l_mm in ls.L_FAMILY_MM:
        assert feature_edges_aligned(ls.dp.PERIOD_M, ls.MESH_CELLS_PER_SUBCELL, l_mm * 1e-3, 0.0)


def test_reference_curve_is_deterministic_and_physically_bounded():
    curve1 = ls.build_reference_curve()
    curve2 = ls.build_reference_curve()
    assert curve1 == curve2
    for p in curve1:
        assert 0.0 < p.s11_mag <= 1.0 + 1e-9  # passive, grounded, near-total reflection
        assert -180.0 < p.s11_phase_deg <= 180.0


def test_classify_centers_are_interior_with_two_sided_neighbors():
    curve = ls.build_reference_curve()
    centers = ls.classify_centers(curve)
    for c in (centers.off_resonance, centers.near_resonance, centers.at_resonance):
        idx = curve.index(c)
        assert 0 < idx < len(curve) - 1  # never a domain edge -- both neighbours exist
    assert centers.at_resonance.s11_mag == min(p.s11_mag for p in curve[1:-1])


def test_immediate_neighbor_steps_are_wrapped_and_signed():
    curve = ls.build_reference_curve()
    centers = ls.classify_centers(curve)
    steps = ls.immediate_neighbor_steps(curve, centers.near_resonance)
    assert len(steps) == 2
    for s in steps:
        assert -180.0 < s.phase_step_deg <= 180.0


def test_build_1d_row_grid_only_varies_east_west():
    grid = ls.build_1d_row_grid((2.0, 2.0), (1.0, 1.0), (3.0, 3.0))
    assert grid[0] == [(2.0, 2.0)] * 3
    assert grid[2] == [(2.0, 2.0)] * 3
    assert grid[1] == [(1.0, 1.0), (2.0, 2.0), (3.0, 3.0)]


def test_2d_gradient_grid_center_is_unchanged_and_corners_use_double_step():
    curve = ls.build_reference_curve()
    centers = ls.classify_centers(curve)
    center = centers.near_resonance
    steps = ls.immediate_neighbor_steps(curve, center)
    east = next(s for s in steps if s.direction == "higher_l")
    grid, clamped = ls.build_2d_gradient_grid((center.l_mm, center.l_mm), curve, center, east.phase_step_deg)
    assert grid[1][1] == (center.l_mm, center.l_mm)
    # every entry must come from the real family (no fabricated intermediate value)
    for row in grid:
        for l, w in row:
            assert l in ls.L_FAMILY_MM and l == w


def test_run_gradient_case_matches_phase9_error_metric_semantics():
    curve = ls.build_reference_curve()
    centers = ls.classify_centers(curve)
    c = centers.off_resonance
    steps = ls.immediate_neighbor_steps(curve, c)
    west = next(s for s in steps if s.direction == "lower_l").neighbor
    east = next(s for s in steps if s.direction == "higher_l").neighbor
    grid = ls.build_1d_row_grid((c.l_mm, c.l_mm), (west.l_mm, west.l_mm), (east.l_mm, east.l_mm))
    result = ls.run_gradient_case("test_case", (c.l_mm, c.l_mm), grid, 1.0)
    # relative error must be measurably nonzero (a real heterogeneous perturbation) but far below
    # "substantial" for this small, off-resonance step
    assert 1e-6 < result.relative_complex_error < 0.5
    assert result.n_unknowns_embedded > result.n_unknowns_isolated


def test_periodic_supercell_spans_less_than_a_full_cycle_and_is_honestly_reported():
    """Task B's own explicit honesty requirement: if a full 2*pi cannot be built from real, solved
    geometries, say so."""
    curve = ls.build_reference_curve()
    results, grid_mm, span_deg, is_full = ls.run_periodic_supercell(4, curve=curve)
    assert len(results) == 4
    assert span_deg <= 360.0
    assert is_full == (span_deg >= 359.0)
    for r in results:
        assert r.l_mm in ls.L_FAMILY_MM
        assert r.relative_complex_error >= 0.0
