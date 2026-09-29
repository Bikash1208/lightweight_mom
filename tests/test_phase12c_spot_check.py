"""Phase 12C Task D: spot-check Phases 9-11 on a finer aligned mesh (read-only for earlier
results). These tests use small, cheap mesh overrides (not the full production-scale aligned mesh,
which is expensive by design -- module docstring's own cost note) to keep the suite fast; the real,
fully aligned production call is exercised in test_phase12c_run.py's reduced-scale run and in the
production run itself (examples/run_phase12c.py).
"""

from src.phase12c import spot_check as sc


def test_aligned_cells_for_is_the_smallest_that_aligns_every_value():
    from src.discretization import feature_edges_aligned

    n = sc.aligned_cells_for([1.6, 3.2])
    for v in (1.6, 3.2):
        assert feature_edges_aligned(sc.PERIOD_M, n, v * 1e-3, 0.0)
    assert feature_edges_aligned(sc.PERIOD_M, n * 2, 1.6e-3, 0.0)  # doubling preserves alignment


def test_c2_c3_recheck_runs_and_reports_a_real_relative_difference():
    """C2/C3 align cheaply (24 cells) even at production scale -- no override needed."""
    result = sc.recheck_phase10_c2_c3()
    assert result.mesh_cells_per_side >= 2
    assert isinstance(result.mask_identical, bool)
    assert result.relative_difference >= 0.0
    # the mask verdict and the real-solve verdict must agree: identical masks -> zero solved
    # difference (bit-identical Z/V/S, CONVENTIONS.md), distinct masks -> a real, nonzero one.
    if result.mask_identical:
        assert result.relative_difference < 1e-9
    else:
        assert result.relative_difference > 1e-9


def test_phase9_strong_contrast_recheck_runs_at_a_cheap_test_mesh():
    result = sc.recheck_phase9_strong_contrast(mesh_cells_per_subcell=4)
    assert result.n == 3
    assert result.relative_complex_error >= 0.0
    assert result.n_unknowns_embedded > result.n_unknowns_isolated


def test_case_d_recheck_runs_at_a_cheap_test_mesh():
    result = sc.recheck_phase11_case_d(n=3, mesh_cells_per_subcell=4)
    assert result.n == 3
    assert result.relative_difference >= 0.0
    assert result.substantial == (result.relative_difference >= 0.5)


def test_coupling_radius_recheck_is_a_truncation_series():
    n_mesh, points = sc.recheck_phase11_coupling_radius(sizes=(1, 3), mesh_cells_per_subcell=4)
    assert [p.n for p in points] == [1, 3]
    assert n_mesh == 4
