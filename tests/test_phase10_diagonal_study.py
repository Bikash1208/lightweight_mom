"""Sec. 13-15, 47: diagonal-single tests, cardinal-vs-diagonal comparison, correlated-environment
test."""

from src.materials import Metal, Substrate
from src.phase10 import diagonal_study as ds
from src.phase10 import neighbor_sensitivity as ns

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9
CENTER = (2.8, 2.8)
CONTRAST = (3.4, 3.4)


def test_diagonal_single_tests_cover_all_four_diagonals():
    cases = ds.run_diagonal_single_tests(CENTER, CONTRAST, P, SUB, MET, F, mesh_cells_per_subcell=8)
    assert set(cases.keys()) == {"D1_NE", "D1_NW", "D1_SE", "D1_SW"}


def test_diagonal_contributions_relative_to_a_shared_baseline():
    five_cell_cases = ns.run_five_cell_test_matrix(CENTER, CONTRAST, P, SUB, MET, F)
    s0 = five_cell_cases["5A_homogeneous"].px
    diag_cases = ds.run_diagonal_single_tests(CENTER, CONTRAST, P, SUB, MET, F)
    contribs = ds.diagonal_contributions(diag_cases, s0)
    assert set(contribs.keys()) == {"NE", "NW", "SE", "SW"}
    for c in contribs.values():
        assert c.delta_s_abs >= 0


def test_compare_cardinal_vs_diagonal_reports_a_finite_ratio_for_a_real_perturbation():
    five_cell_cases = ns.run_five_cell_test_matrix(CENTER, CONTRAST, P, SUB, MET, F)
    contribs = ns.neighbor_contributions(five_cell_cases)
    s0 = five_cell_cases["5A_homogeneous"].px
    diag_cases = ds.run_diagonal_single_tests(CENTER, CONTRAST, P, SUB, MET, F)
    diag_contribs = ds.diagonal_contributions(diag_cases, s0)
    cvd = ds.compare_cardinal_vs_diagonal("moderate", contribs, diag_contribs)
    assert cvd.ratio_diagonal_to_cardinal > 0
    assert cvd.mean_abs_cardinal_delta > 0 and cvd.mean_abs_diagonal_delta > 0


def test_correlated_environment_patterns_differ_for_distinct_a_b():
    result = ds.run_correlated_environment_test(CENTER, (1.0, 1.0), (4.4, 4.4), P, SUB, MET, F)
    assert result.pattern_1_px != result.pattern_2_px
    assert result.difference_abs == abs(result.difference)


def test_run_full_9cell_uniform_differs_from_homogeneous():
    five_cell_cases = ns.run_five_cell_test_matrix(CENTER, CONTRAST, P, SUB, MET, F)
    s0 = five_cell_cases["5A_homogeneous"].px
    s9, n_unknowns = ds.run_full_9cell_uniform(CENTER, CONTRAST, P, SUB, MET, F)
    assert s9 != s0
    assert n_unknowns > 0
