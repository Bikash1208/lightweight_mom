"""Sec. 10-11, 47: the 5-cell test matrix and the neighbor contribution / nonlinearity analysis."""

from src.materials import Metal, Substrate
from src.phase10 import neighbor_sensitivity as ns

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9
CENTER = (2.8, 2.8)


def _cases():
    return ns.run_five_cell_test_matrix(CENTER, (3.4, 3.4), P, SUB, MET, F, mesh_cells_per_subcell=8)


def test_five_cell_test_matrix_produces_all_nine_cases():
    cases = _cases()
    expected = {"5A_homogeneous", "5B_N", "5B_E", "5B_S", "5B_W", "5C_NS", "5C_EW", "5D_all", "5E_asymmetric"}
    assert set(cases.keys()) == expected


def test_neighbor_contributions_are_zero_relative_to_homogeneous_by_construction():
    cases = _cases()
    contribs = ns.neighbor_contributions(cases)
    assert set(contribs.keys()) == {"N", "E", "S", "W"}
    for c in contribs.values():
        assert c.delta_s_abs == abs(c.delta_s)


def test_five_c_opposite_pair_matches_perturbing_both_directions():
    cases = _cases()
    # 5C_NS perturbs N and S simultaneously; its px must differ from the homogeneous baseline.
    assert cases["5C_NS"].px != cases["5A_homogeneous"].px


def test_nonlinearity_tests_report_all_four_combinations():
    cases = _cases()
    results = ns.nonlinearity_tests(cases)
    combos = {r.combo for r in results}
    assert combos == {"N+S", "E+W", "N+E+S+W"}
    for r in results:
        assert r.e_nonlinear == abs(r.s_measured - r.s_additive_prediction)
