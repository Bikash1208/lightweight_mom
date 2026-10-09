"""Phase 22 A: automatic coded-cell states. Fast tests use an analytic stand-in for the MoM solve
(linear phase vs side, plus a mesh-size error that the Richardson step must remove); the real MoM
acceptance case (RT5880, 6.4 mm cell, 10.54 GHz) is a slow test."""

import numpy as np
import pytest

from src.hologram.auto_states import (
    LABEL_NEAR_RESONANCE,
    LABEL_UNREACHABLE_LEVEL,
    auto_states,
    effective_cells,
    empty_state_gamma,
    find_square_for_phase,
    gamma_for_square,
)
from src.materials import Substrate

F = 10.54e9
P = 6.4e-3
SUB = Substrate("RT5880", 2.2, 0.0009, 1.575e-3)


class FakeSolver:
    """phase(L) = -42 - 260 (L - 5.85) deg, limited to [-110, 150] like a real square (one 360-degree branch),
    plus a first-order mesh error +50 deg/mm * dx: the 0.2 mm mesh reads +10 deg, the 0.1 mm mesh +5 deg, so
    Richardson (v64 + (v64 - v32)) must recover the exact curve."""

    def __init__(self):
        self.calls = []

    def __call__(self, side_cells, n, substrate, frequency_hz, period_m):
        self.calls.append((side_cells, n))
        dx_mm = period_m / n * 1e3
        side_mm = side_cells * dx_mm
        phase = float(np.clip(-42.0 - 260.0 * (side_mm - 5.85), -110.0, 150.0)) + 50.0 * dx_mm
        return 0.996 * np.exp(1j * np.radians(phase))


def test_the_mesh_only_builds_even_width_squares():
    """A requested 5.85 mm square on the 0.1 mm mesh really is 5.8 mm (58 cells): why the auto search
    interpolates between mesh-aligned solves instead of solving at the snapped side."""
    assert effective_cells(5.85, 64, P) == 58
    assert effective_cells(5.8, 64, P) == 58
    assert effective_cells(5.9, 64, P) == 60
    assert effective_cells(6.0, 32, P) % 2 == 0


def test_empty_state_is_the_exact_grounded_slab():
    g = empty_state_gamma(SUB, F)
    assert abs(g) == pytest.approx(1.0, abs=2e-3)
    assert np.degrees(np.angle(g)) == pytest.approx(138.14, abs=0.05)


def test_find_square_recovers_the_curve_and_removes_the_mesh_error(tmp_path):
    solver = FakeSolver()
    side, gamma, info = find_square_for_phase(-42.0, SUB, F, P, cache_dir=tmp_path, solver=solver)
    assert side == pytest.approx(5.85, abs=1e-6)
    assert np.degrees(np.angle(gamma)) == pytest.approx(-42.0, abs=0.3)
    assert abs(gamma) == pytest.approx(0.996, abs=1e-3)
    assert info["sensitivity_deg_per_0p01mm"] == pytest.approx(2.6, abs=0.1)
    assert info["mesh_uncertainty_deg"] == pytest.approx(5.0, abs=0.2)  # the 0.2 -> 0.1 mm Richardson shift
    assert LABEL_NEAR_RESONANCE in info["labels"]
    assert all(n in (32, 64) for _, n in solver.calls)
    assert len(solver.calls) == 11 + 3  # 11-point coarse scan + 3 fine solves


def test_side_is_snapped_to_the_fabrication_grid(tmp_path):
    side, gamma, _ = find_square_for_phase(-30.0, SUB, F, P, L_grid_mm=0.05, cache_dir=tmp_path, solver=FakeSolver())
    assert side == pytest.approx(round(side / 0.05) * 0.05)
    # the reported Gamma is the curve at the SNAPPED side, not at the unsnapped optimum
    assert np.degrees(np.angle(gamma)) == pytest.approx(-42.0 - 260.0 * (side - 5.85), abs=0.3)


def test_second_run_is_served_from_the_cache(tmp_path):
    first, second = FakeSolver(), FakeSolver()
    a = auto_states(1, SUB, F, P, cache_dir=tmp_path, solver=first)
    b = auto_states(1, SUB, F, P, cache_dir=tmp_path, solver=second)
    assert a.n_solves == 14 and a.n_cached == 0
    assert b.n_solves == 0 and b.n_cached == 14 and second.calls == []
    assert b.states[1].geometry == a.states[1].geometry


def test_one_bit_pair_is_empty_plus_a_square_about_180_degrees_away(tmp_path):
    result = auto_states(1, SUB, F, P, cache_dir=tmp_path, solver=FakeSolver())
    s0, s1 = result.states
    assert s0.geometry == {"kind": "empty"} and s0.phase_deg == pytest.approx(138.14, abs=0.05)
    assert s1.geometry["kind"] == "rectangle" and s1.geometry["L_mm"] == s1.geometry["W_mm"]
    choice = result.two_state_choice()
    assert choice.achieved_separation_deg == pytest.approx(180.0, abs=5.0)
    assert choice.state_1.geometry["L_mm"] == pytest.approx(5.85, abs=0.05)


def test_n_bit_states_are_spaced_by_the_level_step(tmp_path):
    result = auto_states(2, SUB, F, P, cache_dir=tmp_path, solver=FakeSolver())
    assert [s.level for s in result.states] == [0, 1, 2, 3]
    phases = [s.phase_deg for s in result.states]
    for k in (1, 2, 3):
        got = np.degrees(np.angle(np.exp(1j * np.radians(phases[k] - phases[0] - 90.0 * k))))
        assert abs(got) < 6.0 or LABEL_UNREACHABLE_LEVEL in result.states[k].labels
    candidates = result.library_candidates()
    assert len(candidates) == 4 and candidates[0][1] == 0.0  # the empty cell has no patch
    assert np.degrees(np.angle(candidates[0][4])) == pytest.approx(0.0, abs=1e-9)


def test_fixed_side_computes_that_cells_gamma(tmp_path):
    gamma, info = gamma_for_square(5.8, SUB, F, P, cache_dir=tmp_path, solver=FakeSolver())
    assert np.degrees(np.angle(gamma)) == pytest.approx(-42.0 - 260.0 * (5.8 - 5.85), abs=0.3)
    result = auto_states(1, SUB, F, P, fixed_L_mm=5.8, cache_dir=tmp_path, solver=FakeSolver())
    assert result.states[1].geometry["L_mm"] == 5.8


def test_invalid_requests_are_refused():
    with pytest.raises(ValueError):
        auto_states(0, SUB, F, P, solver=FakeSolver())
    with pytest.raises(ValueError):
        auto_states(2, SUB, F, P, fixed_L_mm=5.8, solver=FakeSolver())


@pytest.mark.slow
def test_acceptance_reference_cell_with_the_real_mom(tmp_path):
    """RT5880 1.575 mm, 6.4 mm cell, 10.54 GHz, 1-bit: empty (+138.1 deg) + a square near 5.85 mm
    (about -42 deg +- 5), separation 180 +- 5 deg."""
    result = auto_states(1, SUB, F, P, cache_dir=tmp_path)
    s0, s1 = result.states
    assert s0.phase_deg == pytest.approx(138.14, abs=0.1)
    assert s1.geometry["L_mm"] == pytest.approx(5.85, abs=0.1)
    assert s1.phase_deg == pytest.approx(-42.0, abs=5.0)
    assert abs(s1.gamma) == pytest.approx(0.996, abs=0.01)
    assert result.two_state_choice().achieved_separation_deg == pytest.approx(180.0, abs=5.0)
    assert result.n_solves == 14
    again = auto_states(1, SUB, F, P, cache_dir=tmp_path)
    assert again.n_solves == 0 and again.solve_seconds < 5.0
