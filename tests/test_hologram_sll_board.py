"""Phase 23 C: the job-result <-> optimiser glue on a small hexagon board: the cell-wise model
reproduces the design's own pattern, optimisation improves held-out SLL within the constraints, and
the optimised result drives the exact export."""

import numpy as np
import pytest

pytest.importorskip("fastapi")

from src.app.routers.hmsa import _run_synthesize_ports_job  # noqa: E402
from src.export.exact_board import exact_board_from_ports_result  # noqa: E402
from src.hologram.sll_board import (  # noqa: E402
    HeldOutGrid,
    _FarField,
    apply_levels,
    build_problem,
    cells_changed,
    evaluate,
    optimise,
)
from src.hologram.sll_optimizer import OptimiserConfig  # noqa: E402
from tests.test_app_phase23_hex4 import F, _Job, _request  # noqa: E402

COARSE = HeldOutGrid(theta_step_deg=3.0, phi_step_deg=6.0, cut_step_u=0.01, patch_radius_deg=3.0, patch_step_deg=0.5)


@pytest.fixture(scope="module")
def board():
    res = _run_synthesize_ports_job(_Job("p23sll"), _request({"preset": "hex4", "side_mm": 50.0, "gap_mm": 0.0}))
    beams = [[(30.0, 45.0 + 90.0 * k)] for k in range(4)]
    prob = build_problem(res, F, 2.2, 0.0009, 1.575, 6.4, beams, du_far=0.05, du_near=0.03)
    return res, prob


def test_problem_reproduces_the_designs_own_aperture_and_pattern(board):
    res, prob = board
    for p_res, bp in zip(res["ports"], prob.ports, strict=True):
        assert int(bp.owned.sum()) == p_res["n_cells"]
        field_ap = np.asarray(p_res["aperture_field_real"]) + 1j * np.asarray(p_res["aperture_field_imag"])
        from src.hologram.sll_board import _aperture

        rebuilt = _aperture(prob, bp, bp.levels0[bp.owned])
        assert np.max(np.abs(rebuilt - field_ap)) < 1e-12
    # the optimiser's cell-wise field equals the direct far field at the optimisation directions
    d = prob.directions
    ff = _FarField(prob, {p.port_id: p.state for p in prob.opt_ports}, {p.port_id: 0.0 for p in prob.opt_ports}, F, {0})
    et, ep = ff(d.theta_deg[:300], d.phi_deg[:300])
    from src.pattern.vector_aperture_pattern import copol_component

    direct = copol_component(et, ep, "rhcp")
    op = prob.opt_ports[0]
    s = op.s_states[0][np.arange(op.state.size), op.state]
    cell = (s @ op.c_co[0])[:300]
    assert np.max(np.abs(cell - direct)) <= 1e-5 * np.max(np.abs(direct))  # complex64 storage


def test_optimiser_improves_held_out_sll_within_constraints(board):
    res, prob = board
    before = evaluate(prob, grid=COARSE, scenarios=["port1"])
    cfg = OptimiserConfig(max_flips=60, q_schedule=(6.0, 12.0, 24.0), seed=0)
    levels, alpha, results = optimise(prob, "single", cfg)
    after = evaluate(prob, levels, alpha, grid=COARSE, scenarios=["port1"])
    b0, b1 = before["port1"]["beams"][0], after["port1"]["beams"][0]
    assert cells_changed(prob, levels)[0] > 0
    assert b1["sll_global_db"] <= b0["sll_global_db"] + 1e-9
    assert b1["directivity_dbi"] >= b0["directivity_dbi"] - 0.5 - 0.1  # budget + coarse-grid slack
    assert results[0].objective_after <= results[0].objective_before
    # fixed cells (under the feed) never move
    op = prob.opt_ports[0]
    assert np.all(levels[0][~op.free] == op.state[~op.free])
    # reproducible
    levels2, _, _ = optimise(prob, "single", cfg)
    assert all(np.array_equal(levels[k], levels2[k]) for k in levels)


def test_all_ports_mode_chooses_port_phases_and_result_exports(board):
    res, prob = board
    cfg = OptimiserConfig(max_flips=20, q_schedule=(6.0, 12.0), optimise_alpha=True, alpha_step_deg=30.0)
    levels, alpha, results = optimise(prob, "all", cfg)
    assert alpha[0] == 0.0 and set(alpha) == {0, 1, 2, 3}
    assert results[0].after["all"].beam_peak_db and all(
        v >= -cfg.max_gain_loss_db - 1e-6 for v in results[0].after["all"].beam_peak_db.values()
    )
    new = apply_levels(res, prob, levels, alpha, evaluation=None)
    board_exact = exact_board_from_ports_result(new, 6.4, 1.575)
    n = sum(p["patch_count"] for p in new["ports"])
    assert len(board_exact.patches) == n
    assert new["sidelobe_optimisation"]["alpha_deg"][0] == 0.0
    # patch count = cells in state 1 of the optimised levels
    assert n == sum(int(np.sum(v == 1)) for v in levels.values())


def test_per_frequency_states_robust_scenarios_and_phase_error_draws(board, monkeypatch):
    """Cells re-solved per frequency (mocked: no MoM in the fast suite), scenario instances, seeded draws."""
    import types

    import src.hologram.auto_states as auto_mod
    from src.hologram.sll_board import evaluate_phase_error_draws

    res, _prob = board
    calls = []

    def fake_auto_states(n_bits, substrate, f, period_m, fixed_L_mm=None, cache_dir=None, **_):
        calls.append((f, fixed_L_mm))
        shift = (f - F) / 1e9 * 60.0  # degrees per GHz of the state-1 phase
        g0 = np.exp(1j * np.radians(138.0))
        g1 = np.exp(1j * np.radians(-41.0 - shift))
        return types.SimpleNamespace(states=[types.SimpleNamespace(gamma=g0), types.SimpleNamespace(gamma=g1)])

    monkeypatch.setattr(auto_mod, "auto_states", fake_auto_states)
    beams = [[(30.0, 45.0 + 90.0 * k)] for k in range(4)]
    f_lo, f_hi = 10.33e9, 10.75e9
    prob = build_problem(
        res,
        F,
        2.2,
        0.0009,
        1.575,
        6.4,
        beams,
        du_far=0.05,
        du_near=0.03,
        extra_frequencies_hz=(f_lo, f_hi),
        n_phase_draws=5,
        phase_error_deg=10.0,
        draw_seed=3,
    )
    assert [c[0] for c in calls] == [f_lo, f_hi] and all(c[1] == 5.85 for c in calls)
    assert prob.frequencies_hz == (F, f_lo, f_hi) and set(prob.freq_psi) == {F, f_lo, f_hi}
    assert all(len(p.c_co) == 3 and len(p.s_states) == 3 for p in prob.opt_ports)
    names = [s.name for s in prob.scenarios]
    assert "port1" in names and "port1@10.33GHz" in names and "all@10.75GHz" in names
    assert sum(n.startswith("all@draw") for n in names) == 5
    # seeded, bounded, port 1 is the phase reference
    prob2 = build_problem(
        res,
        F,
        2.2,
        0.0009,
        1.575,
        6.4,
        beams,
        du_far=0.05,
        du_near=0.03,
        extra_frequencies_hz=(f_lo, f_hi),
        n_phase_draws=5,
        phase_error_deg=10.0,
        draw_seed=3,
    )
    assert prob.draws == prob2.draws
    assert all(d[0] == 0.0 and all(abs(v) <= 10.0 for v in d.values()) for d in prob.draws)
    # an evaluation at f_hi uses the re-solved cell phases, not the design-frequency ones
    held = evaluate(prob, frequency_hz=f_hi, grid=COARSE, scenarios=["port1"])["port1"]["beams"][0]
    prob_held = build_problem(res, F, 2.2, 0.0009, 1.575, 6.4, beams, du_far=0.05, du_near=0.03)
    held_phases = evaluate(prob_held, frequency_hz=f_hi, grid=COARSE, scenarios=["port1"])["port1"]["beams"][0]
    assert held["directivity_dbi"] != pytest.approx(held_phases["directivity_dbi"], abs=1e-6)
    # zero-error draw == the nominal all-ports evaluation
    zero = evaluate_phase_error_draws(prob, None, None, grid=COARSE, errors=[{0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0}])[0]
    nominal = evaluate(prob, grid=COARSE, scenarios=["all"])["all"]
    for a, b in zip(zero["beams"], nominal["beams"], strict=True):
        assert a["directivity_dbi"] == pytest.approx(b["directivity_dbi"], abs=1e-9)
        assert a["sll_global_db"] == pytest.approx(b["sll_global_db"], abs=1e-9)
    # robust optimisation: instances = 3 frequencies (+ draws in all-ports mode), nominal ones constrained
    cfg = OptimiserConfig(max_flips=10, q_schedule=(6.0, 12.0), max_pointing_error_deg=0.3)
    levels, alpha, runs = optimise(prob, "single", cfg, robust=True)
    assert set(runs[0].after) == {"port1", "port1@10.33GHz", "port1@10.75GHz"} or len(runs) == 4
    levels_all, _a, runs_all = optimise(prob, "all", cfg, robust=True)
    assert len(runs_all[0].after) == 3 + 5
