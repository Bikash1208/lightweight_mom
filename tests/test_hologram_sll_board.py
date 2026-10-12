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
    s = op.s_states[np.arange(op.state.size), op.state]
    cell = (s @ op.c_co)[:300]
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
