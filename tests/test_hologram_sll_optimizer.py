"""Phase 23 C: the sidelobe optimiser on a toy aperture with a known optimum, constraint handling,
monotonicity and reproducibility."""

import itertools

import numpy as np
import pytest

from src.hologram.sll_optimizer import OptBeam, OptimiserConfig, OptPort, Scenario, optimise_sidelobes
from src.pattern.cellwise_pattern import DirectionSet


def _linear_array(n=8, spacing_wl=0.5, d_u=0.01, patch_u=0.06, patch_step=0.005):
    """A uniform 1-D array steered to broadside: `c[n, u] = exp(j 2 pi spacing n u)`, binary states +1/-1."""
    u = np.unique(np.concatenate([np.arange(-1.0, 1.0 + 1e-9, d_u), np.arange(-patch_u, patch_u + 1e-9, patch_step)]))
    u = np.round(u, 9)
    idx = np.arange(n)
    c = np.exp(1j * 2.0 * np.pi * spacing_wl * np.outer(idx - (n - 1) / 2.0, u))
    theta = np.degrees(np.arcsin(np.abs(u)))
    phi = np.where(u >= 0, 0.0, 180.0)
    tag = np.where(np.abs(u) <= patch_u + 1e-9, 0, -1)
    d = DirectionSet(theta, phi, u, np.zeros_like(u), tag)
    cone = float(np.degrees(np.arcsin(2.0 / (n * 2.0 * spacing_wl))))  # first null of the uniform array: u = 1/(n*d)
    patch = np.nonzero(tag >= 0)[0]
    port = OptPort(
        port_id=0,
        c_co=c.astype(np.complex128),
        c_cross=np.zeros((n, patch.size), dtype=np.complex128),
        s_states=np.stack([np.ones(n, dtype=complex), -np.ones(n, dtype=complex)], axis=1),
        state=np.zeros(n, dtype=int),
        free=np.ones(n, dtype=bool),
    )
    beam = OptBeam(0, 0, 0.0, 0.0, cone)
    return port, beam, d, u


def _psl_db(states, c, u, cone_u):
    s = np.where(np.asarray(states) == 0, 1.0, -1.0)
    field = s @ c
    p = np.abs(field) ** 2
    main = np.abs(u) <= cone_u
    return 10 * np.log10(p[~main].max() / p[np.abs(u) <= 0.06].max()), 10 * np.log10(p[np.abs(u) <= 0.06].max())


def test_recovers_the_known_optimum_of_a_small_array():
    n = 8
    port, beam, d, u = _linear_array(n)
    port.state[:] = (0, 0, 1, 0, 0, 1, 0, 0)  # a poor starting sequence
    cone_u = np.sin(np.radians(beam.cone_deg))
    start_psl, start_peak = _psl_db(port.state, port.c_co, u, cone_u)
    # brute force over all 2^8 sequences with the same gain budget (peak within 3 dB of the start's)
    best = (np.inf, None)
    for bits in itertools.product((0, 1), repeat=n):
        psl, peak = _psl_db(bits, port.c_co, u, cone_u)
        if peak >= start_peak - 3.0 and psl < best[0] - 1e-12:
            best = (psl, bits)
    assert best[0] < start_psl - 2.0  # there is real room to improve
    cfg = OptimiserConfig(
        max_gain_loss_db=3.0,
        max_pointing_shift_deg=1.0,
        max_level_spread_db=3.0,
        anneal_steps=2000,
        accepts_per_round=1,
        seed=1,
    )
    res = optimise_sidelobes([port], [beam], [Scenario("only", (0,), (0,))], d, cfg)
    got, peak = _psl_db(res.state[0], port.c_co, u, cone_u)
    assert peak >= start_peak - 3.0 - 1e-6
    assert got == pytest.approx(best[0], abs=0.05)  # the optimum is found (the search space is exhaustible)


def test_objective_never_increases_and_constraints_hold():
    port, beam, d, u = _linear_array(12)
    cfg = OptimiserConfig(max_gain_loss_db=0.5, max_pointing_shift_deg=0.3, accepts_per_round=1)
    res = optimise_sidelobes([port], [beam], [Scenario("only", (0,), (0,))], d, cfg)
    assert res.objective_after <= res.objective_before + 1e-12
    objs = [h[2] for h in res.history]
    # within one q the q-norm only falls; across q it is a different norm, so compare per q
    by_q = {}
    for _flips, q, obj in res.history:
        by_q.setdefault(q, []).append(obj)
    for vals in by_q.values():
        assert all(b <= a * (1 + 1e-9) for a, b in zip(vals, vals[1:], strict=False))
    assert objs
    rep = res.after["only"]
    assert rep.beam_peak_db[0] >= -cfg.max_gain_loss_db - 1e-6
    assert rep.beam_shift_deg[0] <= cfg.max_pointing_shift_deg + 1e-6
    assert rep.sll_db <= res.before["only"].sll_db + 1e-9


def test_fixed_cells_are_never_flipped_and_a_zero_budget_changes_nothing():
    port, beam, d, _u = _linear_array(10)
    port.free[[0, 4, 9]] = False
    res = optimise_sidelobes([port], [beam], [Scenario("only", (0,), (0,))], d, OptimiserConfig(max_gain_loss_db=3.0))
    assert all(res.state[0][k] == port.state[k] for k in (0, 4, 9))
    port2, beam2, d2, _ = _linear_array(10)
    res0 = optimise_sidelobes([port2], [beam2], [Scenario("only", (0,), (0,))], d2, OptimiserConfig(max_flips=0))
    assert res0.flips[0] == [] and res0.objective_after == pytest.approx(res0.objective_before)


def test_same_seed_same_result():
    cfg = OptimiserConfig(max_gain_loss_db=2.0, anneal_steps=300, seed=7)
    runs = []
    for _ in range(2):
        port, beam, d, _u = _linear_array(10)
        runs.append(optimise_sidelobes([port], [beam], [Scenario("only", (0,), (0,))], d, cfg))
    assert np.array_equal(runs[0].state[0], runs[1].state[0])
    assert runs[0].objective_after == runs[1].objective_after
    # a different seed may anneal differently but never ends worse than where the greedy phase stopped
    port, beam, d, _u = _linear_array(10)
    other = optimise_sidelobes(
        [port],
        [beam],
        [Scenario("only", (0,), (0,))],
        d,
        OptimiserConfig(max_gain_loss_db=2.0, anneal_steps=300, seed=8),
    )
    assert other.objective_after <= runs[0].objective_before
