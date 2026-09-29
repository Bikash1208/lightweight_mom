"""Sec. 23, 47: the isolated-cell and additive-correction baselines."""

import numpy as np

from src.materials import Metal, Substrate
from src.phase10 import baselines
from src.phase10.dataset import ALL_POSITIONS, EnvironmentSample
from src.phase10.isolated_reference import IsolatedReferenceCache

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9


def _sample(center, neighbors, s_complex):
    return EnvironmentSample("s", "C3", center, neighbors, "test", F, s_complex, "test")


def test_baseline_1_matches_the_real_isolated_solve():
    cache = IsolatedReferenceCache(P, SUB, MET, 8)
    b1 = baselines.build_baseline_1(cache)
    center = (2.8, 2.8)
    s = _sample(center, {p: center for p in ALL_POSITIONS}, 0j)
    expected = cache.get(center, F)
    assert b1(s) == expected


def test_baseline_1_generalizes_to_an_unseen_center():
    """No crash, no lookup-table miss -- a real solve for ANY center (Sec. 21)."""
    cache = IsolatedReferenceCache(P, SUB, MET, 8)
    b1 = baselines.build_baseline_1(cache)
    s = _sample((1.8666, 2.8), {p: (1.8666, 2.8) for p in ALL_POSITIONS}, 0j)
    result = b1(s)
    assert isinstance(result, complex)


def test_baseline_2_recovers_a_purely_additive_synthetic_relation():
    """Construct training data that is EXACTLY additive-linear in the neighbor offsets around a
    known S_C, and confirm the fit recovers it (near-)exactly."""
    cache = IsolatedReferenceCache(P, SUB, MET, 8)
    center = (2.8, 2.8)
    s0 = cache.get(center, F)
    rng = np.random.default_rng(0)
    true_beta_re = rng.normal(size=16) * 1e-9
    true_beta_im = rng.normal(size=16) * 1e-9

    def make_sample(deltas):
        neighbors = {}
        feats = []
        for pos, (dl, dw) in zip(ALL_POSITIONS, deltas):
            neighbors[pos] = (center[0] + dl, center[1] + dw)
            feats.extend((dl, dw))
        feats = np.asarray(feats)
        s = complex(s0.real + feats @ true_beta_re, s0.imag + feats @ true_beta_im)
        return _sample(center, neighbors, s)

    train = []
    for _ in range(60):
        deltas = [(rng.uniform(-0.5, 0.5), rng.uniform(-0.5, 0.5)) for _ in ALL_POSITIONS]
        train.append(make_sample(deltas))

    model = baselines.fit_baseline_2(train, cache)
    test_deltas = [(rng.uniform(-0.5, 0.5), rng.uniform(-0.5, 0.5)) for _ in ALL_POSITIONS]
    test_sample = make_sample(test_deltas)
    predicted = model.predict(test_sample)
    assert abs(predicted - test_sample.s_complex) < 1e-6 * max(abs(test_sample.s_complex), 1e-30) + 1e-20
