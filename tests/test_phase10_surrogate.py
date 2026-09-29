"""Sec. 22-24, 27, 35-37, 47: the ridge surrogate, its domain check, and (de)serialization."""

import numpy as np

from src.materials import Metal, Substrate
from src.phase10 import surrogate as surr_mod
from src.phase10.dataset import ALL_POSITIONS, EnvironmentSample
from src.phase10.isolated_reference import IsolatedReferenceCache

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9


def _sample(center, neighbors, s_complex=0j, frequency_hz=F):
    return EnvironmentSample("s", "C3", center, neighbors, "test", frequency_hz, s_complex, "test")


def test_nonlinear_features_has_the_documented_length():
    s = _sample((2.8, 2.8), {p: (3.0, 3.0) for p in ALL_POSITIONS})
    feats = surr_mod.nonlinear_features(s)
    assert feats.shape == (surr_mod.N_FEATURES,)


def test_nonlinear_features_are_all_zero_for_a_homogeneous_sample():
    s = _sample((2.8, 2.8), {p: (2.8, 2.8) for p in ALL_POSITIONS})
    feats = surr_mod.nonlinear_features(s)
    assert np.allclose(feats, 0.0)


def test_domain_spec_rejects_out_of_range_geometry():
    d = surr_mod.DomainSpec(frequency_hz=F)
    ok = _sample((2.8, 2.8), {p: (2.8, 2.8) for p in ALL_POSITIONS})
    bad_center = _sample((0.5, 2.8), {p: (2.8, 2.8) for p in ALL_POSITIONS})
    bad_neighbor = _sample((2.8, 2.8), {**{p: (2.8, 2.8) for p in ALL_POSITIONS}, "N": (5.0, 5.0)})
    bad_freq = _sample((2.8, 2.8), {p: (2.8, 2.8) for p in ALL_POSITIONS}, frequency_hz=10.0e9)
    assert d.contains(ok)
    assert not d.contains(bad_center)
    assert not d.contains(bad_neighbor)
    assert not d.contains(bad_freq)


def test_predict_returns_out_of_domain_sentinel_not_a_number():
    cache = IsolatedReferenceCache(P, SUB, MET, 8)
    model = surr_mod.EnvironmentSurrogate(
        cache, np.zeros(surr_mod.N_FEATURES), np.zeros(surr_mod.N_FEATURES), 0.0, surr_mod.DomainSpec(frequency_hz=F)
    )
    bad = _sample((0.1, 2.8), {p: (2.8, 2.8) for p in ALL_POSITIONS})
    assert model.predict(bad) == surr_mod.OUT_OF_DOMAIN


def test_save_and_load_round_trips_predictions_exactly(tmp_path):
    cache = IsolatedReferenceCache(P, SUB, MET, 8)
    rng = np.random.default_rng(0)
    model = surr_mod.EnvironmentSurrogate(
        cache,
        rng.normal(size=surr_mod.N_FEATURES),
        rng.normal(size=surr_mod.N_FEATURES),
        0.01,
        surr_mod.DomainSpec(frequency_hz=F),
    )
    s = _sample((2.8, 2.8), {p: (3.0, 2.6) for p in ALL_POSITIONS})
    before = model.predict(s)

    path = tmp_path / "model.json"
    model.save(path)
    loaded = surr_mod.EnvironmentSurrogate.load(path)
    after = loaded.predict(s)
    assert before == after


def test_fit_surrogate_is_reproducible():
    cache1 = IsolatedReferenceCache(P, SUB, MET, 8)
    cache2 = IsolatedReferenceCache(P, SUB, MET, 8)
    rng = np.random.default_rng(1)
    center = (2.8, 2.8)

    def make(i):
        neighbors = {p: (center[0] + rng.uniform(-0.3, 0.3), center[1] + rng.uniform(-0.3, 0.3)) for p in ALL_POSITIONS}
        return _sample(center, neighbors, complex(rng.normal(), rng.normal()) * 1e-9)

    train = [make(i) for i in range(30)]
    val = [make(i) for i in range(10)]
    m1, lam1, err1 = surr_mod.fit_surrogate(train, val, F, cache1, ridge_candidates=(0.1,))
    m2, lam2, err2 = surr_mod.fit_surrogate(train, val, F, cache2, ridge_candidates=(0.1,))
    assert lam1 == lam2
    assert np.allclose(m1.beta_re, m2.beta_re)
    assert np.allclose(m1.beta_im, m2.beta_im)
