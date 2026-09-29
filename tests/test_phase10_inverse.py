"""Sec. 38-40, 47: first-stage, MoM-verified geometry inversion."""

import numpy as np

from src.materials import Metal, Substrate
from src.phase10 import inverse as inv_mod
from src.phase10 import surrogate as surr_mod
from src.phase10.dataset import ALL_POSITIONS
from src.phase10.isolated_reference import IsolatedReferenceCache

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9


def _dummy_surrogate():
    cache = IsolatedReferenceCache(P, SUB, MET, 8)
    return surr_mod.EnvironmentSurrogate(cache, np.zeros(surr_mod.N_FEATURES), np.zeros(surr_mod.N_FEATURES),
                                          0.0, surr_mod.DomainSpec(frequency_hz=F))


def test_invert_target_reaches_a_target_that_is_exactly_achievable():
    """The target IS the homogeneous-center response at (2.8, 2.8) in a homogeneous environment,
    and (2.8, 2.8) is itself one of the candidates -- the search (over a handful of the 18-value
    sweep grid points) must find a candidate the real MoM solver confirms reproduces it almost
    exactly (not necessarily (2.8, 2.8) itself: this dummy zero-coefficient "surrogate" reduces to
    the isolated-cell response alone, and two different geometries can legitimately have similar
    isolated responses -- the pipeline's job is to find and MoM-verify a genuinely close match,
    which is what is actually checked here)."""
    surrogate = _dummy_surrogate()
    center = (2.8, 2.8)
    neighbors = {p: center for p in ALL_POSITIONS}
    cache = IsolatedReferenceCache(P, SUB, MET, 8)
    target = cache.get(center, F)   # surrogate predicts exactly S_C for a homogeneous environment
    result = inv_mod.invert_target(target, neighbors, F, surrogate, P, SUB, MET, mesh_cells_per_subcell=8,
                                    candidate_lw_values_mm=[2.8])
    assert result.status == "REACHED"
    assert result.best_candidate_lw_mm == (2.8, 2.8)
    assert result.mom_relative_error < 1e-6


def test_invert_target_reports_unreachable_for_an_absurdly_large_target():
    surrogate = _dummy_surrogate()
    center = (2.8, 2.8)
    neighbors = {p: center for p in ALL_POSITIONS}
    absurd_target = complex(1.0, 0.0)   # dipole moments in this project are ~1e-9 to 1e-7 A.m
    result = inv_mod.invert_target(absurd_target, neighbors, F, surrogate, P, SUB, MET,
                                    mesh_cells_per_subcell=8, candidate_lw_values_mm=[1.0, 2.8, 4.4])
    assert result.status == inv_mod.TARGET_UNREACHABLE
