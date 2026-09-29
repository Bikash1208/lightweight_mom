"""Sec. 28-29, 47: environment-aware aperture -> far-field comparison, reusing Phase 8E's own
(unmodified) planar-aperture radiation integral."""

import numpy as np

from src.materials import Metal, Substrate
from src.phase10 import beam_validation as bv
from src.phase10 import generalization as gen
from src.phase10 import surrogate as surr_mod
from src.phase10.isolated_reference import IsolatedReferenceCache

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9


def _dummy_surrogate():
    cache = IsolatedReferenceCache(P, SUB, MET, 8)
    return surr_mod.EnvironmentSurrogate(cache, np.zeros(surr_mod.N_FEATURES), np.zeros(surr_mod.N_FEATURES),
                                          0.0, surr_mod.DomainSpec(frequency_hz=F))


def test_build_aperture_from_cases_uses_physical_subcell_positions():
    surrogate = _dummy_surrogate()
    grid_mm, cases = gen.run_gradient_test(P, SUB, MET, F, surrogate, mesh_cells_per_subcell=8)
    x, y, values, excluded = bv.build_aperture_from_cases(cases, P, use_surrogate=False)
    assert len(x) == len(cases)
    assert excluded == 0
    span = x.max() - x.min()
    assert span > 0 and span < gen.GRID_N * P   # positions span less than the whole outer grid


def test_run_beam_validation_matches_identical_apertures_when_surrogate_equals_mom():
    """A trivial 'surrogate' that returns the exact MoM value for every case must give an
    essentially perfect far-field match."""
    surrogate = _dummy_surrogate()
    grid_mm, cases = gen.run_gradient_test(P, SUB, MET, F, surrogate, mesh_cells_per_subcell=8)
    from src.phase10.generalization import SpatialCase
    identical_cases = [SpatialCase(c.i, c.j, c.center_lw_mm, c.neighbors_lw_mm, c.px_mom, c.px_mom)
                       for c in cases]
    result = bv.run_beam_validation(identical_cases, P, F)
    assert result["comparison"].complex_error_mean < 1e-20
