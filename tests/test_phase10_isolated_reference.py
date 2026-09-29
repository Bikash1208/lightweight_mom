"""Sec. 47: the shared isolated-reference cache -- memoization and generalization to a center
geometry never seen anywhere else (Sec. 21's whole point for Baseline 2/the surrogate)."""

from src.materials import Metal, Substrate
from src.phase10.isolated_reference import IsolatedReferenceCache

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9


def test_get_is_memoized():
    cache = IsolatedReferenceCache(P, SUB, MET, 8)
    v1 = cache.get((2.8, 2.8), F)
    assert len(cache._cache) == 1
    v2 = cache.get((2.8, 2.8), F)
    assert v1 == v2
    assert len(cache._cache) == 1


def test_get_handles_a_non_grid_continuous_geometry():
    """The gradient test uses non-0.2mm-grid-aligned L values; the isolated solve must handle
    them without any special-casing (it is a general MoM solve, not a lookup table)."""
    cache = IsolatedReferenceCache(P, SUB, MET, 8)
    v = cache.get((1.8666, 2.8), F)
    assert isinstance(v, complex)
    assert abs(v) > 0
