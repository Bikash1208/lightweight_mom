"""Sec. 17-25: grid builders, LPA error calculation, classification -- and the real, measured
6-scenario suite against the actual solver (fast: <10s at this module's reduced mesh resolution)."""

import pytest

from src.materials import Metal, Substrate
from src.phase9.local_periodicity import (
    CENTER_LW_MM,
    SUBSTANTIAL_DEVIATION_RELATIVE,
    build_3x3_grid,
    classify_local_periodicity,
    compute_lpa_error,
    extend_to_5x5,
    run_local_periodicity_suite,
)
from src.phase9.dipole_extraction import DipoleMoment

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3


def test_build_3x3_grid_places_center_at_1_1():
    for kind in ("weak", "moderate", "strong", "symmetric"):
        grid = build_3x3_grid(kind)
        assert grid[1][1] == pytest.approx((CENTER_LW_MM[0] * 1e-3, CENTER_LW_MM[1] * 1e-3))


def test_build_3x3_grid_rejects_unknown_kind():
    with pytest.raises(ValueError):
        build_3x3_grid("nonexistent")


def test_symmetric_grid_is_actually_symmetric_under_180_degree_rotation():
    grid = build_3x3_grid("symmetric")
    for i in range(3):
        for j in range(3):
            assert grid[i][j] == grid[2 - i][2 - j]


def test_extend_to_5x5_keeps_the_inner_3x3_unchanged():
    inner = build_3x3_grid("moderate")
    outer_ring = ((1.0, 1.0),) * 16
    grid5 = extend_to_5x5(inner, outer_ring)
    for i in range(3):
        for j in range(3):
            assert grid5[i + 1][j + 1] == inner[i][j]


def test_extend_to_5x5_fills_every_outer_position():
    inner = build_3x3_grid("weak")
    outer_ring = ((2.0, 2.0),) * 16
    grid5 = extend_to_5x5(inner, outer_ring)
    assert all(grid5[i][j] is not None for i in range(5) for j in range(5))


def test_compute_lpa_error_zero_for_identical_dipoles():
    dipole = DipoleMoment(1 + 2j, 0j)
    ce, me, pe, rel = compute_lpa_error(dipole, dipole)
    assert ce == pytest.approx(0.0, abs=1e-12)
    assert me == pytest.approx(0.0, abs=1e-12)
    assert pe == pytest.approx(0.0, abs=1e-9)
    assert rel == pytest.approx(0.0, abs=1e-12)


def test_classify_local_periodicity_validated_when_all_within_noise_floor():
    class _Fake:
        def __init__(self, rel):
            self.relative_complex_error = rel

    results = [_Fake(1e-10), _Fake(1e-9)]
    assert classify_local_periodicity(results) == "VALIDATED"


def test_classify_local_periodicity_invalidated_when_one_case_is_substantial():
    class _Fake:
        def __init__(self, rel):
            self.relative_complex_error = rel

    results = [_Fake(1e-8), _Fake(SUBSTANTIAL_DEVIATION_RELATIVE + 0.1)]
    assert classify_local_periodicity(results) == "INVALIDATED"


def test_classify_local_periodicity_unvalidated_for_the_ambiguous_middle_ground():
    class _Fake:
        def __init__(self, rel):
            self.relative_complex_error = rel

    results = [_Fake(0.05), _Fake(0.1)]
    assert classify_local_periodicity(results) == "UNVALIDATED"


def test_run_local_periodicity_suite_produces_all_six_required_scenarios():
    isolated, results = run_local_periodicity_suite(SUB, MET, P, 12.5e9, mesh_cells_per_subcell=6, floquet_factor=2)
    names = {r.name for r in results}
    assert names == {"3x3_weak", "3x3_moderate", "3x3_strong", "3x3_symmetric", "5x5_moderate", "5x5_strong"}
    assert isolated.n_unknowns > 0


def test_run_local_periodicity_suite_strong_contrast_shows_a_real_substantial_deviation():
    """This is the real, measured Phase 9 finding: strong-contrast neighbors produce a deviation
    comparable to or exceeding the signal's own magnitude."""
    isolated, results = run_local_periodicity_suite(SUB, MET, P, 12.5e9, mesh_cells_per_subcell=6, floquet_factor=2)
    strong = next(r for r in results if r.name == "3x3_strong")
    assert strong.relative_complex_error > SUBSTANTIAL_DEVIATION_RELATIVE
