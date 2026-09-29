"""Phase 8A Sec. 19-26: circular phase matching and coverage-sufficiency, exercised with synthetic
required-phase arrays (brief Sec. 44 Tests 6, 9, 10, 11) since no reference-field-to-required-phase
generator exists in this project to produce real ones (`feasibility_required_phase.py`)."""

import pytest

from src.feasibility_matching import is_phase_coverage_sufficient, match_required_phase

# A synthetic library whose S11xx phase occupies a narrow, contiguous arc -- deliberately shaped
# like the real Phase 7 library's actual ~25 deg span (not a full-circle synthetic library), so
# the partial-vs-full-range distinction (Sec. 24) is tested honestly.
_NARROW_TABLE = [
    {"L_mm": 1.0, "W_mm": 1.0, "phase_deg": 60.0},
    {"L_mm": 1.2, "W_mm": 1.0, "phase_deg": 65.0},
    {"L_mm": 1.4, "W_mm": 1.0, "phase_deg": 70.0},
    {"L_mm": 1.6, "W_mm": 1.0, "phase_deg": 75.0},
    {"L_mm": 1.8, "W_mm": 1.0, "phase_deg": 80.0},
    {"L_mm": 2.0, "W_mm": 1.0, "phase_deg": 85.0},
]


def test_match_required_phase_rejects_non_positive_tolerance():
    with pytest.raises(ValueError):
        match_required_phase([70.0], _NARROW_TABLE, 0.0)


def test_match_required_phase_finds_exact_candidate():
    results = match_required_phase([70.0], _NARROW_TABLE, tolerance_deg=1.0)
    assert results[0].matched
    assert results[0].nearest_library_phase_deg == pytest.approx(70.0)
    assert results[0].phase_error_deg == pytest.approx(0.0, abs=1e-9)


# --------------------------------------------------------- Test 9: correct circular phase error


def test_circular_phase_error_near_the_wrap_boundary_is_short_not_long():
    # required=170, library phase available at -170: true difference is 20 deg the short way
    # around, not 340 deg the naive way.
    table = [{"L_mm": 1.0, "W_mm": 1.0, "phase_deg": -170.0}]
    results = match_required_phase([170.0], table, tolerance_deg=30.0)
    assert results[0].matched
    assert abs(results[0].phase_error_deg) == pytest.approx(20.0, abs=1e-6)


# --------------------------------------------------------------- Test 6: no extrapolation


def test_phase_far_outside_library_range_is_unmatched_not_extrapolated():
    # The narrow library spans 60-85 deg; 250 deg is nowhere near it.
    results = match_required_phase([250.0], _NARROW_TABLE, tolerance_deg=5.0)
    assert results[0].matched is False
    assert results[0].candidate_L_mm is not None  # nearest candidate is still *reported*...
    assert abs(results[0].phase_error_deg) > 5.0  # ...but correctly rejected as out of tolerance


def test_is_phase_coverage_sufficient_reports_unmatched_samples_explicitly():
    result = is_phase_coverage_sufficient([250.0], _NARROW_TABLE, tolerance_deg=5.0)
    assert result.sufficient is False
    assert result.n_unmatched == 1
    assert result.unmatched_samples == (250.0,)


# ------------------------------------------------------------- Test 10: partial-phase feasible


def test_required_phase_entirely_inside_narrow_library_region_is_sufficient():
    required = [61.0, 68.0, 74.0, 79.0, 84.0]  # all inside the library's 60-85 deg span
    result = is_phase_coverage_sufficient(required, _NARROW_TABLE, tolerance_deg=2.0)
    assert result.sufficient is True
    assert result.n_unmatched == 0
    assert result.max_error_deg <= 2.0


# --------------------------------------------------------------- Test 11: full-range failure


def test_required_phase_spanning_the_full_circle_is_insufficient():
    import numpy as np

    required = list(np.linspace(-180.0, 180.0, 24, endpoint=False))  # a genuine 360 deg requirement
    result = is_phase_coverage_sufficient(required, _NARROW_TABLE, tolerance_deg=5.0)
    assert result.sufficient is False
    assert result.n_unmatched > 0


def test_full_range_failure_matches_expected_phase_coverage_insufficient_scenario():
    # Mirrors the real Phase 7 library's own measured behaviour (brief Sec. 23): a library with
    # ~25 deg of coverage cannot satisfy a nominal full 0-360 deg requirement.
    import numpy as np

    required = list(np.linspace(-180.0, 180.0, 36, endpoint=False))
    result = is_phase_coverage_sufficient(required, _NARROW_TABLE, tolerance_deg=5.0)
    assert not result.sufficient
    assert result.n_unmatched >= result.n_required - 6  # only samples near 60-85 deg can match
