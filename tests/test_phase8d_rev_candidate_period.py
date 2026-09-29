"""Sec. 11-13: HMSARequirements schema, the candidate-period sampling heuristic, and the
CANDIDATE-vs-VALIDATED distinction."""

import math

import pytest

from src.phase8c import data_source
from src.phase8d_rev import library_entry as _library_entry
from src.phase8d_rev import universal_database as _udb
from src.phase8d_rev.candidate_period import (
    HMSARequirements,
    candidate_periods_from_sampling,
    required_transverse_spatial_frequency_per_m,
    validate_period,
)


def test_hmsa_requirements_defaults_have_no_fixed_period():
    req = HMSARequirements()
    assert req.fixed_period_mm is None


def test_required_transverse_spatial_frequency_needs_at_least_two_beams():
    with pytest.raises(ValueError):
        required_transverse_spatial_frequency_per_m(((0.0, 0.0),), 12.5e9)


def test_required_transverse_spatial_frequency_zero_for_identical_beams():
    k = required_transverse_spatial_frequency_per_m(((0.0, 0.0), (0.0, 0.0)), 12.5e9)
    assert k == 0.0


def test_required_transverse_spatial_frequency_nonzero_for_different_beams():
    k = required_transverse_spatial_frequency_per_m(((0.0, 0.0), (30.0, 0.0)), 12.5e9)
    assert k > 0.0


def test_candidate_periods_from_sampling_rejects_nonpositive_k():
    with pytest.raises(ValueError):
        candidate_periods_from_sampling(0.0)


def test_candidate_periods_from_sampling_are_all_candidate_heuristic_never_validated():
    candidates = candidate_periods_from_sampling(100.0, n_candidates=3)
    assert len(candidates) == 3
    assert all(c["status"] == "CANDIDATE_HEURISTIC" for c in candidates)
    # Nyquist bound: P <= pi/k_max; the largest-fraction candidate should equal that bound.
    p_max_mm = math.pi / 100.0 * 1e3
    assert math.isclose(candidates[0]["p_mm"], p_max_mm, rel_tol=1e-9)


@pytest.mark.skipif(not data_source.DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present")
def test_validate_period_reports_validated_for_the_registered_benchmark_period(phase8c_full_database):
    entry = _library_entry.from_database(phase8c_full_database, "lib1", "rectangular_patch",
                                          _library_entry.REFERENCE_VALIDATION_LIBRARY)
    udb = _udb.UniversalDatabase()
    udb.register(entry)
    result = validate_period(udb, "rectangular_patch", entry.p_x_mm, entry.p_y_mm, frequency_hz=12.5e9)
    assert result.status == "VALIDATED"
    assert entry.library_id in result.matching_library_ids


@pytest.mark.skipif(not data_source.DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present")
def test_validate_period_reports_candidate_only_for_an_unregistered_period(phase8c_full_database):
    entry = _library_entry.from_database(phase8c_full_database, "lib1", "rectangular_patch",
                                          _library_entry.REFERENCE_VALIDATION_LIBRARY)
    udb = _udb.UniversalDatabase()
    udb.register(entry)
    result = validate_period(udb, "rectangular_patch", 3.2, 3.2, frequency_hz=12.5e9)
    assert result.status == "CANDIDATE_ONLY"
    assert result.matching_library_ids == ()


@pytest.mark.skipif(not data_source.DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present")
def test_validate_period_reports_grating_safety_independently_of_validated_status(phase8c_full_database):
    entry = _library_entry.from_database(phase8c_full_database, "lib1", "rectangular_patch",
                                          _library_entry.REFERENCE_VALIDATION_LIBRARY)
    udb = _udb.UniversalDatabase()
    udb.register(entry)
    result = validate_period(udb, "rectangular_patch", 3.2, 3.2, frequency_hz=12.5e9)
    assert result.status == "CANDIDATE_ONLY"
    assert result.grating_safe is True   # 3.2 mm at 12.5 GHz is well below grating onset either way
