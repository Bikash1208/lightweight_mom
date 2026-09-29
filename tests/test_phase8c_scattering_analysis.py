"""Phase 8C Sec. 17-21, 41, 55: full-library scattering extraction, circular coverage (with the
[-179,+179] near-wrap case), and complex-vs-phase state separation."""

import math

import pytest

from src.phase8c.scattering_analysis import (
    ScatteringRow,
    candidate_coverage,
    distinct_complex_and_phase_states,
    extract_component_table,
)


def _row(l, w, real, imag):
    z = complex(real, imag)
    mag = abs(z)
    mag_db = 20.0 * math.log10(mag) if mag > 0.0 else -300.0
    return ScatteringRow(l, w, real, imag, mag, mag_db, math.degrees(math.atan2(imag, real)))


# ---------------------------------------------------------------------- Test: circular near-wrap

def test_circular_coverage_treats_minus179_plus179_as_close_not_358_apart():
    rows = [_row(1, 1, *_polar(1.0, -179.0)), _row(1, 2, *_polar(1.0, 179.0))]
    cov = candidate_coverage(rows, "S11xx")
    # These two points are 2 deg apart the short way, so the largest gap should be ~358 deg
    # (everything else empty) and coverage only ~2 deg -- NOT the naive 358 deg reading of
    # max(phase)-min(phase) being mistaken for "coverage".
    assert cov.largest_gap_deg == pytest.approx(358.0, abs=0.01)
    assert cov.coverage_deg == pytest.approx(2.0, abs=0.01)


def _polar(mag, phase_deg):
    import math

    rad = math.radians(phase_deg)
    return mag * math.cos(rad), mag * math.sin(rad)


# ------------------------------------------------------------------------- Test: duplicate states

def test_distinct_states_groups_bit_identical_complex_values_deterministically():
    rows = [_row(1, 1, 0.5, 0.1), _row(1.2, 1, 0.5, 0.1), _row(1.4, 1, 0.5, 0.2)]
    result = distinct_complex_and_phase_states(rows)
    assert result.n_geometries == 3
    assert result.n_distinct_complex_states == 2   # first two are identical complex values


# --------------------------------------------------------------- Test: phase/magnitude separation

def test_same_phase_different_magnitude_are_not_the_same_complex_state():
    rows = [_row(1, 1, *_polar(0.5, 30.0)), _row(1.2, 1, *_polar(0.9, 30.0))]
    result = distinct_complex_and_phase_states(rows)
    assert result.n_distinct_phase_states == 1     # same phase
    assert result.n_distinct_complex_states == 2    # but different magnitude -> different complex state


def test_candidate_coverage_reports_magnitude_statistics_not_only_phase():
    rows = [_row(1, 1, *_polar(0.2, 0.0)), _row(1.2, 1, *_polar(0.8, 10.0)), _row(1.4, 1, *_polar(0.5, 20.0))]
    cov = candidate_coverage(rows, "S11xx")
    assert cov.magnitude_min == pytest.approx(0.2, abs=1e-6)
    assert cov.magnitude_max == pytest.approx(0.8, abs=1e-6)
    assert cov.magnitude_median == pytest.approx(0.5, abs=1e-6)


# -------------------------------------------------------------------------------- real-data tests

def test_extract_component_table_real_s11_and_s21(phase8c_full_database):
    s11 = extract_component_table(phase8c_full_database, "S11_xx", 100)
    s21 = extract_component_table(phase8c_full_database, "S21_xx", 100)
    assert len(s11) == 324 and len(s21) == 324


def test_real_s21_has_different_coverage_than_s11_not_selected_by_size(phase8c_full_database):
    # This test only checks the numbers differ (a real, measured fact); it must NOT be read as
    # "therefore S21 is preferred" (Sec. 1's explicit prohibition) -- no selection happens here.
    s11 = extract_component_table(phase8c_full_database, "S11_xx", 100)
    s21 = extract_component_table(phase8c_full_database, "S21_xx", 100)
    cov11 = candidate_coverage(s11, "S11xx")
    cov21 = candidate_coverage(s21, "S21xx")
    assert cov11.coverage_deg != pytest.approx(cov21.coverage_deg, abs=0.01)


def test_real_distinct_states_s11_and_s21_are_not_assumed_equal(phase8c_full_database):
    s11 = extract_component_table(phase8c_full_database, "S11_xx", 100)
    s21 = extract_component_table(phase8c_full_database, "S21_xx", 100)
    r11 = distinct_complex_and_phase_states(s11)
    r21 = distinct_complex_and_phase_states(s21)
    assert r11.n_geometries == r21.n_geometries == 324
    # (Sec. 41: do not assume the 169-state S11 result applies equally to S21 -- checked, not assumed.)
