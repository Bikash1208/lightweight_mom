"""Sec. 15: complex-response nearest-neighbour matching, on synthetic data."""

import math

from src.phase8d_rev.response_matching import complex_distance, nearest_response_match


def test_complex_distance_is_plain_euclidean_distance():
    assert math.isclose(complex_distance(0 + 0j, 3 + 4j), 5.0)


def test_nearest_response_match_orders_by_complex_error_ascending():
    required = 1.0 + 0j
    candidates = [
        ("far", 1.0, 1.0, "lib", -1.0 + 0j),
        ("near", 2.0, 2.0, "lib", 0.9 + 0.1j),
        ("exact", 3.0, 3.0, "lib", 1.0 + 0j),
    ]
    matches = nearest_response_match(required, candidates)
    assert [m.record_id for m in matches] == ["exact", "near", "far"]
    assert matches[0].complex_error == 0.0
    assert matches[0].amplitude_error == 0.0
    assert matches[0].phase_error_deg == 0.0


def test_nearest_response_match_reports_amplitude_and_phase_error_separately():
    # same magnitude, different phase -> zero amplitude error, nonzero phase error.
    required = 1.0 + 0j
    candidates = [
        ("same_mag_diff_phase", 1.0, 1.0, "lib", complex(math.cos(math.radians(90)), math.sin(math.radians(90))))
    ]
    match = nearest_response_match(required, candidates)[0]
    assert math.isclose(match.amplitude_error, 0.0, abs_tol=1e-12)
    assert math.isclose(match.phase_error_deg, 90.0, abs_tol=1e-9)
    assert match.complex_error > 0.0  # the complex distance is NOT zero even though amplitude matches


def test_nearest_response_match_returns_every_candidate_not_a_single_winner():
    required = 0.5 + 0.5j
    candidates = [(f"c{i}", float(i), float(i), "lib", complex(i, i)) for i in range(5)]
    matches = nearest_response_match(required, candidates)
    assert len(matches) == 5
