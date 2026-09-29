"""Part 3: reflectarray required-phase-map, feed illumination, and library-based cell selection."""

import numpy as np
import pytest

from src.reflectarray import designer as des


def test_broadside_beam_directly_over_a_centered_feed_is_symmetric():
    """theta0=0 (broadside): the required phase should be symmetric about the array centre
    (a purely radial function of distance from the feed's foot-point), since only R_feed(x,y)
    contributes (the steering term vanishes at theta0=0)."""
    x = np.array([[-0.05, 0.0, 0.05]] * 3)
    y = np.array([[-0.05] * 3, [0.0] * 3, [0.05] * 3])
    phase = des.required_phase_deg(x, y, 10e9, feed_xyz_m=(0.0, 0.0, 0.2), theta0_deg=0.0, phi0_deg=0.0)
    # corners equidistant from the feed's foot-point (0,0) must have identical required phase
    assert phase[0, 0] == pytest.approx(phase[0, 2], abs=1e-9)
    assert phase[0, 0] == pytest.approx(phase[2, 0], abs=1e-9)
    assert phase[0, 0] == pytest.approx(phase[2, 2], abs=1e-9)


def test_feed_at_center_gives_zero_illumination_angle_at_center():
    taper = des.feed_illumination_taper(np.array([0.0]), np.array([0.0]), feed_xyz_m=(0.0, 0.0, 0.2), q=2.0)
    assert taper[0] == pytest.approx(1.0)  # directly below the feed: theta_feed = 0, cos^q = 1


def test_illumination_taper_decreases_away_from_boresight():
    center = des.feed_illumination_taper(np.array([0.0]), np.array([0.0]), (0.0, 0.0, 0.2), q=2.0)[0]
    edge = des.feed_illumination_taper(np.array([0.1]), np.array([0.0]), (0.0, 0.0, 0.2), q=2.0)[0]
    assert edge < center


def test_assign_cells_picks_nearest_phase_candidate():
    x = np.array([[0.0, 0.008]])
    y = np.array([[0.0, 0.0]])
    required = np.array([[0.0, 90.0]])
    candidates = [
        ("c1", 2.0, 2.0, "lib", complex(1.0, 0.0)),  # phase 0 deg
        ("c2", 3.0, 3.0, "lib", complex(0.0, 1.0)),  # phase 90 deg
        ("c3", 4.0, 4.0, "lib", complex(-1.0, 0.0)),  # phase 180 deg
    ]
    assignments = des.assign_cells_from_library(x, y, required, candidates)
    assert assignments[0].chosen_l_mm == 2.0
    assert assignments[1].chosen_l_mm == 3.0
    assert assignments[0].phase_error_deg == pytest.approx(0.0, abs=1e-6)


def test_max_phase_gradient_detects_a_step():
    grid = np.array([[0.0, 10.0, 170.0]])  # last step wraps: |170-10|=160, or 200 the other way -> 160
    grad = des.max_phase_gradient_deg_per_cell(grid)
    assert grad == pytest.approx(160.0, abs=1e-6)


def test_design_reflectarray_labels_steep_gradient_when_exceeded():
    x, y = np.meshgrid(np.linspace(-0.02, 0.02, 3), np.linspace(-0.02, 0.02, 3), indexing="xy")
    candidates = [
        ("c1", 2.0, 2.0, "lib", complex(1.0, 0.0)),
        ("c2", 3.0, 3.0, "lib", np.exp(1j * np.radians(90.0))),
        ("c3", 4.0, 4.0, "lib", np.exp(1j * np.radians(-90.0))),
    ]
    design = des.design_reflectarray(
        x, y, 30e9, feed_xyz_m=(0.0, 0.0, 0.02), theta0_deg=60.0, phi0_deg=0.0, library_candidates=candidates
    )
    assert "LOCAL_PERIODICITY_APPROXIMATION" in design.labels
    assert design.max_phase_gradient_deg_per_cell >= 0.0


def test_phase_coverage_ok_is_false_above_the_threshold_even_with_a_candidate():
    """PHASE13B A4: a candidate existing is not enough -- large phase error must be flagged."""
    x = np.array([[0.0]])
    y = np.array([[0.0]])
    required = np.array([[90.0]])
    # only a 0deg candidate exists: 90deg error, well above the 30deg default threshold
    candidates = [("c1", 2.0, 2.0, "lib", complex(1.0, 0.0))]
    assignments = des.assign_cells_from_library(x, y, required, candidates)
    assert assignments[0].chosen_l_mm == 2.0  # still the (only) chosen candidate
    assert assignments[0].phase_error_deg == pytest.approx(90.0, abs=1e-6)
    assert assignments[0].phase_coverage_ok is False


def test_phase_error_stats_reflects_the_threshold():
    x = np.array([[0.0, 0.008]])
    y = np.array([[0.0, 0.0]])
    required = np.array([[0.0, 90.0]])
    candidates = [("c1", 2.0, 2.0, "lib", complex(1.0, 0.0))]  # 0deg: exact for cell 0, 90deg off for cell 1
    assignments = des.assign_cells_from_library(x, y, required, candidates)
    stats = des.phase_error_stats(assignments)
    assert stats.n_cells == 2
    assert stats.n_unreachable == 1
    assert stats.max_deg == pytest.approx(90.0, abs=1e-6)
    assert stats.max_deg >= stats.mean_deg >= 0.0


def test_design_reflectarray_labels_phase_coverage_shortfall_when_library_is_too_coarse():
    """Sparse 3-candidate library (0/+90/-90deg only, missing the whole ~180deg neighbourhood) over
    a real steered design: some cells must land far enough from any candidate to exceed the
    default 30deg threshold, so the design-level shortfall label must appear."""
    x, y = np.meshgrid(np.linspace(-0.02, 0.02, 5), np.linspace(-0.02, 0.02, 5), indexing="xy")
    candidates = [
        ("c1", 2.0, 2.0, "lib", complex(1.0, 0.0)),
        ("c2", 3.0, 3.0, "lib", np.exp(1j * np.radians(90.0))),
        ("c3", 4.0, 4.0, "lib", np.exp(1j * np.radians(-90.0))),
    ]
    design = des.design_reflectarray(
        x, y, 30e9, feed_xyz_m=(0.0, 0.0, 0.02), theta0_deg=60.0, phi0_deg=0.0, library_candidates=candidates
    )
    assert design.phase_error_stats.n_unreachable > 0
    assert any(label.startswith("PHASE_COVERAGE_SHORTFALL(") for label in design.labels)


def test_assign_cells_threshold_is_overridable():
    x = np.array([[0.0]])
    y = np.array([[0.0]])
    required = np.array([[10.0]])
    candidates = [("c1", 2.0, 2.0, "lib", complex(1.0, 0.0))]  # 10deg error
    strict = des.assign_cells_from_library(x, y, required, candidates, threshold_deg=5.0)
    loose = des.assign_cells_from_library(x, y, required, candidates, threshold_deg=15.0)
    assert strict[0].phase_coverage_ok is False
    assert loose[0].phase_coverage_ok is True
