"""Part 0 task 2: re-stating the Phase 12B/12C/12D resonance-region mesh error at each phase's own
mesh, by scaling the (separately validated, tests/test_phase13_convergence.py) G4 power law."""

import pytest

from src.phase13 import convergence as conv
from src.phase13 import resonance_restatement as rr


def _fabricated_fit():
    """A `FitResult` matching PHASE12D_ADDENDUM.md's own measured G4 numbers (v_inf=-174.65,
    order=0.96, and `a` reconstructed so that evaluating at dx=0.2mm reproduces the measured
    5.9 deg correction) -- avoids a real 7-point Richardson solve just to test the scaling logic."""
    dx_0p2 = 0.2e-3
    order = 0.96
    correction_at_0p2 = 5.91  # -174.65 - (-168.74)
    a = correction_at_0p2 / dx_0p2**order
    return conv.FitResult(conv.STATUS_CONVERGED, -174.65, 0.5, order, -168.74, a)


def test_restate_uses_only_convergence_status_and_coefficient():
    fit = _fabricated_fit()
    results = rr.restate_resonance_region_error(reference_fit=fit)
    by_label = {r.label: r for r in results}
    assert set(by_label) == set(rr.MESHES_TO_RESTATE)


def test_restated_error_at_the_0p2mm_config_mesh_reproduces_the_calibration_point():
    """The design-point config mesh IS 40 cells / 0.2mm -- scaling the law back to its own
    calibration dx must reproduce the input correction almost exactly."""
    fit = _fabricated_fit()
    results = {r.label: r for r in rr.restate_resonance_region_error(reference_fit=fit)}
    config_result = results["phase12c_design_point_config (single-cell library mesh)"]
    assert config_result.cells == 40
    assert config_result.dx_mm == pytest.approx(0.2, abs=1e-9)
    assert config_result.scaled_correction_deg == pytest.approx(5.91, abs=0.01)


def test_coarser_meshes_have_larger_scaled_correction():
    """Phase 12C's own locality-study mesh (20 cells, 0.4mm) is coarser than the config's (40
    cells, 0.2mm); a first-order-convergent error must be LARGER there, not smaller."""
    fit = _fabricated_fit()
    results = {r.label: r for r in rr.restate_resonance_region_error(reference_fit=fit)}
    config_c = results["phase12c_design_point_config (single-cell library mesh)"].scaled_correction_deg
    locality_c = results["phase12c.locality_study (3x3/5x5 LPA suite mesh)"].scaled_correction_deg
    step_c = results["phase12d.step_curve (Task A step-size study mesh)"].scaled_correction_deg
    assert locality_c > config_c  # 0.4mm coarser than 0.2mm
    assert config_c < step_c < locality_c  # 0.2mm < 0.333mm < 0.4mm dx, monotonic in dx


def test_not_converged_reference_fit_raises_rather_than_guessing():
    bad_fit = conv.FitResult(conv.STATUS_NOT_MONOTONIC_NO_FIT, None, None, None, -168.0, None)
    with pytest.raises(ValueError):
        rr.restate_resonance_region_error(reference_fit=bad_fit)


@pytest.mark.slow
def test_fresh_richardson_call_end_to_end():
    results = rr.restate_resonance_region_error()
    assert len(results) == 3
    for r in results:
        assert r.scaled_correction_deg > 0.0
