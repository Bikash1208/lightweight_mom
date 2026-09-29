"""Part 0 task 1: Richardson extrapolation, checked against PHASE12D_ADDENDUM.md's own hand-
measured numbers at G4 (L = 6.0 mm, resonance) and a fast, cheap synthetic check off-resonance."""

import numpy as np
import pytest

from src.phase12c import design_point as dp
from src.phase13 import convergence as conv


def test_solve_at_resolution_reproduces_addendum_g4_40cell_point():
    """PHASE12D_ADDENDUM.md: G4 at 40 cells, floquet 40, 12.5 GHz -> -168.74 deg."""
    uc = dp.unit_cell_at(6.0e-3)
    point = conv.solve_at_resolution(uc, 12.5e9, dp.REFERENCE_PLANES, cells=40)
    assert point.dx_m == pytest.approx(8.0e-3 / 40)
    assert point.s11_phase_deg == pytest.approx(-168.74, abs=0.05)


@pytest.mark.slow
def test_richardson_extrapolation_reproduces_addendum_g4_phi_inf():
    """Full 7-point family (40..96 cells): addendum reports p ~ 0.96, phi_inf ~ -174.65 deg."""
    uc = dp.unit_cell_at(6.0e-3)
    result = conv.richardson_extrapolate_s11xx(uc, 12.5e9, dp.REFERENCE_PLANES)
    assert result.phase.status == conv.STATUS_CONVERGED
    assert result.phase.order == pytest.approx(0.96, abs=0.1)
    assert result.phase.value_inf == pytest.approx(-174.65, abs=1.0)
    assert result.phase.value_inf_error is not None and result.phase.value_inf_error < 5.0
    # the correction from the finest mesh actually run (96 cells, -172.11 deg) is still > 0.5 deg
    assert not result.is_converged_within(0.5)


def test_not_monotonic_gives_no_fit_not_a_bad_fit():
    """Synthetic non-monotonic sequence: must refuse to fit, per 'fit nothing unless monotonic'."""
    fit = conv._fit_power_law([0.1, 0.05, 0.025, 0.0125], [1.0, 2.0, 1.5, 1.7])
    assert fit.status == conv.STATUS_NOT_MONOTONIC_NO_FIT
    assert fit.value_inf is None


def test_too_few_points():
    fit = conv._fit_power_law([0.1, 0.05], [1.0, 1.2])
    assert fit.status == conv.STATUS_TOO_FEW_POINTS


def test_power_law_fit_recovers_known_synthetic_law():
    dx = np.array([0.2, 0.1, 0.05, 0.025, 0.0125])
    true_inf, true_a, true_p = -170.0, 40.0, 1.0
    v = true_inf + true_a * dx**true_p
    fit = conv._fit_power_law(dx, v)
    assert fit.status == conv.STATUS_CONVERGED
    assert fit.value_inf == pytest.approx(true_inf, abs=1e-4)
    assert fit.order == pytest.approx(true_p, abs=1e-3)
    # value_inf_error is |v_inf - finest measured value| (PHASE13B: curve_fit's covariance is not
    # a credible uncertainty), not a fit-quality residual -- even for this exact synthetic law it
    # is the size of the last step, not ~0.
    finest_v = true_inf + true_a * dx[-1] ** true_p
    assert fit.value_inf_error == pytest.approx(abs(true_inf - finest_v), abs=1e-4)


@pytest.mark.slow
def test_off_resonance_case_converges_well_within_the_finest_fixed_mesh():
    """Off-resonance (L = 2.0 mm, the same centre Phase 12D used): first-order error should be
    much smaller than at resonance, since there is no near-singular edge current there."""
    uc = dp.unit_cell_at(2.0e-3)
    result = conv.richardson_extrapolate_s11xx(uc, 12.5e9, dp.REFERENCE_PLANES, cell_counts=conv.FAST_CELL_COUNTS)
    assert result.phase.status == conv.STATUS_CONVERGED
    assert result.is_converged_within(1.0)
    # PHASE13B: value_inf_error must be the extrapolation-correction size, defined even for the
    # app's 4-point FAST_CELL_COUNTS family (where a leave-one-out spread would not be).
    assert result.phase.value_inf_error == pytest.approx(abs(result.phase.value_inf - conv._unwrap_to_match(result)))
