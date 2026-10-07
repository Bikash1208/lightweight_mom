"""Integrated polarization metrics for a requested beam (XPD over the main-beam region, polarization
efficiency, co-pol side-lobe level excluding the main beam, cross-pol peak).

The designs use the ideal phase library of `test_hologram_vortex_sign` (every phase available), so
the checks measure the vector model and the polarization compensation, not a solved cell library."""

import numpy as np
import pytest

from src.dispersion.analytic_tm0 import tm0_beta
from src.green_function import LayeredGreenFunction
from src.hologram.hmsa_synthesis import synthesize_formulation_b
from src.materials import Substrate
from src.pattern.vector_aperture_pattern import (
    aperture_current_from_scalar,
    beam_polarization_metrics,
    compute_vector_pattern,
)
from src.phase8b.object_wave import Beam
from src.phase8b.reference_wave import monopole_reference_field
from src.phase8e.reference_field_input import phase_only_reference_field
from tests.test_hologram_vortex_sign import _ideal_library

FREQ = 10.54e9
PERIOD = 8e-3
SUB = Substrate("user", 2.2, 0.0009, 1.575e-3)
K_SWG = tm0_beta(FREQ, 2.2, 1.575e-3)


@pytest.fixture(scope="module")
def pencil_patterns():
    xs = (np.arange(37) - 18) * PERIOD + PERIOD / 2.0
    X, Y = np.meshgrid(xs, xs)
    reference = phase_only_reference_field(X, Y, np.degrees(np.angle(monopole_reference_field(X, Y, K_SWG))))
    green = LayeredGreenFunction(SUB, FREQ, backing="pec")
    patterns = {}
    for pol in ("hp", "rhcp"):
        beam = Beam(theta_deg=20.0, phi_deg=0.0, weight=1.0 + 0j, desired_polarization=pol, oam_charge=0)
        design = synthesize_formulation_b(reference, [beam], FREQ, _ideal_library(), PERIOD, feed_xy_m=(0.0, 0.0))
        dx = float(design.x_m[0, 1] - design.x_m[0, 0])
        jx, jy = aperture_current_from_scalar(design.x_m, design.y_m, design.aperture_field, feed_xy_m=(0.0, 0.0))
        patterns[pol] = compute_vector_pattern(design.x_m, design.y_m, jx, jy, FREQ, green, dx, dx)
    return patterns


def test_hp_pencil_integrated_xpd_is_below_10_db(pencil_patterns):
    metrics = beam_polarization_metrics(pencil_patterns["hp"], "hp", 20.0, 0.0)
    assert metrics.integrated_xpd_db < 10.0


def test_hp_pencil_polarization_efficiency_is_about_half(pencil_patterns):
    metrics = beam_polarization_metrics(pencil_patterns["hp"], "hp", 20.0, 0.0)
    assert metrics.polarization_efficiency == pytest.approx(0.5, abs=0.10)


def test_hp_pencil_point_xpd_at_the_beam_centre_stays_high(pencil_patterns):
    metrics = beam_polarization_metrics(pencil_patterns["hp"], "hp", 20.0, 0.0)
    assert metrics.point_xpd_db > 30.0


def test_cross_pol_peak_and_copol_side_lobe_are_reported_separately(pencil_patterns):
    metrics = beam_polarization_metrics(pencil_patterns["hp"], "hp", 20.0, 0.0)
    assert metrics.cross_peak_db is not None and metrics.sll_copol_db is not None
    assert metrics.region_fraction < 0.05


def test_rhcp_pencil_efficiency_is_about_half_too(pencil_patterns):
    metrics = beam_polarization_metrics(pencil_patterns["rhcp"], "rhcp", 20.0, 0.0)
    assert metrics.polarization_efficiency == pytest.approx(0.5, abs=0.10)


def test_rhcp_pencil_copol_sll_is_a_real_sidelobe_not_the_main_lobes_own_far_shoulder(pencil_patterns):
    """Regression (PHASE16 follow-ups): this beam is requested at phi=0, straddling this grid's
    own phi=0/360 seam. An earlier flood-fill main-lobe exclusion split the lobe into two
    disconnected islands there, reporting its own far shoulder (~-1.3 dB) as the side-lobe level
    instead of the true first sidelobe (well below -10 dB for a clean circularly-polarized pencil
    on this aperture)."""
    metrics = beam_polarization_metrics(pencil_patterns["rhcp"], "rhcp", 20.0, 0.0)
    assert metrics.sll_copol_db < -10.0
