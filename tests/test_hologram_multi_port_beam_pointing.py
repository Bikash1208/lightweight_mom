"""A 1x3 RHCP board at 25 deg (10.54 GHz, 4 mm cells, 28 mm guard gap, 120 mm active ports): each
port's own main beam against the request. The ideal phase library (`test_hologram_vortex_sign`)
isolates the synthesis and the vector model from library quantisation."""

import numpy as np
import pytest

from src.dispersion.analytic_tm0 import tm0_beta
from src.green_function import LayeredGreenFunction
from src.hologram.hmsa_synthesis import synthesize_formulation_b
from src.hologram.multi_port_synthesis import port_reference_field, synthesize_ports
from src.hologram.ports import layout_preset_1xn
from src.materials import Substrate
from src.pattern.vector_aperture_pattern import aperture_current_from_scalar, compute_vector_pattern
from src.phase8b.object_wave import Beam
from tests.test_hologram_vortex_sign import _ideal_library

FREQ = 10.54e9
CELL = 4e-3
GAP = 28e-3
REQUEST_THETA_DEG = 25.0
SUB = Substrate("user", 2.2, 0.0009, 1.575e-3)
K_SWG = tm0_beta(FREQ, 2.2, 1.575e-3)


def _worst_pointing_error_deg(active_side_m):
    regions = layout_preset_1xn(3, active_side_m, gap_m=GAP)
    beams = {
        r.port_id: [Beam(theta_deg=REQUEST_THETA_DEG, phi_deg=0.0, weight=1.0 + 0j, desired_polarization="rhcp")]
        for r in regions
    }
    library = _ideal_library()
    multi = synthesize_ports(
        regions, beams,
        lambda ref, bs, feed: synthesize_formulation_b(ref, bs, FREQ, library, CELL, feed_xy_m=feed),
        CELL, k_swg=K_SWG, min_lambda_sw=1.0,
    )
    green = LayeredGreenFunction(SUB, FREQ, backing="pec")
    worst = 0.0
    for port in multi.ports:
        design = port.design
        jx, jy = aperture_current_from_scalar(design.x_m, design.y_m, design.aperture_field, feed_xy_m=port.feed_xy_m)
        pat = compute_vector_pattern(design.x_m, design.y_m, jx, jy, FREQ, green, CELL, CELL, 60.0, 241, 73)
        worst = max(worst, abs(pat.main_beam_theta_deg - REQUEST_THETA_DEG))
    return worst


def test_1x3_rhcp_board_at_25_deg_points_within_1_5_deg():
    assert _worst_pointing_error_deg(120e-3) <= 1.5


def test_corrected_beam_is_aimed_past_the_measured_miss():
    from src.hologram.pointing_correction import BeamAim, corrected_beams

    beam = Beam(theta_deg=REQUEST_THETA_DEG, phi_deg=0.0, weight=1.0 + 0j, desired_polarization="rhcp")
    aim = BeamAim(REQUEST_THETA_DEG, 0.0, 23.25, 0.0, 0.0)
    (aimed,) = corrected_beams([beam], [aim], [beam])
    # re-aim in direction-cosine space: sin(theta_new) = 2 sin(target) - sin(achieved)
    expected = float(np.degrees(np.arcsin(2.0 * np.sin(np.radians(25.0)) - np.sin(np.radians(23.25)))))
    assert aimed.theta_deg == pytest.approx(expected, abs=1e-6)
    assert aimed.phi_deg == pytest.approx(0.0)


def test_pointing_correction_converges_for_a_steep_beam_not_on_the_phi_axis():
    """Regression for a live-browser stall (PHASE16 follow-ups): a beam at (35deg, 180deg) only
    moved from a 1.43deg miss to a 1.41deg miss over a pass, because the shared hemisphere grid's
    own phi spacing -- several deg once projected to a real angle at this theta -- was coarser than
    the miss pointing correction was trying to resolve. The fine local measurement fixes this."""
    from src.hologram.pointing_correction import correct_beam_aim, measure_beam_aims

    regions = layout_preset_1xn(1, 100e-3)
    region = regions[0]
    library = _ideal_library()
    green = LayeredGreenFunction(SUB, FREQ, backing="pec")
    beams = [Beam(theta_deg=35.0, phi_deg=180.0, weight=1.0 + 0j, desired_polarization="rhcp")]
    feed = (region.feed_x_m, region.feed_y_m)

    reference = port_reference_field(region, CELL / 4.0, K_SWG)

    def synthesize(bs):
        return synthesize_formulation_b(reference, bs, FREQ, library, CELL, feed_xy_m=feed)

    def measure(design, bs):
        return measure_beam_aims(design, bs, green, CELL, FREQ, feed)

    _design, _used, history = correct_beam_aim(beams, synthesize, measure, iterations=2)
    assert history[0][0].miss_deg > 1.0
    assert history[-1][0].miss_deg < 0.3


def test_pointing_correction_brings_the_92_mm_port_onto_its_request():
    from src.hologram.pointing_correction import correct_beam_aim, measure_beam_aims

    regions = layout_preset_1xn(1, 92e-3)
    region = regions[0]
    library = _ideal_library()
    green = LayeredGreenFunction(SUB, FREQ, backing="pec")
    beams = [Beam(theta_deg=REQUEST_THETA_DEG, phi_deg=0.0, weight=1.0 + 0j, desired_polarization="rhcp")]
    feed = (region.feed_x_m, region.feed_y_m)

    reference = port_reference_field(region, CELL / 4.0, K_SWG)

    def synthesize(bs):
        return synthesize_formulation_b(reference, bs, FREQ, library, CELL, feed_xy_m=feed)

    def measure(design, bs):
        return measure_beam_aims(design, bs, green, CELL, FREQ, feed)

    _design, _used, history = correct_beam_aim(beams, synthesize, measure, iterations=2)
    assert history[0][0].miss_deg > 1.0
    assert history[-1][0].miss_deg < 0.5
