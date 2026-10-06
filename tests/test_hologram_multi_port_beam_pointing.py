"""A 1x3 RHCP board at 25 deg (10.54 GHz, 4 mm cells, 28 mm guard gap): each port's own main beam
against the request. The ideal phase library (`test_hologram_vortex_sign`) isolates the synthesis
and the vector model from library quantisation.

The app's 1xN preset takes the guard gap out of each nominal side, so a 142 mm nominal port (the
new default, 5 lambda0) keeps 114 mm of aperture; a 120 mm nominal port keeps 92 mm. The finite
aperture pulls the beam inward, so the 92 mm case misses the 1.5 deg target (measured 1.75 deg).
That case is kept as a strict expected failure rather than dropped."""

import pytest

from src.dispersion.analytic_tm0 import tm0_beta
from src.green_function import LayeredGreenFunction
from src.hologram.hmsa_synthesis import synthesize_formulation_b
from src.hologram.multi_port_synthesis import synthesize_ports
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


def _worst_pointing_error_deg(nominal_side_m):
    regions = layout_preset_1xn(3, nominal_side_m, gap_m=GAP)
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


def test_default_port_size_points_the_rhcp_beam_within_1_5_deg():
    assert _worst_pointing_error_deg(142e-3) <= 1.5


@pytest.mark.xfail(strict=True, reason="120 mm nominal ports keep 92 mm of aperture; measured 1.75 deg pointing error")
def test_120_mm_ports_point_the_rhcp_beam_within_1_5_deg():
    assert _worst_pointing_error_deg(120e-3) <= 1.5
