"""Phase 15 Part C validation (Mishra, Barad, Sarkar, IEEE TVT 2026): reproduce the paper's own
5-port layout QUALITATIVELY -- four single-beam ports at (15,0),(30,90),(45,180),(60,270) deg and
one quad-beam port radiating all four, with LHCP from this project's own vector model -- on a small
synthetic board (not the paper's own 26 GHz dimensions: this test only needs the RATIO of port size
to wavelength and to unit-cell period the paper itself uses, not its absolute frequency; this
project's hologram synthesis is scale-invariant in that ratio). No gain-agreement claim is made
(this project has no measurement of the paper's own unit cell)."""

import math

import numpy as np
import pytest

from src.dispersion.analytic_tm0 import tm0_beta
from src.frequency import free_space_wavelength
from src.green_function import LayeredGreenFunction
from src.hologram import hmsa_synthesis as hs
from src.hologram.multi_port_synthesis import synthesize_ports
from src.hologram.ports import layout_preset_2x2_plus_1_large
from src.materials import Substrate
from src.pattern.vector_aperture_pattern import (
    aperture_current_from_scalar,
    compute_vector_pattern,
    copol_component,
)
from src.phase8b.object_wave import Beam

FREQ_HZ = 10e9
SUBSTRATE = Substrate("rt5880-like", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
K_SWG = tm0_beta(FREQ_HZ, SUBSTRATE.eps_r, SUBSTRATE.thickness)
WAVELENGTH_M = free_space_wavelength(FREQ_HZ)
PERIOD_M = WAVELENGTH_M / 3.0
SMALL_SIDE_M = 5.0 * WAVELENGTH_M

_BEAM_DIRECTIONS = [(15.0, 0.0), (30.0, 90.0), (45.0, 180.0), (60.0, 270.0)]
_CANDIDATES = [(f"c{i}", 2.0, 2.0, "lib", np.exp(1j * math.radians(i * 5.0))) for i in range(72)]


def _per_port(ref, beams, feed_xy_m):
    return hs.synthesize_formulation_b(ref, beams, FREQ_HZ, _CANDIDATES, period_m=PERIOD_M, feed_xy_m=feed_xy_m)


def _green_function():
    return LayeredGreenFunction(SUBSTRATE, FREQ_HZ, backing="pec")


def _peak_direction(design, feed_xy_m):
    jx, jy = aperture_current_from_scalar(design.x_m, design.y_m, design.aperture_field, feed_xy_m=feed_xy_m)
    dx = float(design.x_m[0, 1] - design.x_m[0, 0])
    dy = float(design.y_m[1, 0] - design.y_m[0, 0])
    pat = compute_vector_pattern(
        design.x_m, design.y_m, jx, jy, FREQ_HZ, _green_function(), dx, dy, theta_max_deg=70.0, n_theta=36, n_phi=73
    )
    return pat


@pytest.mark.slow
def test_mishra_5port_layout_steers_each_port_to_its_own_commanded_direction():
    regions = layout_preset_2x2_plus_1_large(SMALL_SIDE_M)
    beams = {
        i: [Beam(theta_deg=t, phi_deg=p, desired_polarization="lhcp")] for i, (t, p) in enumerate(_BEAM_DIRECTIONS)
    }
    beams[4] = [Beam(theta_deg=t, phi_deg=p, desired_polarization="lhcp", weight=0.5) for t, p in _BEAM_DIRECTIONS]
    multi = synthesize_ports(regions, beams, _per_port, period_m=PERIOD_M, k_swg=K_SWG)

    for port_result, (theta_cmd, phi_cmd) in zip(multi.ports[:4], _BEAM_DIRECTIONS, strict=True):
        pat = _peak_direction(port_result.design, port_result.feed_xy_m)
        assert pat.main_beam_theta_deg == pytest.approx(theta_cmd, abs=8.0)
        phi_err = abs((pat.main_beam_phi_deg - phi_cmd + 180.0) % 360.0 - 180.0)
        assert phi_err < 10.0, f"port {port_result.port_id}: expected phi~{phi_cmd}, got {pat.main_beam_phi_deg}"

        # LHCP from the vector model: the beam's own co-pol (LHCP) component should carry
        # materially more power than its cross-pol (RHCP) component at the commanded direction.
        i_t = np.argmin(np.abs(pat.theta_deg[:, 0] - theta_cmd))
        j_p = np.argmin(np.abs((pat.phi_deg[0, :] - phi_cmd + 180.0) % 360.0 - 180.0))
        co = copol_component(np.array([pat.e_theta[i_t, j_p]]), np.array([pat.e_phi[i_t, j_p]]), "lhcp")[0]
        cross = copol_component(np.array([pat.e_theta[i_t, j_p]]), np.array([pat.e_phi[i_t, j_p]]), "rhcp")[0]
        xpd_db = 10.0 * np.log10(np.abs(co) ** 2 / np.abs(cross) ** 2)
        assert xpd_db > 3.0, f"port {port_result.port_id}: LHCP/RHCP ratio only {xpd_db:.1f}dB"


@pytest.mark.slow
def test_mishra_quad_beam_port_radiates_near_all_four_commanded_directions():
    regions = layout_preset_2x2_plus_1_large(SMALL_SIDE_M)
    beams = {
        i: [Beam(theta_deg=t, phi_deg=p, desired_polarization="lhcp")] for i, (t, p) in enumerate(_BEAM_DIRECTIONS)
    }
    beams[4] = [Beam(theta_deg=t, phi_deg=p, desired_polarization="lhcp", weight=0.5) for t, p in _BEAM_DIRECTIONS]
    multi = synthesize_ports(regions, beams, _per_port, period_m=PERIOD_M, k_swg=K_SWG)
    quad_port = multi.ports[4]
    pat = _peak_direction(quad_port.design, quad_port.feed_xy_m)

    power = pat.power_db
    for theta_cmd, phi_cmd in _BEAM_DIRECTIONS:
        i_t = np.argmin(np.abs(pat.theta_deg[:, 0] - theta_cmd))
        j_p = np.argmin(np.abs((pat.phi_deg[0, :] - phi_cmd + 180.0) % 360.0 - 180.0))
        # each commanded direction should be a local-ish high-power region of the quad-beam
        # pattern (within 10dB of the pattern's own overall peak) -- "radiates all four beams",
        # qualitatively, not a claim of four equal-height, perfectly isolated main lobes.
        assert power[i_t, j_p] > -10.0, f"({theta_cmd},{phi_cmd}): only {power[i_t, j_p]:.1f}dB below peak"
