"""Formulation A with several beams on one port: the summed object wave is normalised to unit peak
amplitude (`hmsa_synthesis.object_wave_scale`), so Psi stays in [-1, 1].

Beam-level balance is not met with equal weights: the per-beam levels span about 6.4 dB, from the
real-valued modulation, not from saturation. That check is a strict expected failure until the
weights are balanced (a decision, since balancing changes the weights the user set)."""

import numpy as np
import pytest

from src.dispersion.hfss_table import DispersionSample, DispersionTable
from src.green_function import LayeredGreenFunction
from src.hologram import hmsa_synthesis as hs
from src.materials import Substrate
from src.pattern.vector_aperture_pattern import aperture_current_from_scalar, compute_vector_pattern
from src.phase8b.object_wave import Beam
from src.phase8e.reference_field_input import spherical_feed_reference_field

FREQ = 10.54e9
PERIOD = 4e-3
BEAMS = [(15.0, 0.0), (30.0, 90.0), (45.0, 180.0), (60.0, 270.0)]
SUB = Substrate("user", 2.2, 0.0009, 1.575e-3)


def _table():
    samples = tuple(
        DispersionSample(geometry_param=g, beta_rad_per_m=None, x_s_ohm=x)
        for g, x in [(1.0, 50.0), (2.0, 100.0), (3.0, 150.0), (4.0, 200.0)]
    )
    return DispersionTable(source_path="synthetic", frequency_hz=FREQ, samples=samples)


def _four_beam_design():
    c = np.linspace(-0.06, 0.06, 31)
    x, y = np.meshgrid(c, c, indexing="xy")
    ref = spherical_feed_reference_field(x, y, feed_height_m=0.05, frequency_hz=FREQ, amplitude_at_1m=0.05)
    beams = [Beam(theta_deg=t, phi_deg=p, weight=1.0 + 0j, desired_polarization="rhcp") for t, p in BEAMS]
    return hs.synthesize_formulation_a(
        ref, beams, FREQ, x_avg_ohm=125.0, modulation_index=0.3, dispersion_table=_table(), period_m=PERIOD
    )


def _local_peak_levels_db(design):
    jx, jy = aperture_current_from_scalar(design.x_m, design.y_m, design.aperture_field, feed_xy_m=(0.0, 0.0))
    green = LayeredGreenFunction(SUB, FREQ, backing="pec")
    pat = compute_vector_pattern(design.x_m, design.y_m, jx, jy, FREQ, green, PERIOD, PERIOD, 60.0, 121, 73)
    power = np.abs(pat.e_theta) ** 2 + np.abs(pat.e_phi) ** 2
    peak = float(np.max(power))
    levels = []
    for t, p in BEAMS:
        t0, p0 = np.radians(t), np.radians(p)
        tt, pp = np.radians(pat.theta_deg), np.radians(pat.phi_deg)
        cos_ang = np.cos(tt) * np.cos(t0) + np.sin(tt) * np.sin(t0) * np.cos(pp - p0)
        near = np.degrees(np.arccos(np.clip(cos_ang, -1.0, 1.0))) <= 12.0
        levels.append(10.0 * np.log10(np.max(np.where(near, power, 0.0)) / peak))
    return levels


def test_single_unit_beam_is_not_rescaled():
    c = np.linspace(-0.02, 0.02, 5)
    x, y = np.meshgrid(c, c, indexing="xy")
    e_obj = np.exp(1j * np.zeros_like(x))
    assert hs.object_wave_scale(e_obj) == 1.0


def test_four_beam_single_port_unreachable_is_below_5_percent():
    design = _four_beam_design()
    assert design.unreachable_fraction < 0.05
    assert any(label.startswith("OBJECT_WAVE_NORMALISED") for label in design.labels)


@pytest.mark.xfail(strict=True, reason="equal weights give a 6.4 dB spread across four beams; balancing is a decision")
def test_four_beam_single_port_beams_are_within_3_db():
    levels = _local_peak_levels_db(_four_beam_design())
    assert max(levels) - min(levels) <= 3.0
