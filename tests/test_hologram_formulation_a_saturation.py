"""Formulation A with several beams on one port: the summed object wave is normalised to unit peak
amplitude (`hmsa_synthesis.object_wave_scale`), so Psi stays in [-1, 1].

With equal requested levels the raw per-beam levels span about 6.4 dB (the real-valued modulation,
not saturation). The default level matching brings them within 3 dB."""

import numpy as np

from src.dispersion.hfss_table import DispersionSample, DispersionTable
from src.green_function import LayeredGreenFunction
from src.hologram import hmsa_synthesis as hs
from src.hologram.pointing_correction import measure_beam_aims, refine_design
from src.materials import Substrate
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


def _synthesize_for(ref, table):
    def synthesize(beams):
        return hs.synthesize_formulation_a(
            ref, beams, FREQ, x_avg_ohm=125.0, modulation_index=0.3, dispersion_table=table, period_m=PERIOD
        )

    return synthesize


def _reference():
    c = np.linspace(-0.06, 0.06, 31)
    x, y = np.meshgrid(c, c, indexing="xy")
    return spherical_feed_reference_field(x, y, feed_height_m=0.05, frequency_hz=FREQ, amplitude_at_1m=0.05)


def _four_beam_beams():
    return [Beam(theta_deg=t, phi_deg=p, weight=1.0 + 0j, desired_polarization="rhcp") for t, p in BEAMS]


def _four_beam_design():
    return _synthesize_for(_reference(), _table())(_four_beam_beams())


def _level_spread_db(design, beams):
    green = LayeredGreenFunction(SUB, FREQ, backing="pec")
    aims = measure_beam_aims(design, beams, green, PERIOD, FREQ, (0.0, 0.0))
    levels = np.array([a.level_db for a in aims])
    return float(levels.max() - levels.min())


def test_four_beam_single_port_unreachable_is_below_5_percent():
    design = _four_beam_design()
    assert design.unreachable_fraction < 0.05
    assert any(label.startswith("OBJECT_WAVE_NORMALISED") for label in design.labels)


def test_four_beam_single_port_beams_are_within_3_db_after_level_matching():
    beams = _four_beam_beams()
    green = LayeredGreenFunction(SUB, FREQ, backing="pec")
    synthesize = _synthesize_for(_reference(), _table())

    def measure(design, beams_now):
        return measure_beam_aims(design, beams_now, green, PERIOD, FREQ, (0.0, 0.0))

    design, _used, report = refine_design(beams, synthesize, measure, pointing=False, balance=True)
    assert report["balance"]["max_error_after_db"] <= report["balance"]["max_error_before_db"]
    assert _level_spread_db(design, _used) <= 3.0


def test_single_unit_beam_is_not_rescaled():
    c = np.linspace(-0.02, 0.02, 5)
    x, y = np.meshgrid(c, c, indexing="xy")
    e_obj = np.exp(1j * np.zeros_like(x))
    assert hs.object_wave_scale(e_obj) == 1.0
