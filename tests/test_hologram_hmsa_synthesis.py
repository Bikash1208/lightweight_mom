"""Part 3: HMSA hologram synthesis, both formulations, end to end on small synthetic fields."""

import numpy as np
import pytest

from src.dispersion.hfss_table import DispersionSample, DispersionTable
from src.hologram import hmsa_synthesis as hs
from src.phase8b.object_wave import Beam
from src.phase8e.reference_field_input import phase_only_reference_field, spherical_feed_reference_field

FREQ = 10e9


def _grid(n=5, extent_m=0.02):
    coords = np.linspace(-extent_m, extent_m, n)
    return np.meshgrid(coords, coords, indexing="xy")


def _dispersion_table():
    samples = tuple(DispersionSample(geometry_param=g, beta_rad_per_m=None, x_s_ohm=x)
                     for g, x in [(1.0, 50.0), (2.0, 100.0), (3.0, 150.0), (4.0, 200.0)])
    return DispersionTable(source_path="synthetic", frequency_hz=FREQ, samples=samples)


def test_formulation_b_phase_only_field_end_to_end():
    x, y = _grid()
    ref = phase_only_reference_field(x, y, np.zeros_like(x))   # uniform-phase reference
    beams = [Beam(theta_deg=0.0, phi_deg=0.0)]   # broadside object beam
    candidates = [
        ("c1", 2.0, 2.0, "lib", complex(1.0, 0.0)),
        ("c2", 3.0, 3.0, "lib", np.exp(1j * np.radians(45.0))),
        ("c3", 4.0, 4.0, "lib", np.exp(1j * np.radians(-45.0))),
    ]
    design = hs.synthesize_formulation_b(ref, beams, FREQ, candidates, period_m=8e-3)
    assert design.formulation == "B_phase_difference"
    assert "LOCAL_PERIODICITY_APPROXIMATION" in design.labels
    assert len(design.cells) == x.size
    assert all(c.reachable for c in design.cells)


def test_formulation_a_needs_complex_field():
    x, y = _grid()
    ref = phase_only_reference_field(x, y, np.zeros_like(x))
    beams = [Beam(theta_deg=0.0, phi_deg=0.0)]
    with pytest.raises(ValueError, match="complex reference field"):
        hs.formulation_a_control_map(ref, beams, FREQ, x_avg_ohm=100.0, modulation_index=0.5)


def test_formulation_a_end_to_end_with_complex_field_and_dispersion_table():
    x, y = _grid()
    # `interference_modulation`'s Psi is dimensionless only for UNIT-amplitude fields (its own
    # docstring: "performs no implicit normalization"); amplitude_at_1m is chosen so |E_ref| ~ 1
    # at these small radii (feed_height=0.05m), so Psi stays in the physically-intended [-1, 1].
    ref = spherical_feed_reference_field(x, y, feed_height_m=0.05, frequency_hz=FREQ, amplitude_at_1m=0.05)
    beams = [Beam(theta_deg=0.0, phi_deg=0.0)]
    table = _dispersion_table()
    design = hs.synthesize_formulation_a(ref, beams, FREQ, x_avg_ohm=125.0, modulation_index=0.3,
                                          dispersion_table=table, period_m=8e-3)
    assert design.formulation == "A_impedance_modulation"
    assert "EXPERIMENTAL" in design.labels
    assert len(design.cells) == x.size
    # X_avg=125, M=0.3 -> X in [125*0.7, 125*1.3] = [87.5, 162.5], well within table's [50,200] range
    assert design.unreachable_fraction == 0.0


def test_out_of_range_x_s_marks_cells_unreachable():
    x, y = _grid()
    ref = spherical_feed_reference_field(x, y, feed_height_m=0.05, frequency_hz=FREQ)
    beams = [Beam(theta_deg=0.0, phi_deg=0.0)]
    table = _dispersion_table()   # range [50, 200] ohm
    design = hs.synthesize_formulation_a(ref, beams, FREQ, x_avg_ohm=1000.0, modulation_index=0.1,
                                          dispersion_table=table, period_m=8e-3)
    assert design.unreachable_fraction > 0.0
    assert any(not c.reachable for c in design.cells)
