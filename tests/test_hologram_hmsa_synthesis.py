"""Part 3: HMSA hologram synthesis, both formulations, end to end on small synthetic fields."""

from dataclasses import replace

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
    samples = tuple(
        DispersionSample(geometry_param=g, beta_rad_per_m=None, x_s_ohm=x)
        for g, x in [(1.0, 50.0), (2.0, 100.0), (3.0, 150.0), (4.0, 200.0)]
    )
    return DispersionTable(source_path="synthetic", frequency_hz=FREQ, samples=samples)


def test_formulation_b_phase_only_field_end_to_end():
    x, y = _grid()
    ref = phase_only_reference_field(x, y, np.zeros_like(x))  # uniform-phase reference
    beams = [Beam(theta_deg=0.0, phi_deg=0.0)]  # broadside object beam
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
    design = hs.synthesize_formulation_a(
        ref, beams, FREQ, x_avg_ohm=125.0, modulation_index=0.3, dispersion_table=table, period_m=8e-3
    )
    assert design.formulation == "A_impedance_modulation"
    assert "EXPERIMENTAL" in design.labels
    assert len(design.cells) == x.size
    # X_avg=125, M=0.3 -> X in [125*0.7, 125*1.3] = [87.5, 162.5], well within table's [50,200] range
    assert design.unreachable_fraction == 0.0


def test_out_of_range_x_s_marks_cells_unreachable():
    x, y = _grid()
    ref = spherical_feed_reference_field(x, y, feed_height_m=0.05, frequency_hz=FREQ)
    beams = [Beam(theta_deg=0.0, phi_deg=0.0)]
    table = _dispersion_table()  # range [50, 200] ohm
    design = hs.synthesize_formulation_a(
        ref, beams, FREQ, x_avg_ohm=1000.0, modulation_index=0.1, dispersion_table=table, period_m=8e-3
    )
    assert design.unreachable_fraction > 0.0
    assert any(not c.reachable for c in design.cells)


def test_synthesis_resamples_onto_the_cell_grid_not_the_raw_field_grid():
    """PHASE13B A2: a fine raw field grid (11x11, 2mm pitch, 20mm aperture) synthesized at an
    8mm period must produce cells for the PERIOD grid (2x2 cells fit in 20mm), not 121 (one per
    raw sample)."""
    n = 11
    extent_m = 0.01  # +-10mm -> 20mm aperture
    x, y = _grid(n=n, extent_m=extent_m)
    ref = phase_only_reference_field(x, y, np.zeros_like(x))
    ref = replace(ref, dx_m=0.002, dy_m=0.002)  # a real .fld import always carries its own pitch
    beams = [Beam(theta_deg=0.0, phi_deg=0.0)]
    candidates = [("c1", 2.0, 2.0, "lib", complex(1.0, 0.0))]
    design = hs.synthesize_formulation_b(ref, beams, FREQ, candidates, period_m=8e-3)
    assert design.cells != ()
    assert len(design.cells) == 4  # floor(20mm / 8mm) = 2 per side -> 2x2
    assert design.x_m.shape == (2, 2)


def test_phase_error_stats_and_shortfall_label_for_formulation_b():
    x, y = _grid()
    ref = phase_only_reference_field(x, y, np.zeros_like(x))
    beams = [Beam(theta_deg=0.0, phi_deg=0.0)]
    # only a 0deg candidate: broadside (uniform required phase 0deg) is an exact match everywhere,
    # so force a mismatch with a non-broadside object beam instead.
    beams = [Beam(theta_deg=60.0, phi_deg=0.0)]
    candidates = [("c1", 2.0, 2.0, "lib", complex(1.0, 0.0))]
    design = hs.synthesize_formulation_b(ref, beams, FREQ, candidates, period_m=8e-3)
    assert design.phase_error_stats is not None
    assert design.phase_error_stats["n_cells"] == len(design.cells)
    assert design.phase_error_stats["max_deg"] >= design.phase_error_stats["mean_deg"] >= 0.0
    if design.phase_error_stats["n_unreachable"] > 0:
        assert any(label.startswith("PHASE_COVERAGE_SHORTFALL(") for label in design.labels)
        assert any(not c.reachable for c in design.cells)
