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


def test_formulation_b_aperture_field_reconstructs_the_object_beam_when_unquantized():
    """Regression test: an earlier version of this code built the Pattern-prediction aperture
    field as `exp(j*achieved_modulation_phase)` alone, silently dropping the reference wave's own
    phase -- correct only for the special case of a UNIFORM-ZERO reference (every other test in
    this file happens to use exactly that degenerate case, which is why none of them caught it).
    With a genuinely non-uniform reference phase and a library fine enough to match the required
    modulation almost exactly (near-zero quantization error), the RADIATED aperture field
    (design.aperture_field = E_ref * modulation) must reconstruct the pure object-beam phase
    pattern -- exp(j*obj_phase) -- not the modulation phasor alone."""
    x, y = _grid(n=7, extent_m=0.03)
    # a non-trivial, non-uniform reference phase (a linear ramp, not the degenerate uniform-zero
    # case every other test in this file uses)
    ref_phase_deg = 40.0 * x / 0.01
    ref = phase_only_reference_field(x, y, ref_phase_deg)
    beams = [Beam(theta_deg=15.0, phi_deg=0.0)]

    # A dense, fine-grained phase library (360 candidates, ~1deg spacing) so the achieved phase is
    # essentially exact -- isolates the aperture-field construction from quantization error.
    candidates = [(f"c{i}", 1.0, 1.0, "lib", np.exp(1j * np.radians(i))) for i in range(-180, 180)]
    design = hs.synthesize_formulation_b(ref, beams, FREQ, candidates, period_m=8e-3)

    from src.phase8b.object_wave import combine_beams_scalar

    e_obj = combine_beams_scalar(design.x_m, design.y_m, beams, FREQ)
    expected = np.exp(1j * np.angle(e_obj))  # the pure desired-beam phase pattern

    achieved_deg = np.array([c.achieved_value for c in design.cells]).reshape(design.x_m.shape)
    # sanity: quantization really is near-exact with this dense a library
    assert np.max([c.phase_error_deg for c in design.cells]) < 2.0

    phase_diff = np.angle(design.aperture_field * np.conj(expected))
    assert np.max(np.abs(np.degrees(phase_diff))) < 3.0  # matches to within quantization error

    # the OLD (buggy) construction -- the modulation phasor alone -- must NOT match nearly as well,
    # proving this test actually distinguishes the fix from the regression it guards against.
    old_aperture_field = np.exp(1j * np.radians(achieved_deg))
    old_phase_diff = np.angle(old_aperture_field * np.conj(expected))
    assert np.max(np.abs(np.degrees(old_phase_diff))) > 30.0


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


def test_formulation_b_natural_pencil_default_is_bit_identical_to_pre_phase15():
    """PHASE15 Part A's compatibility invariant: leaving every beam at its new defaults
    (`desired_polarization="natural"`, `oam_charge=0`) and `multiplexing="full_aperture"` must
    reproduce exactly what `synthesize_formulation_b` computed before Part A existed."""
    x, y = _grid()
    ref = phase_only_reference_field(x, y, np.zeros_like(x))
    beams = [Beam(theta_deg=15.0, phi_deg=30.0)]
    candidates = [("c1", 2.0, 2.0, "lib", complex(1.0, 0.0))]
    design_default = hs.synthesize_formulation_b(ref, beams, FREQ, candidates, period_m=8e-3)
    design_explicit = hs.synthesize_formulation_b(
        ref,
        beams,
        FREQ,
        candidates,
        period_m=8e-3,
        feed_xy_m=(0.0, 0.0),
        cp_phase_offset_deg=0.0,
        multiplexing="full_aperture",
    )
    assert np.array_equal(design_default.control_map, design_explicit.control_map)
    assert np.array_equal(design_default.aperture_field, design_explicit.aperture_field)


def test_formulation_b_polarization_compensation_changes_the_control_map():
    x, y = _grid()
    ref = phase_only_reference_field(x, y, np.zeros_like(x))
    natural = [Beam(theta_deg=15.0, phi_deg=0.0, desired_polarization="natural")]
    rhcp = [Beam(theta_deg=15.0, phi_deg=0.0, desired_polarization="rhcp")]
    candidates = [(f"c{i}", 1.0, 1.0, "lib", np.exp(1j * np.radians(i))) for i in range(-180, 180)]
    d_natural = hs.synthesize_formulation_b(ref, natural, FREQ, candidates, period_m=8e-3)
    d_rhcp = hs.synthesize_formulation_b(ref, rhcp, FREQ, candidates, period_m=8e-3)
    assert not np.allclose(d_natural.control_map, d_rhcp.control_map)


def test_formulation_b_sector_multiplexing_restricts_each_beam_to_its_own_region():
    x, y = _grid(n=9, extent_m=0.03)
    ref = phase_only_reference_field(x, y, np.zeros_like(x))
    beams = [Beam(theta_deg=15.0, phi_deg=0.0), Beam(theta_deg=15.0, phi_deg=180.0)]
    candidates = [(f"c{i}", 1.0, 1.0, "lib", np.exp(1j * np.radians(i))) for i in range(-180, 180)]
    full = hs.synthesize_formulation_b(ref, beams, FREQ, candidates, period_m=8e-3, multiplexing="full_aperture")
    sectors = hs.synthesize_formulation_b(ref, beams, FREQ, candidates, period_m=8e-3, multiplexing="sectors")
    assert not np.allclose(full.control_map, sectors.control_map)


def test_synthesize_formulation_a_accepts_feed_and_multiplexing_kwargs():
    x, y = _grid()
    ref = spherical_feed_reference_field(x, y, feed_height_m=0.05, frequency_hz=FREQ, amplitude_at_1m=0.05)
    beams = [Beam(theta_deg=0.0, phi_deg=0.0, desired_polarization="hp")]
    table = _dispersion_table()
    design = hs.synthesize_formulation_a(
        ref,
        beams,
        FREQ,
        x_avg_ohm=125.0,
        modulation_index=0.3,
        dispersion_table=table,
        period_m=8e-3,
        feed_xy_m=(0.0, 0.0),
        cp_phase_offset_deg=0.0,
        multiplexing="full_aperture",
    )
    assert design.formulation == "A_impedance_modulation"


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
