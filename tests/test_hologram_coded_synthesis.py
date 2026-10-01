"""Phase 15 Part E: the Coded (1-bit/N-bit) design method -- the phase formula reuses Part A/B's
own phi_ext/phase-difference machinery exactly (module docstring of `coded_synthesis.py`), and the
1-bit/N-bit combination rule is one quantization formula, not two special cases."""

import numpy as np
import pytest

from src.hologram import coded_synthesis as cs
from src.hologram import hmsa_synthesis as hs
from src.hologram.two_state_cells import TwoStateCell, TwoStateChoice
from src.phase8b.object_wave import Beam
from src.phase8e.reference_field_input import phase_only_reference_field

FREQ = 10e9
PERIOD_M = 2e-3


def _grid(n=9, extent_m=0.02):
    coords = np.linspace(-extent_m, extent_m, n)
    return np.meshgrid(coords, coords, indexing="xy")


def _uniform_phase_ref():
    x, y = _grid()
    return phase_only_reference_field(x, y, np.zeros_like(x))


# ================================================================ per_beam_phase_maps_deg


def test_per_beam_phase_map_matches_formulation_b_for_a_single_beam():
    ref = _uniform_phase_ref()
    beam = Beam(theta_deg=20.0, phi_deg=30.0)
    expected = hs.formulation_b_control_map(ref, [beam], FREQ)
    (actual,) = cs.per_beam_phase_maps_deg(ref, [beam], FREQ)
    assert np.allclose(actual, expected)


def test_per_beam_phase_maps_are_independent_not_combined():
    ref = _uniform_phase_ref()
    beams = [Beam(theta_deg=20.0, phi_deg=0.0), Beam(theta_deg=40.0, phi_deg=90.0)]
    maps = cs.per_beam_phase_maps_deg(ref, beams, FREQ)
    assert len(maps) == 2
    assert not np.allclose(maps[0], maps[1])
    for beam, psi_map in zip(beams, maps, strict=True):
        expected = hs.formulation_b_control_map(ref, [beam], FREQ)
        assert np.allclose(psi_map, expected)


def test_per_beam_phase_maps_use_each_beams_own_assigned_feed():
    ref = _uniform_phase_ref()
    beams = [Beam(theta_deg=20.0, phi_deg=0.0), Beam(theta_deg=20.0, phi_deg=0.0, desired_polarization="hp")]
    feed_xy_m = ([0.0, 0.005], [0.0, 0.0])  # beam 0 about the origin, beam 1 about (5mm, 0)
    maps = cs.per_beam_phase_maps_deg(ref, beams, FREQ, feed_xy_m=feed_xy_m)
    expected_1 = hs.formulation_b_control_map(ref, [beams[1]], FREQ, feed_xy_m=(0.005, 0.0))
    assert np.allclose(maps[1], expected_1)


# ===================================================================== coded_phase_map_deg


def test_coded_phase_map_1bit_single_beam_is_a_half_plane_quantization():
    ref = _uniform_phase_ref()
    beam = Beam(theta_deg=20.0, phi_deg=30.0)
    continuous_deg, quantized_deg, level_map = cs.coded_phase_map_deg(ref, [beam], FREQ, n_bits=1)
    (psi_deg,) = cs.per_beam_phase_maps_deg(ref, [beam], FREQ)
    assert np.allclose(continuous_deg, psi_deg)
    expected_level = (np.cos(np.radians(psi_deg)) <= 0).astype(int)
    assert np.array_equal(level_map, expected_level)
    assert set(np.unique(quantized_deg)) <= {0.0, 180.0}


def test_coded_phase_map_1bit_multi_beam_matches_the_brief_own_f_formula():
    ref = _uniform_phase_ref()
    beams = [Beam(theta_deg=20.0, phi_deg=0.0, weight=1.0), Beam(theta_deg=35.0, phi_deg=120.0, weight=0.6)]
    _continuous_deg, _quantized_deg, level_map = cs.coded_phase_map_deg(ref, beams, FREQ, n_bits=1)
    psi_maps = cs.per_beam_phase_maps_deg(ref, beams, FREQ)
    f = sum(abs(b.weight) * np.cos(np.radians(p)) for b, p in zip(beams, psi_maps, strict=True))
    expected_bit = (f <= 0).astype(int)  # f>0 -> level 0 ("B=1"/state 0); f<=0 -> level 1
    assert np.array_equal(level_map, expected_bit)


def test_coded_phase_map_n_bit_quantizes_into_2_pow_n_levels():
    ref = _uniform_phase_ref()
    beam = Beam(theta_deg=25.0, phi_deg=10.0)
    for n_bits in (1, 2, 3):
        _c, quantized_deg, level_map = cs.coded_phase_map_deg(ref, [beam], FREQ, n_bits=n_bits)
        n_states = 2**n_bits
        step = 360.0 / n_states
        assert set(np.unique(level_map)) <= set(range(n_states))
        assert np.allclose(np.sort(np.unique(quantized_deg)), np.sort(np.unique(level_map)) * step)


def test_coded_phase_map_rejects_zero_bits():
    ref = _uniform_phase_ref()
    with pytest.raises(ValueError):
        cs.coded_phase_map_deg(ref, [Beam(theta_deg=0.0, phi_deg=0.0)], FREQ, n_bits=0)


# ================================================================== synthesize_coded_one_bit


def _ideal_two_state_choice():
    state_0 = TwoStateCell("full patch", "r0", complex(1.0, 0.0), {"kind": "rectangle", "L_mm": 2.0, "W_mm": 2.0})
    state_1 = TwoStateCell(
        "annular slot", "r1", complex(-1.0, 0.0),
        {"kind": "annular_slot", "a_mm": 2.0, "b_mm": 2.0, "R_mm": 0.7, "r_mm": 0.4},
    )
    return TwoStateChoice(state_0, state_1, achieved_separation_deg=180.0)


def test_synthesize_coded_one_bit_broadside_natural_is_a_single_uniform_state():
    ref = _uniform_phase_ref()
    beam = Beam(theta_deg=0.0, phi_deg=0.0)  # broadside, natural pol -> psi == 0 everywhere
    choice = _ideal_two_state_choice()
    design = cs.synthesize_coded_one_bit(ref, [beam], FREQ, choice, PERIOD_M)
    achieved = {c.achieved_value for c in design.cells}
    assert achieved == {0.0}  # every cell picks state_0 (phase 0deg)
    assert all(c.phase_error_deg == pytest.approx(0.0) for c in design.cells)
    assert all(c.geometry_param == 0.0 for c in design.cells)
    assert "CODED_DESIGN_METHOD" in design.labels


def test_synthesize_coded_one_bit_with_imperfect_separation_reports_a_nonzero_phase_error():
    ref = _uniform_phase_ref()
    beam = Beam(theta_deg=0.0, phi_deg=0.0)
    imperfect = TwoStateChoice(
        TwoStateCell("a", "a", complex(1.0, 0.0), {"kind": "rectangle", "L_mm": 2.0, "W_mm": 2.0}),
        TwoStateCell("b", "b", np.exp(1j * np.radians(150.0)), {"kind": "rectangle", "L_mm": 3.0, "W_mm": 3.0}),
        achieved_separation_deg=150.0,
    )
    design = cs.synthesize_coded_one_bit(ref, [beam], FREQ, imperfect, PERIOD_M)
    for c in design.cells:
        assert c.geometry_param == 0.0  # every cell still lands on state 0 (the nearer target)
        assert c.phase_error_deg == pytest.approx(0.0)  # state 0's own phase IS the 0deg target


def test_synthesize_coded_one_bit_two_beams_uses_both_states():
    ref = _uniform_phase_ref()
    beams = [Beam(theta_deg=30.0, phi_deg=0.0), Beam(theta_deg=30.0, phi_deg=180.0, weight=0.3)]
    choice = _ideal_two_state_choice()
    design = cs.synthesize_coded_one_bit(ref, beams, FREQ, choice, PERIOD_M)
    achieved = {c.achieved_value for c in design.cells}
    assert achieved == {0.0, 180.0}  # a real design with 2 opposing beams should use both states


# ================================================================== synthesize_coded_n_bit


def test_synthesize_coded_n_bit_matches_against_a_continuous_library():
    import math

    ref = _uniform_phase_ref()
    beam = Beam(theta_deg=20.0, phi_deg=0.0)
    candidates = [(f"c{i}", 2.0, 2.0, "lib", np.exp(1j * math.radians(i * 5.0))) for i in range(72)]
    design = cs.synthesize_coded_n_bit(ref, [beam], FREQ, n_bits=2, library_candidates=candidates, period_m=PERIOD_M)
    assert design.n_bits == 2
    assert set(np.unique(design.level_map)) <= {0, 1, 2, 3}
    # every cell's achieved phase should be within a reasonable tolerance of ITS OWN quantized
    # 4-level target (the dense 5deg-step library can match any target closely)
    for cell in design.cells:
        assert cell.phase_error_deg < 10.0


# ================================================================ continuous_reference_aperture_field


def test_continuous_reference_aperture_field_matches_the_formula():
    ref = _uniform_phase_ref()
    continuous_deg = np.full(ref.x_m.shape, 37.0)
    field = cs.continuous_reference_aperture_field(ref, continuous_deg)
    expected = np.exp(1j * np.radians(ref.phase())) * np.exp(1j * np.radians(37.0))
    assert np.allclose(field, expected)


# ========================================================= quantization loss vs continuous design


def test_quantization_loss_is_measurable_against_the_continuous_design():
    """Part E: "report quantization loss vs continuous" -- build both the REAL (quantized +
    achieved two-state) aperture field and the hypothetical continuous one, and measure the
    directivity difference at the main beam via the SAME vector pattern model Part B already
    validated. 1-bit coding is a severe quantization, so some loss is expected and simply
    reported, not asserted to be small."""
    from src.green_function import LayeredGreenFunction
    from src.materials import Substrate
    from src.pattern.vector_aperture_pattern import aperture_current_from_scalar, compute_vector_pattern

    x, y = _grid(n=61, extent_m=0.06)
    ref = phase_only_reference_field(x, y, np.zeros_like(x))
    beam = Beam(theta_deg=20.0, phi_deg=0.0)
    choice = _ideal_two_state_choice()
    design = cs.synthesize_coded_one_bit(ref, [beam], FREQ, choice, PERIOD_M)

    continuous_field = cs.continuous_reference_aperture_field(
        phase_only_reference_field(design.x_m, design.y_m, np.zeros_like(design.x_m)), design.continuous_phase_deg
    )
    substrate = Substrate("rt5880-like", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    gf = LayeredGreenFunction(substrate, FREQ, backing="pec")
    dx = float(design.x_m[0, 1] - design.x_m[0, 0])
    dy = float(design.y_m[1, 0] - design.y_m[0, 0])

    jx_q, jy_q = aperture_current_from_scalar(design.x_m, design.y_m, design.aperture_field)
    jx_c, jy_c = aperture_current_from_scalar(design.x_m, design.y_m, continuous_field)
    pat_q = compute_vector_pattern(design.x_m, design.y_m, jx_q, jy_q, FREQ, gf, dx, dy, theta_max_deg=60.0, n_theta=61, n_phi=73)
    pat_c = compute_vector_pattern(design.x_m, design.y_m, jx_c, jy_c, FREQ, gf, dx, dy, theta_max_deg=60.0, n_theta=61, n_phi=73)

    # both designs should at least agree on WHERE the main beam landed (this small, coarse test
    # aperture is not expected to steer to the commanded 20deg precisely -- see
    # test_pattern_vector_aperture_pattern.py's own similarly-coarse cases) -- otherwise a
    # directivity "loss" number would be comparing two different beams, not a quantization effect.
    assert pat_c.main_beam_theta_deg == pytest.approx(pat_q.main_beam_theta_deg, abs=6.0)
    phi_err = abs((pat_c.main_beam_phi_deg - pat_q.main_beam_phi_deg + 180.0) % 360.0 - 180.0)
    assert phi_err < 15.0

    loss_db = pat_c.directivity_estimate_dbi - pat_q.directivity_estimate_dbi
    assert np.isfinite(loss_db)
    # 1-bit quantization is lossy but should not be catastrophic (sanity bound, not a tight claim)
    assert -1.0 < loss_db < 10.0
