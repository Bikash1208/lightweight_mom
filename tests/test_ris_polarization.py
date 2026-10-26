"""Phase 20 R3 item 2: circular-basis reflection (handedness relative to each wave's own propagation direction)."""

import math

import numpy as np
import pytest

from src.excitation import ReferencePlanes
from src.green_function import LayeredGreenFunction
from src.materials import Substrate
from src.ris.polarization import L, R, cartesian_reflection, circular_reflection, jones_to_te_tm
from src.sparameters import oblique_stack_s_parameters


def _te_tm(r_cart):
    return jones_to_te_tm(r_cart, 0.0)


def test_metal_plate_flips_handedness():
    c = circular_reflection(-np.eye(2))
    assert abs(c[R, R]) < 1e-15 and abs(c[L, L]) < 1e-15
    assert c[L, R] == pytest.approx(-1.0) and c[R, L] == pytest.approx(-1.0)


def test_grounded_slab_at_normal_incidence_flips_handedness_and_at_oblique_is_elliptic_but_lossless():
    sub = Substrate("s", 2.2, 0.0, 1.575e-3)
    g = LayeredGreenFunction(sub, 10e9, backing="pec")
    planes = ReferencePlanes(0.0, None, sub.thickness)
    s0 = oblique_stack_s_parameters(g, planes, 0.0, 0.0).s11
    c0 = circular_reflection(s0)
    assert abs(c0[R, R]) < 1e-12 and abs(c0[L, R]) == pytest.approx(1.0, abs=1e-12)
    s30 = oblique_stack_s_parameters(g, planes, math.radians(30), 0.0).s11
    c30 = circular_reflection(s30)
    assert (
        abs(c30[R, R]) > 1e-3
    )  # TE and TM phases differ at oblique incidence: the reflected wave is no longer purely opposite-handed
    # unitarity (lossless, single reflected beam): the circular transform preserves the 2-norm of the incident-polarisation column
    assert np.linalg.norm(c30[:, R]) == pytest.approx(1.0, abs=1e-9)


@pytest.mark.parametrize("alpha_deg", [0.0, 15.0, 45.0, 90.0, 135.0])
def test_half_wave_plate_rotated_by_alpha_gives_phase_two_alpha(alpha_deg):
    """R = Rot(a) diag(1, -1) Rot(-a): the handedness-preserving terms are r_RR = exp(+j 2a) and r_LL = exp(-j 2a); no flipped terms."""
    a = math.radians(alpha_deg)
    rot = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
    r_cart = rot @ np.diag([1.0, -1.0]) @ rot.T
    c = circular_reflection(_te_tm(r_cart))
    assert c[R, R] == pytest.approx(np.exp(2j * a), abs=1e-12) and c[L, L] == pytest.approx(np.exp(-2j * a), abs=1e-12)
    assert abs(c[L, R]) < 1e-12 and abs(c[R, L]) < 1e-12


def test_cartesian_roundtrip_and_azimuth_rotation():
    rng = np.random.default_rng(3)
    m = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    for phi in (0.0, 0.4, 2.0):
        assert jones_to_te_tm(cartesian_reflection(m, phi), phi) == pytest.approx(m, abs=1e-12)
    # TE/TM at phi = 0: TM = x, TE = y
    r = np.diag([2.0, 3.0]).astype(complex)  # x, y
    s = jones_to_te_tm(r, 0.0)
    assert s[0, 0] == 3.0 and s[1, 1] == 2.0  # (TE, TM) = (y, x)
    with pytest.raises(ValueError):
        circular_reflection(np.zeros((3, 3)))


def test_te_tm_from_circular_is_the_exact_inverse():
    from src.ris.polarization import te_tm_from_circular

    rng = np.random.default_rng(3)
    s = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    assert np.allclose(te_tm_from_circular(circular_reflection(s)), s, atol=1e-14)


def test_structure_frame_turns_with_azimuth_and_pb_phase_is_explicit():
    from src.ris.polarization import cartesian_reflection, from_structure_frame, to_structure_frame

    # a half-wave plate turned by alpha reflects with r_RR = exp(+2j alpha) in the fixed structure frame (R3 convention); the incidence-plane frame at
    # azimuth phi sees it as exp(2j (alpha - phi)): to_structure_frame restores the first
    alpha = math.radians(37.0)
    jones = np.array([[math.cos(2 * alpha), math.sin(2 * alpha)], [math.sin(2 * alpha), -math.cos(2 * alpha)]])  # [[cos2a, sin2a], [sin2a, -cos2a]]
    ref = circular_reflection(jones_to_te_tm(jones, 0.0))
    for phi in (0.0, math.radians(50.0), math.radians(-120.0)):
        s = jones_to_te_tm(jones, phi)  # the same structure seen from azimuth phi (normal incidence)
        c_o = circular_reflection(s)
        assert np.allclose(to_structure_frame(c_o, phi), ref, atol=1e-12)
        assert np.allclose(from_structure_frame(to_structure_frame(c_o, phi), phi), c_o, atol=1e-14)
    assert cartesian_reflection(jones_to_te_tm(jones, 0.7), 0.7) == pytest.approx(jones, abs=1e-12)
