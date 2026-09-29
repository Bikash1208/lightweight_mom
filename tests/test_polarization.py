"""Circular basis, co/cross-polar conversion and handedness under exp(+j w t).

The handedness checks do not call the code under test to decide what "right-handed" means: they
build the physical space-time field E(z, t) = Re{e exp(j (w t -+ k z))}, verify which way it
travels, and evaluate the rotation about the propagation direction with a 3-D cross product.
"""

import numpy as np
import pytest

from src.polarization import (
    CIRCULAR_BASIS,
    CIRCULAR_COMPONENTS,
    CO_COMPONENTS,
    CROSS_COMPONENTS,
    E_L,
    E_R,
    LINEAR_COMPONENTS,
    amplitudes_from_circular,
    amplitudes_to_circular,
    block_from_circular,
    block_to_circular,
    co_cross_power,
    component_handedness,
    handedness,
    matrix4_from_circular,
    matrix4_to_circular,
    rotation_sense_about_z,
    wave_direction,
)

RNG = np.random.default_rng(7)


def random_block(n=None):
    shape = (2, 2) if n is None else (n, 2, 2)
    return RNG.normal(size=shape) + 1j * RNG.normal(size=shape)


def rot(psi):
    c, s = np.cos(psi), np.sin(psi)
    return np.array([[c, -s], [s, c]])


# ------------------------------------------------------------------------ the basis itself


def test_basis_vectors_are_exactly_the_specified_ones():
    assert np.allclose(E_R, np.array([1.0, 1.0j]) / np.sqrt(2.0), atol=0, rtol=1e-15)
    assert np.allclose(E_L, np.array([1.0, -1.0j]) / np.sqrt(2.0), atol=0, rtol=1e-15)
    assert np.array_equal(CIRCULAR_BASIS[:, 0], E_R) and np.array_equal(CIRCULAR_BASIS[:, 1], E_L)


def test_basis_is_orthonormal_and_unitary():
    t = CIRCULAR_BASIS
    assert np.allclose(t.conj().T @ t, np.eye(2), atol=1e-15)
    assert np.allclose(t @ t.conj().T, np.eye(2), atol=1e-15)
    assert abs(np.vdot(E_R, E_L)) < 1e-15  # e_R^H e_L = 0
    assert np.linalg.norm(E_R) == pytest.approx(1.0, abs=1e-15)


def test_basis_arrays_are_read_only():
    with pytest.raises(ValueError):
        E_R[0] = 0.0
    with pytest.raises(ValueError):
        CIRCULAR_BASIS[0, 0] = 0.0


def test_conjugate_of_R_is_L():
    """Time reversal / conjugation swaps the labels: the fact behind the reciprocity relation."""
    assert np.array_equal(np.conj(E_R), E_L)


# ------------------------------------------------------------ closed forms and unitarity


def test_block_conversion_matches_hand_derived_closed_forms():
    b = random_block()
    xx, xy, yx, yy = b[0, 0], b[0, 1], b[1, 0], b[1, 1]
    c = block_to_circular(b)
    assert np.isclose(c[0, 0], (xx + yy + 1j * (xy - yx)) / 2, atol=1e-14)  # RR
    assert np.isclose(c[1, 0], (xx - yy + 1j * (xy + yx)) / 2, atol=1e-14)  # LR (received L, incident R)
    assert np.isclose(c[0, 1], (xx - yy - 1j * (xy + yx)) / 2, atol=1e-14)  # RL
    assert np.isclose(c[1, 1], (xx + yy - 1j * (xy - yx)) / 2, atol=1e-14)  # LL


def test_component_index_tables_are_received_then_incident():
    assert LINEAR_COMPONENTS == {"xx": (0, 0), "yx": (1, 0), "xy": (0, 1), "yy": (1, 1)}
    assert CIRCULAR_COMPONENTS == {"RR": (0, 0), "LR": (1, 0), "RL": (0, 1), "LL": (1, 1)}
    b = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=complex)  # xx=1, xy=2, yx=3, yy=4
    assert b[LINEAR_COMPONENTS["yx"]] == 3.0 and b[LINEAR_COMPONENTS["xy"]] == 2.0


def test_conversion_round_trips_and_stacks():
    b = random_block(5)
    assert np.allclose(block_from_circular(block_to_circular(b)), b, atol=1e-14)
    assert np.allclose(block_to_circular(block_from_circular(b)), b, atol=1e-14)
    assert np.allclose(block_to_circular(b)[3], block_to_circular(b[3]), atol=1e-15)  # leading axes


def test_conversion_is_a_unitary_similarity_so_it_preserves_norms_and_singular_values():
    b = random_block()
    c = block_to_circular(b)
    assert np.linalg.norm(c) == pytest.approx(np.linalg.norm(b), rel=1e-14)
    assert np.allclose(np.linalg.svd(c, compute_uv=False), np.linalg.svd(b, compute_uv=False), atol=1e-13)
    assert np.trace(c) == pytest.approx(np.trace(b), abs=1e-13)  # similarity: same trace
    assert np.linalg.det(c) == pytest.approx(np.linalg.det(b), abs=1e-13)


def test_four_by_four_conversion_is_the_blockwise_conversion():
    m = RNG.normal(size=(4, 4)) + 1j * RNG.normal(size=(4, 4))
    c = matrix4_to_circular(m)
    for i, j in ((0, 0), (0, 1), (1, 0), (1, 1)):
        assert np.allclose(
            c[2 * i : 2 * i + 2, 2 * j : 2 * j + 2],
            block_to_circular(m[2 * i : 2 * i + 2, 2 * j : 2 * j + 2]),
            atol=1e-14,
        )
    assert np.allclose(matrix4_from_circular(c), m, atol=1e-14)


def test_conversion_rejects_wrong_shapes():
    for bad in (np.zeros((3, 3)), np.zeros(4), np.zeros((2, 3))):
        with pytest.raises(ValueError):
            block_to_circular(bad)
    with pytest.raises(ValueError):
        matrix4_to_circular(np.zeros((2, 2)))


def test_amplitude_conversion_and_known_states():
    assert np.allclose(amplitudes_to_circular(E_R), [1.0, 0.0], atol=1e-15)
    assert np.allclose(amplitudes_to_circular(E_L), [0.0, 1.0], atol=1e-15)
    assert np.allclose(amplitudes_to_circular([1.0, 0.0]), [1 / np.sqrt(2)] * 2, atol=1e-15)  # x = (R + L)/sqrt2
    a = amplitudes_to_circular(np.array([1.0, 1.0]) / np.sqrt(2.0))  # 45-degree linear
    assert abs(a[0]) == pytest.approx(abs(a[1])) == pytest.approx(1 / np.sqrt(2))
    v = RNG.normal(size=2) + 1j * RNG.normal(size=2)
    assert np.allclose(amplitudes_from_circular(amplitudes_to_circular(v)), v, atol=1e-14)
    assert np.linalg.norm(amplitudes_to_circular(v)) == pytest.approx(np.linalg.norm(v))
    with pytest.raises(ValueError):
        amplitudes_to_circular([1.0, 2.0, 3.0])


def test_wave_plates_have_the_textbook_effect_in_this_basis():
    half = np.diag([1.0, -1.0]).astype(complex)  # half-wave plate: pure R<->L conversion
    c = block_to_circular(half)
    assert np.allclose(c, [[0, 1], [1, 0]], atol=1e-15)
    quarter = np.diag([1.0, -1.0j])  # y delayed by 90 deg under exp(+jwt)
    out = quarter @ (np.array([1.0, 1.0]) / np.sqrt(2.0))  # 45-degree linear in
    assert np.allclose(out, E_L, atol=1e-15)  # (x - j y)/sqrt2 out
    assert np.allclose(block_to_circular(np.eye(2)), np.eye(2), atol=1e-15)  # identity stays diagonal


# --------------------------------------------------------------------- co / cross polar


def test_co_and_cross_component_tables():
    assert CO_COMPONENTS == {"linear": ("xx", "yy"), "circular": ("RR", "LL")}
    assert CROSS_COMPONENTS == {"linear": ("yx", "xy"), "circular": ("LR", "RL")}
    for basis, table in (("linear", LINEAR_COMPONENTS), ("circular", CIRCULAR_COMPONENTS)):
        for name in CO_COMPONENTS[basis]:
            assert table[name][0] == table[name][1]  # same label received and incident
        for name in CROSS_COMPONENTS[basis]:
            assert table[name][0] != table[name][1]


def test_co_cross_power_on_known_matrices_and_conservation():
    b = np.array([[0.5, 0.2j], [0.3, -0.4]], dtype=complex)  # xx, xy, yx, yy
    co, cross = co_cross_power(b, "linear")
    assert co == pytest.approx(0.5**2 + 0.4**2) and cross == pytest.approx(0.2**2 + 0.3**2)
    for _ in range(20):
        r = random_block()
        cl, xl = co_cross_power(r, "linear")
        cc, xc = co_cross_power(block_to_circular(r), "circular")
        assert cl + xl == pytest.approx(cc + xc, rel=1e-13)  # total power is basis independent


def test_co_polar_in_one_basis_is_cross_polar_in_the_other_for_the_wave_plate():
    """diag(1,-1): purely co-polar linearly, purely cross-polar circularly."""
    half = np.diag([1.0, -1.0]).astype(complex)
    assert co_cross_power(half, "linear")[1] == 0.0
    assert co_cross_power(block_to_circular(half), "circular")[0] == pytest.approx(0.0, abs=1e-15)


def test_co_cross_power_accepts_stacks_and_rejects_bad_basis():
    b = random_block(4)
    co, cross = co_cross_power(b, "linear")
    assert co.shape == cross.shape == (4,)
    with pytest.raises(ValueError):
        co_cross_power(b, "elliptical")


# ------------------------------------------------- geometric (Pancharatnam-Berry) phase law


@pytest.mark.parametrize("psi_deg", [0.0, 15.0, 45.0, 90.0, 137.0, -60.0])
def test_rotating_a_structure_by_psi_gives_the_geometric_phase_law(psi_deg):
    """Structure rotated counter-clockwise (seen from +z) by psi: S' = R(psi) S R(psi)^T.
    Then RR, LL unchanged; LR' = exp(+2j psi) LR; RL' = exp(-2j psi) RL."""
    psi = np.radians(psi_deg)
    s = random_block()
    c0 = block_to_circular(s)
    c1 = block_to_circular(rot(psi) @ s @ rot(psi).T)
    assert np.isclose(c1[0, 0], c0[0, 0], atol=1e-13)
    assert np.isclose(c1[1, 1], c0[1, 1], atol=1e-13)
    assert np.isclose(c1[1, 0], np.exp(2j * psi) * c0[1, 0], atol=1e-13)
    assert np.isclose(c1[0, 1], np.exp(-2j * psi) * c0[0, 1], atol=1e-13)


def test_basis_vectors_are_eigenvectors_of_rotation_with_opposite_phase():
    psi = 0.7
    assert np.allclose(rot(psi) @ E_R, np.exp(-1j * psi) * E_R, atol=1e-15)
    assert np.allclose(rot(psi) @ E_L, np.exp(+1j * psi) * E_L, atol=1e-15)


# ------------------------------------------------------------------ rotation and handedness


def field_at(e, t, z, direction, omega=2.0 * np.pi, c=1.0):
    """Physical field Re{e exp(j (w t -+ k z))}: '+z' uses exp(-jkz) (project convention)."""
    k = omega / c
    phase = omega * t - k * z if direction == "+z" else omega * t + k * z
    return np.real(np.outer(np.exp(1j * phase), e) if np.ndim(t) else np.exp(1j * phase) * e)


def test_wave_directions_follow_from_the_time_convention():
    """Under exp(+j w t), exp(-j k z) is a +z wave: the pattern at (z + d, t + d/c) equals (z, t)."""
    e = E_R
    for direction, sign in (("+z", +1), ("-z", -1)):
        z, t, d = 0.31, 0.17, 0.05
        assert np.allclose(field_at(e, t + d, z + sign * d, direction), field_at(e, t, z, direction), atol=1e-12)
        assert not np.allclose(field_at(e, t + d, z - sign * d, direction), field_at(e, t, z, direction), atol=1e-3)


def test_R_rotates_clockwise_seen_from_plus_z_by_direct_sampling():
    t = np.linspace(0.0, 1.0, 9)[:-1]  # one period, w = 2 pi
    for label, vec, expected in (("R", E_R, -1), ("L", E_L, +1)):
        e_t = field_at(vec, t, 0.0, "+z")
        angle = np.unwrap(np.arctan2(e_t[:, 1], e_t[:, 0]))
        assert np.all(np.sign(np.diff(angle)) == expected), label
        assert e_t[0] == pytest.approx([1 / np.sqrt(2), 0.0], abs=1e-12)  # starts along +x
        assert e_t[2] == pytest.approx([0.0, -1 / np.sqrt(2) if expected < 0 else 1 / np.sqrt(2)], abs=1e-12)
        assert rotation_sense_about_z(vec) == expected


def test_rotation_sense_matches_brute_force_for_random_ellipses_and_flags_linear():
    t = np.linspace(0.0, 1.0, 4001)[:-1]
    for _ in range(30):
        e = RNG.normal(size=2) + 1j * RNG.normal(size=2)
        e_t = field_at(e, t, 0.0, "+z")
        torque = e_t[:-1, 0] * np.diff(e_t[:, 1]) - e_t[:-1, 1] * np.diff(e_t[:, 0])  # z of E x dE
        assert rotation_sense_about_z(e) == (1 if torque.mean() > 0 else -1)
    assert rotation_sense_about_z([1.0, 0.0]) == 0
    assert rotation_sense_about_z([1.0, 1.0]) == 0
    assert rotation_sense_about_z([1.0, np.exp(0.7j)]) != 0  # 0.7 rad phase difference: elliptical
    with pytest.raises(ValueError):
        rotation_sense_about_z([0.0, 0.0])


def brute_force_handedness(label, direction):
    """RHCP iff n . (E x dE/dt) > 0 for the physical field of a wave travelling along n."""
    e = E_R if label == "R" else E_L
    t = np.linspace(0.0, 1.0, 2001)[:-1]
    e_t = np.column_stack([field_at(e, t, 0.0, direction), np.zeros(t.size)])  # (Ex, Ey, 0)
    de = np.gradient(e_t, t, axis=0)
    n = np.array([0.0, 0.0, 1.0 if direction == "+z" else -1.0])
    return "RHCP" if np.mean(np.cross(e_t, de) @ n) > 0 else "LHCP"


@pytest.mark.parametrize("label", ["R", "L"])
@pytest.mark.parametrize("direction", ["+z", "-z"])
def test_handedness_matches_right_hand_rule_evaluated_on_the_physical_field(label, direction):
    assert handedness(label, direction) == brute_force_handedness(label, direction)


def test_documented_handedness_table():
    assert handedness("R", "-z") == "RHCP" and handedness("L", "-z") == "LHCP"
    assert handedness("R", "+z") == "LHCP" and handedness("L", "+z") == "RHCP"


def test_wave_direction_of_each_port_role():
    assert wave_direction("top", "incident") == "-z"
    assert wave_direction("top", "outgoing") == "+z"
    assert wave_direction("bottom", "incident") == "+z"
    assert wave_direction("bottom", "outgoing") == "-z"
    for args in (("side", "incident"), ("top", "reflected")):
        with pytest.raises(ValueError):
            wave_direction(*args)


def test_component_handedness_of_every_circular_parameter():
    """Handedness follows the propagation direction only; the fixed-basis label is unchanged."""

    def flip(h):
        return "LHCP" if h == "RHCP" else "RHCP"

    for fam, (out_port, in_port) in {"S11": (1, 1), "S12": (1, 2), "S21": (2, 1), "S22": (2, 2)}.items():
        in_dir = "-z" if in_port == 1 else "+z"
        out_dir = "+z" if out_port == 1 else "-z"
        for comp in CIRCULAR_COMPONENTS:
            received, incident = comp[0], comp[1]
            h_in = "RHCP" if incident == "R" else "LHCP"
            h_out = "RHCP" if received == "R" else "LHCP"
            if in_dir == "+z":
                h_in = flip(h_in)
            if out_dir == "+z":
                h_out = flip(h_out)
            assert component_handedness(f"{fam}_{comp}") == (h_in, h_out), (fam, comp)


def test_transmission_keeps_handedness_reflection_reverses_it_for_same_label():
    assert component_handedness("S21_RR") == ("RHCP", "RHCP")
    assert component_handedness("S21_LL") == ("LHCP", "LHCP")
    assert component_handedness("S11_RR") == ("RHCP", "LHCP")
    assert component_handedness("S22_LL") == ("RHCP", "LHCP")  # bottom port: incident +z
    assert component_handedness("S21_LR") == ("RHCP", "LHCP")  # label conversion = handedness flip


@pytest.mark.parametrize("bad", ["S21_xx", "S31_RR", "S21", "T21_RR", "S21_RRR", "S2_RR", ""])
def test_component_handedness_rejects_bad_names(bad):
    with pytest.raises(ValueError):
        component_handedness(bad)


def test_handedness_rejects_bad_arguments():
    for args in (("x", "+z"), ("R", "z"), ("R", "up")):
        with pytest.raises(ValueError):
            handedness(*args)
