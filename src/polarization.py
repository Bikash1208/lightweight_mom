"""Linear <-> circular polarization basis, co/cross-polar bookkeeping and handedness.

Time convention exp(+j w t): a phasor E represents the physical field Re{E exp(+j w t)}.

Basis (fixed by the project, not inferred from any software package)
---------------------------------------------------------------------
    e_R = (e_x + j e_y) / sqrt(2)          e_L = (e_x - j e_y) / sqrt(2)

The basis vectors are *global* x-y vectors, identical at both ports and for incident and
scattered waves: the label R/L is a label of the transverse vector, not of the propagation
direction. The matrix T = [e_R  e_L] is unitary, so total power is basis independent.

Rotation in time (the fact everything else follows from)
--------------------------------------------------------
For a transverse phasor e = (e_x, e_y) the physical field at a fixed point is
E(t) = Re{e exp(+j w t)}, and its rotation rate about +z, averaged over a period, is

    <E_x dE_y/dt - E_y dE_x/dt> = w * Im(e_x conj(e_y)).

For e_R: e_x conj(e_y) = (1/2)(-j), so the sign is NEGATIVE. E(t) = (cos wt, -sin wt)/sqrt(2):
at t = 0 it points along +x and a quarter period later along -y. Hence

    e_R rotates CLOCKWISE seen from +z looking back at the origin (from +x toward -y),
    e_L rotates COUNTER-CLOCKWISE seen from +z.

Handedness of a travelling wave (IEEE/right-hand-rule definition)
------------------------------------------------------------------
"Right-handed" (RHCP): with the right thumb along the propagation direction, the fingers
curl in the direction in which E rotates in time. (The optics convention names the same
waves the other way round; this project uses only the right-hand-rule definition.)
The rotation of e_R is clockwise about +z, i.e. counter-clockwise about -z, so

    wave travelling toward -z (the incident wave at the top port, the transmitted wave at
    the bottom port):   e_R is RHCP,  e_L is LHCP
    wave travelling toward +z (the reflected wave at the top port, the incident wave at
    the bottom port):   e_R is LHCP,  e_L is RHCP

A mirror-like reflection therefore keeps the basis label but reverses the physical
handedness: for a full PEC sheet S11 = -I in the x-y basis, hence S11_RR = -1 and
S11_LR = 0 (labels), i.e. an RHCP incident wave returns as LHCP.

Circular S-parameters
---------------------
    S_circ = T^H S_lin T,   entries named [received][incident], e.g. S21_LR = L-labelled
    component at port 2 for an R-labelled wave incident on port 1.

In terms of the linear entries (xx, xy, yx, yy = [received][incident]):

    RR = (xx + yy + j (xy - yx)) / 2        RL = (xx - yy - j (xy + yx)) / 2
    LR = (xx - yy + j (xy + yx)) / 2        LL = (xx + yy - j (xy - yx)) / 2

Rotating the structure counter-clockwise (seen from +z) by psi multiplies RR and LL by 1
and LR by exp(+2j psi), RL by exp(-2j psi) (geometric / Pancharatnam-Berry phase).
"""

import math

import numpy as np

SQRT_HALF = math.sqrt(0.5)

E_R = np.array([1.0, 1.0j]) * SQRT_HALF
E_L = np.array([1.0, -1.0j]) * SQRT_HALF
E_R.setflags(write=False)
E_L.setflags(write=False)

# Columns are the basis vectors: E_lin = CIRCULAR_BASIS @ (a_R, a_L),  a_circ = CIRCULAR_BASIS^H E_lin.
CIRCULAR_BASIS = np.stack([E_R, E_L], axis=1)
CIRCULAR_BASIS.setflags(write=False)

# Two-letter component suffixes are [received][incident]; value = (row, column) in the 2x2 block.
LINEAR_COMPONENTS = {"xx": (0, 0), "yx": (1, 0), "xy": (0, 1), "yy": (1, 1)}
CIRCULAR_COMPONENTS = {"RR": (0, 0), "LR": (1, 0), "RL": (0, 1), "LL": (1, 1)}
CO_COMPONENTS = {"linear": ("xx", "yy"), "circular": ("RR", "LL")}
CROSS_COMPONENTS = {"linear": ("yx", "xy"), "circular": ("LR", "RL")}

_BLOCK4 = np.kron(np.eye(2), CIRCULAR_BASIS)      # blockdiag(T, T) over (port1 x/y, port2 x/y)


def _check_block(name, a, size):
    a = np.asarray(a, dtype=complex)
    if a.ndim < 2 or a.shape[-2:] != (size, size):
        raise ValueError(f"{name} must have trailing shape ({size}, {size}), got {a.shape}")
    return a


def block_to_circular(block):
    """T^H B T for the trailing (2, 2) block(s): linear (x, y) -> circular (R, L) entries."""
    b = _check_block("block", block, 2)
    return CIRCULAR_BASIS.conj().T @ b @ CIRCULAR_BASIS


def block_from_circular(block):
    """Inverse of ``block_to_circular``: T B T^H."""
    b = _check_block("block", block, 2)
    return CIRCULAR_BASIS @ b @ CIRCULAR_BASIS.conj().T


def matrix4_to_circular(matrix):
    """Same change of basis for the 4x4 two-port matrix [[S11, S12], [S21, S22]]."""
    m = _check_block("matrix", matrix, 4)
    return _BLOCK4.conj().T @ m @ _BLOCK4


def matrix4_from_circular(matrix):
    m = _check_block("matrix", matrix, 4)
    return _BLOCK4 @ m @ _BLOCK4.conj().T


def amplitudes_to_circular(e_xy):
    """(a_R, a_L) = T^H (E_x, E_y) for trailing axis of length 2."""
    e = np.asarray(e_xy, dtype=complex)
    if e.shape[-1] != 2:
        raise ValueError(f"expected trailing axis of length 2, got {e.shape}")
    return e @ CIRCULAR_BASIS.conj()


def amplitudes_from_circular(a_rl):
    """(E_x, E_y) = T (a_R, a_L)."""
    a = np.asarray(a_rl, dtype=complex)
    if a.shape[-1] != 2:
        raise ValueError(f"expected trailing axis of length 2, got {a.shape}")
    return a @ CIRCULAR_BASIS.T


def co_cross_power(block, basis):
    """(sum |co|^2, sum |cross|^2) over the co- and cross-polarized entries of a 2x2 block.

    'co' = received label equals incident label (xx, yy | RR, LL); 'cross' = the other two.
    The sum of both is the squared Frobenius norm, which the basis change preserves.
    """
    if basis not in CO_COMPONENTS:
        raise ValueError(f"basis must be one of {sorted(CO_COMPONENTS)}, got {basis!r}")
    b = _check_block("block", block, 2)
    table = LINEAR_COMPONENTS if basis == "linear" else CIRCULAR_COMPONENTS
    co = sum(np.abs(b[..., table[c][0], table[c][1]]) ** 2 for c in CO_COMPONENTS[basis])
    cross = sum(np.abs(b[..., table[c][0], table[c][1]]) ** 2 for c in CROSS_COMPONENTS[basis])
    return co, cross


# ------------------------------------------------------------------ rotation / handedness

def rotation_sense_about_z(e_xy):
    """+1 / -1 / 0: E(t) = Re{e exp(+j w t)} rotates counter-clockwise / clockwise seen from +z
    / does not rotate (linear). Sign of Im(e_x conj(e_y)), evaluated with a tolerance relative
    to |e|^2 so rounding noise in a linear phasor is not read as rotation."""
    ex, ey = np.asarray(e_xy, dtype=complex)
    scale = abs(ex) ** 2 + abs(ey) ** 2
    if scale == 0.0:
        raise ValueError("zero polarization vector has no rotation sense")
    s = (ex * np.conj(ey)).imag / scale
    return 0 if abs(s) < 1e-12 else (1 if s > 0 else -1)


_DIRECTIONS = ("+z", "-z")


def wave_direction(port, role):
    """Propagation direction of the incident/outgoing wave at a port ('+z' or '-z').

    top port: incident travels -z (toward the pattern), outgoing +z;
    bottom port: incident travels +z, outgoing -z.
    """
    if port not in ("top", "bottom") or role not in ("incident", "outgoing"):
        raise ValueError("port must be 'top' or 'bottom', role 'incident' or 'outgoing'")
    going_up = (port == "top") == (role == "outgoing")
    return "+z" if going_up else "-z"


def handedness(label, direction):
    """'RHCP' or 'LHCP' (right-hand rule, thumb along the propagation direction) of the wave
    whose transverse phasor is e_R (label 'R') or e_L (label 'L') travelling along ``direction``."""
    if label not in ("R", "L") or direction not in _DIRECTIONS:
        raise ValueError(f"label must be 'R' or 'L' and direction '+z' or '-z', got {label!r}, {direction!r}")
    sense = rotation_sense_about_z(E_R if label == "R" else E_L)       # about +z
    about_propagation = sense if direction == "+z" else -sense
    return "RHCP" if about_propagation > 0 else "LHCP"


def component_handedness(name):
    """(incident, outgoing) handedness of a circular S-parameter name such as 'S21_LR'.

    S<out><in>: 'S21' = received at port 2, incident at port 1. Both waves keep their
    fixed-basis label; only the propagation direction decides the physical handedness.
    """
    family, _, comp = name.partition("_")
    if len(family) != 3 or family[0] != "S" or family[1:] not in ("11", "12", "21", "22") \
            or comp not in CIRCULAR_COMPONENTS:
        raise ValueError(f"expected a circular S-parameter name like 'S21_LR', got {name!r}")
    ports = {"1": "top", "2": "bottom"}
    out_port, in_port = ports[family[1]], ports[family[2]]
    received, incident = comp[0], comp[1]
    return (
        handedness(incident, wave_direction(in_port, "incident")),
        handedness(received, wave_direction(out_port, "outgoing")),
    )
