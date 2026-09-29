"""Independent reference for the layered-medium Green's function (TEST ONLY).

Deliberately shares NO derivation with src/green_function.py: it does not use the
transmission-line model, TE/TM impedances or the closed forms. It solves the
full 3-D Maxwell boundary-value problem by plane-wave field matching:

* every region holds plane waves exp(-j K.r) with K = kt*p + s*kappa*z_hat,
* E is perpendicular to K (TE: q; TM: s*kappa*p - kt*z_hat, unnormalised),
* H = K x E / (omega*mu0)              (Faraday's law with d/dr -> -jK),
* tangential E continuous at z = 0 and z = -h, tangential H continuous at
  z = -h, and  z_hat x [H(0+) - H(0-)] = J  at the current sheet z = 0.

Only ``normal_wavenumber`` (the documented kz branch) is shared with the code
under test, and the branch is verified physically here by the energy balance.
"""

import numpy as np

from src.constants import C0, EPS0, MU0
from src.periodic_bc import normal_wavenumber

ZHAT = np.array([0.0, 0.0, 1.0])


def _field_function(terms, h, sheet_side="above"):
    """Return fields(z) -> (E_t, H_t) from plane-wave terms.

    terms: list of dicts {region, amp, E, H, kap, s, zref}; the term contributes
    amp * (E, H)[:2] * exp(-1j*s*kap*(z - zref)) inside its region
    (0: z > 0 [z = 0 per sheet_side], 1: -h <= z < 0 [z = 0 per sheet_side], 2: z < -h).
    """

    def fields(z):
        if z > 0 or (z == 0 and sheet_side == "above"):
            region = 0
        elif z < -h:
            region = 2
        else:
            region = 1
        e = np.zeros(2, dtype=complex)
        hh = np.zeros(2, dtype=complex)
        for t in terms:
            if t["region"] != region:
                continue
            phase = t["amp"] * np.exp(-1j * t["s"] * t["kap"] * (z - t["zref"]))
            e += phase * t["E"][:2]
            hh += phase * t["H"][:2]
        return e, hh

    return fields


def solve_current_sheet(kx, ky, J, f_hz, eps_c, h, branch_sign=1.0):
    """Fields of a Floquet current sheet J = (Jx, Jy) exp(-j(kx x + ky y)) at z = 0.

    Stack: air (z>0) | eps_c (-h<z<0) | air (z<-h). ``branch_sign=-1`` flips the
    sign of every kz (a deliberately wrong branch, used to test the tests).

    Returns a dict with E_t at z=0 ('Et'), and power flows per unit area:
    'P_source', 'P_up', 'P_down', 'P_diss'.
    """
    w = 2 * np.pi * f_hz
    k0 = w / C0
    k1 = k0 * np.sqrt(complex(eps_c))
    kt = float(np.hypot(kx, ky))
    p = np.array([1.0, 0.0, 0.0]) if kt == 0.0 else np.array([kx, ky, 0.0]) / kt
    q = np.cross(ZHAT, p)
    kap0 = branch_sign * normal_wavenumber(k0, kx, ky)
    kap1 = branch_sign * normal_wavenumber(k1, kx, ky)

    def wave(pol, s, kap):
        K = kt * p + s * kap * ZHAT
        E = q.astype(complex) if pol == "TE" else (s * kap * p - kt * ZHAT)
        return E, np.cross(K, E) / (w * MU0)

    # unknown order: [u0 TE, u0 TM, u1 TE, u1 TM, d1 TE, d1 TM, u2 TE, u2 TM]
    cols = [
        (0, wave("TE", +1, kap0)),
        (0, wave("TM", +1, kap0)),  # region 0, up, ref z=0
        (1, wave("TE", +1, kap1)),
        (1, wave("TM", +1, kap1)),  # region 1, up, ref z=0
        (1, wave("TE", -1, kap1)),
        (1, wave("TM", -1, kap1)),  # region 1, down, ref z=0
        (2, wave("TE", -1, kap0)),
        (2, wave("TM", -1, kap0)),  # region 2, down, ref z=-h
    ]
    up1 = np.exp(+1j * kap1 * h)  # up wave  exp(-j kap z) evaluated at z = -h
    dn1 = np.exp(-1j * kap1 * h)  # down wave exp(+j kap z) evaluated at z = -h

    M = np.zeros((8, 8), dtype=complex)
    rhs = np.zeros(8, dtype=complex)
    for c, (region, (E, H)) in enumerate(cols):
        Et, Ht = E[:2], H[:2]
        if c in (0, 1):  # region 0 at z=0: +
            M[0:2, c] += Et
            M[2:4, c] += Ht
        elif c in (2, 3, 4, 5):  # region 1 at z=0: - ;  at z=-h: phase factor
            M[0:2, c] -= Et
            M[2:4, c] -= Ht
            ph = up1 if c in (2, 3) else dn1
            M[4:6, c] += Et * ph
            M[6:8, c] += Ht * ph
        else:  # region 2 at z=-h: -
            M[4:6, c] -= Et
            M[6:8, c] -= Ht
    # z_hat x (H+ - H-) = J  ->  dHx = Jy, dHy = -Jx
    rhs[2], rhs[3] = J[1], -J[0]
    amp = np.linalg.solve(M, rhs)

    E0 = amp[0] * cols[0][1][0] + amp[1] * cols[1][1][0]
    H0 = amp[0] * cols[0][1][1] + amp[1] * cols[1][1][1]
    E2 = amp[6] * cols[6][1][0] + amp[7] * cols[7][1][0]
    H2 = amp[6] * cols[6][1][1] + amp[7] * cols[7][1][1]
    Et0 = E0[:2]

    def sz(E, H):
        return 0.5 * np.real(np.cross(E, np.conj(H)))[2]

    P_source = -0.5 * np.real(np.conj(np.asarray(J, dtype=complex)) @ Et0)
    P_up = sz(E0, H0)
    P_down = -sz(E2, H2)

    # dissipation in the substrate: (omega*eps0/2)*(-Im eps_c)*integral |E|^2 dz
    z = np.linspace(-h, 0.0, 4001)
    Ez = np.zeros((z.size, 3), dtype=complex)
    for c, s in zip((2, 3, 4, 5), (+1, +1, -1, -1)):
        Ez += amp[c] * np.exp(-1j * s * kap1 * z)[:, None] * cols[c][1][0][None, :]
    e2 = np.sum(np.abs(Ez) ** 2, axis=1)
    P_diss = 0.5 * w * EPS0 * (-complex(eps_c).imag) * np.trapezoid(e2, z)

    zref = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -h, -h]
    regions = [0, 0, 1, 1, 1, 1, 2, 2]
    ss = [+1, +1, +1, +1, -1, -1, -1, -1]
    kaps = [kap0, kap0, kap1, kap1, kap1, kap1, kap0, kap0]
    terms = [
        dict(region=regions[c], amp=amp[c], E=cols[c][1][0], H=cols[c][1][1], kap=kaps[c], s=ss[c], zref=zref[c])
        for c in range(8)
    ]

    def fields(z, sheet_side="above"):
        return _field_function(terms, h, sheet_side)(z)

    return dict(Et=Et0, P_source=P_source, P_up=P_up, P_down=P_down, P_diss=P_diss, fields=fields)


def solve_incident_fields(kx, ky, E_inc, f_hz, eps_c, h, port="top"):
    """Unpatterned stack under a unit plane wave; returns ``fields(z) -> (E_t, H_t)``.

    port 'top': wave in air (z>0) travelling toward -z with 3-D E-field vector ``E_inc`` at z = 0.
    port 'bottom': wave in air (z<-h) travelling toward +z with ``E_inc`` at z = -h.
    ``E_inc`` must be perpendicular to the incident wavevector. Independent field matching
    (no sheet, J = 0), the same machinery as ``solve_current_sheet``.
    """
    w = 2 * np.pi * f_hz
    k0 = w / C0
    k1 = k0 * np.sqrt(complex(eps_c))
    kt = float(np.hypot(kx, ky))
    p = np.array([1.0, 0.0, 0.0]) if kt == 0.0 else np.array([kx, ky, 0.0]) / kt
    q = np.cross(ZHAT, p)
    kap0 = normal_wavenumber(k0, kx, ky)
    kap1 = normal_wavenumber(k1, kx, ky)

    def wave(pol, s, kap):
        K = kt * p + s * kap * ZHAT
        E = q.astype(complex) if pol == "TE" else (s * kap * p - kt * ZHAT)
        return E, np.cross(K, E) / (w * MU0)

    s_inc = -1 if port == "top" else +1
    E_te, H_te = wave("TE", s_inc, kap0)
    E_tm, H_tm = wave("TM", s_inc, kap0)
    c_te = np.dot(E_inc, q)
    c_tm = np.dot(E_inc, p) / (s_inc * kap0)  # p-component of the TM vector is s*kap
    E_i = c_te * E_te + c_tm * E_tm
    H_i = c_te * H_te + c_tm * H_tm
    assert np.allclose(E_i, E_inc, atol=1e-12), "E_inc is not transverse to k_inc"

    cols = [
        wave("TE", +1, kap0),
        wave("TM", +1, kap0),  # region 0, up-going (reflected or transmitted)
        wave("TE", +1, kap1),
        wave("TM", +1, kap1),  # region 1 up
        wave("TE", -1, kap1),
        wave("TM", -1, kap1),  # region 1 down
        wave("TE", -1, kap0),
        wave("TM", -1, kap0),  # region 2, down-going (ref z=-h)
    ]
    up1, dn1 = np.exp(+1j * kap1 * h), np.exp(-1j * kap1 * h)
    M = np.zeros((8, 8), dtype=complex)
    rhs = np.zeros(8, dtype=complex)
    for c, (E, H) in enumerate(cols):
        Et, Ht = E[:2], H[:2]
        if c in (0, 1):
            M[0:2, c] += Et
            M[2:4, c] += Ht
        elif c in (2, 3, 4, 5):
            M[0:2, c] -= Et
            M[2:4, c] -= Ht
            ph = up1 if c in (2, 3) else dn1
            M[4:6, c] += Et * ph
            M[6:8, c] += Ht * ph
        else:
            M[4:6, c] -= Et
            M[6:8, c] -= Ht
    if port == "top":
        rhs[0:2], rhs[2:4] = -E_i[:2], -H_i[:2]  # incident term sits in region 0 at z = 0
    else:
        rhs[4:6], rhs[6:8] = +E_i[:2], +H_i[:2]  # incident term sits in region 2 at z = -h
    amp = np.linalg.solve(M, rhs)

    zref = [0.0] * 6 + [-h, -h]
    regions = [0, 0, 1, 1, 1, 1, 2, 2]
    ss = [+1, +1, +1, +1, -1, -1, -1, -1]
    kaps = [kap0, kap0, kap1, kap1, kap1, kap1, kap0, kap0]
    terms = [
        dict(region=regions[c], amp=amp[c], E=cols[c][0], H=cols[c][1], kap=kaps[c], s=ss[c], zref=zref[c])
        for c in range(8)
    ]
    inc_region, inc_zref = (0, 0.0) if port == "top" else (2, -h)
    terms.append(dict(region=inc_region, amp=1.0, E=E_i, H=H_i, kap=kap0, s=s_inc, zref=inc_zref))
    return _field_function(terms, h)


def solve_incident_plane_wave(kx, ky, E_inc, f_hz, eps_c, h):
    """Tangential E (Ex, Ey) at z = 0 of the unpatterned stack for top-port incidence."""
    return solve_incident_fields(kx, ky, E_inc, f_hz, eps_c, h, "top")(0.0)[0]


# ============================================================================================
# Phase 12A: PEC-terminated (grounded) stack -- air (z>0) | eps_c (-h<z<0) | PEC ground (z=-h).
#
# Deliberately shares nothing with src/green_function.py's grounded closed forms: instead of the
# transmission-line/sigma-substitution algebra, this solves the SAME 3-D Maxwell boundary-value
# problem as ``solve_current_sheet``/``solve_incident_fields`` above, just with region 2 (air
# below the ground) removed and its 4 matching equations (E_t, H_t continuous at z=-h) replaced
# by the 2 equations of the PEC boundary condition, E_t(-h) = 0. Everything else (plane-wave
# terms, TE/TM decomposition, H = K x E/(w mu0)) is identical. ``_field_function`` needs no
# change at all for this: with no region-2 terms in ``terms``, its z < -h branch already
# evaluates the (empty) sum over "terms whose region is 2" and returns exactly zero -- precisely
# "no field below the ground" (GREEN_FUNCTION.md Sec. 12A), not a special case added here.
# ============================================================================================


def solve_current_sheet_grounded(kx, ky, J, f_hz, eps_c, h, branch_sign=1.0):
    """Fields of a Floquet current sheet J = (Jx, Jy) exp(-j(kx x + ky y)) at z = 0, PEC-grounded
    stack. Returns 'Et' (E_t at z=0), 'P_source'/'P_up' power flows (no 'P_down': nothing exists
    below the ground), 'P_diss', and 'fields(z, sheet_side)'."""
    w = 2 * np.pi * f_hz
    k0 = w / C0
    k1 = k0 * np.sqrt(complex(eps_c))
    kt = float(np.hypot(kx, ky))
    p = np.array([1.0, 0.0, 0.0]) if kt == 0.0 else np.array([kx, ky, 0.0]) / kt
    q = np.cross(ZHAT, p)
    kap0 = branch_sign * normal_wavenumber(k0, kx, ky)
    kap1 = branch_sign * normal_wavenumber(k1, kx, ky)

    def wave(pol, s, kap):
        K = kt * p + s * kap * ZHAT
        E = q.astype(complex) if pol == "TE" else (s * kap * p - kt * ZHAT)
        return E, np.cross(K, E) / (w * MU0)

    # unknown order: [u0 TE, u0 TM, u1 TE, u1 TM, d1 TE, d1 TM]  (no region-2 unknowns)
    cols = [
        (0, wave("TE", +1, kap0)),
        (0, wave("TM", +1, kap0)),  # region 0, up, ref z=0
        (1, wave("TE", +1, kap1)),
        (1, wave("TM", +1, kap1)),  # region 1, up, ref z=0
        (1, wave("TE", -1, kap1)),
        (1, wave("TM", -1, kap1)),  # region 1, down, ref z=0
    ]
    up1 = np.exp(+1j * kap1 * h)  # up wave  exp(-j kap z) evaluated at z = -h
    dn1 = np.exp(-1j * kap1 * h)  # down wave exp(+j kap z) evaluated at z = -h

    M = np.zeros((6, 6), dtype=complex)
    rhs = np.zeros(6, dtype=complex)
    for c, (region, (E, H)) in enumerate(cols):
        Et, Ht = E[:2], H[:2]
        if c in (0, 1):  # region 0 at z=0: +
            M[0:2, c] += Et
            M[2:4, c] += Ht
        else:  # region 1: - at z=0; E_t(-h) = 0 (PEC ground)
            M[0:2, c] -= Et
            M[2:4, c] -= Ht
            ph = up1 if c in (2, 3) else dn1
            M[4:6, c] += Et * ph
    # z_hat x (H+ - H-) = J at z=0  ->  dHx = Jy, dHy = -Jx; rhs[4:6] = 0 (E_t(-h) = 0)
    rhs[2], rhs[3] = J[1], -J[0]
    amp = np.linalg.solve(M, rhs)

    E0 = amp[0] * cols[0][1][0] + amp[1] * cols[1][1][0]
    H0 = amp[0] * cols[0][1][1] + amp[1] * cols[1][1][1]
    Et0 = E0[:2]

    def sz(E, H):
        return 0.5 * np.real(np.cross(E, np.conj(H)))[2]

    P_source = -0.5 * np.real(np.conj(np.asarray(J, dtype=complex)) @ Et0)
    P_up = sz(E0, H0)

    z = np.linspace(-h, 0.0, 4001)
    Ez = np.zeros((z.size, 3), dtype=complex)
    for c, s in zip((2, 3, 4, 5), (+1, +1, -1, -1)):
        Ez += amp[c] * np.exp(-1j * s * kap1 * z)[:, None] * cols[c][1][0][None, :]
    e2 = np.sum(np.abs(Ez) ** 2, axis=1)
    P_diss = 0.5 * w * EPS0 * (-complex(eps_c).imag) * np.trapezoid(e2, z)

    zref = [0.0] * 6
    regions = [0, 0, 1, 1, 1, 1]
    ss = [+1, +1, +1, +1, -1, -1]
    kaps = [kap0, kap0, kap1, kap1, kap1, kap1]
    terms = [
        dict(region=regions[c], amp=amp[c], E=cols[c][1][0], H=cols[c][1][1], kap=kaps[c], s=ss[c], zref=zref[c])
        for c in range(6)
    ]

    def fields(z, sheet_side="above"):
        return _field_function(terms, h, sheet_side)(z)

    return dict(Et=Et0, P_source=P_source, P_up=P_up, P_diss=P_diss, fields=fields)


def solve_incident_fields_grounded(kx, ky, E_inc, f_hz, eps_c, h):
    """Unpatterned PEC-grounded stack under a unit plane wave incident from the top (there is no
    bottom port); returns ``fields(z) -> (E_t, H_t)``. Independent field matching (no sheet,
    J = 0; E_t(-h) = 0 replaces the region-2 matching of ``solve_incident_fields``)."""
    w = 2 * np.pi * f_hz
    k0 = w / C0
    k1 = k0 * np.sqrt(complex(eps_c))
    kt = float(np.hypot(kx, ky))
    p = np.array([1.0, 0.0, 0.0]) if kt == 0.0 else np.array([kx, ky, 0.0]) / kt
    q = np.cross(ZHAT, p)
    kap0 = normal_wavenumber(k0, kx, ky)
    kap1 = normal_wavenumber(k1, kx, ky)

    def wave(pol, s, kap):
        K = kt * p + s * kap * ZHAT
        E = q.astype(complex) if pol == "TE" else (s * kap * p - kt * ZHAT)
        return E, np.cross(K, E) / (w * MU0)

    s_inc = -1  # top port only: incident wave travels -z
    E_te, H_te = wave("TE", s_inc, kap0)
    E_tm, H_tm = wave("TM", s_inc, kap0)
    c_te = np.dot(E_inc, q)
    c_tm = np.dot(E_inc, p) / (s_inc * kap0)
    E_i = c_te * E_te + c_tm * E_tm
    H_i = c_te * H_te + c_tm * H_tm
    assert np.allclose(E_i, E_inc, atol=1e-12), "E_inc is not transverse to k_inc"

    cols = [
        wave("TE", +1, kap0),
        wave("TM", +1, kap0),  # region 0, up-going (reflected)
        wave("TE", +1, kap1),
        wave("TM", +1, kap1),  # region 1 up
        wave("TE", -1, kap1),
        wave("TM", -1, kap1),  # region 1 down
    ]
    up1, dn1 = np.exp(+1j * kap1 * h), np.exp(-1j * kap1 * h)
    M = np.zeros((6, 6), dtype=complex)
    rhs = np.zeros(6, dtype=complex)
    for c, (E, H) in enumerate(cols):
        Et, Ht = E[:2], H[:2]
        if c in (0, 1):
            M[0:2, c] += Et
            M[2:4, c] += Ht
        else:
            M[0:2, c] -= Et
            M[2:4, c] -= Ht
            ph = up1 if c in (2, 3) else dn1
            M[4:6, c] += Et * ph
    rhs[0:2], rhs[2:4] = -E_i[:2], -H_i[:2]  # incident term sits in region 0 at z = 0
    amp = np.linalg.solve(M, rhs)

    zref = [0.0] * 6
    regions = [0, 0, 1, 1, 1, 1]
    ss = [+1, +1, +1, +1, -1, -1]
    kaps = [kap0, kap0, kap1, kap1, kap1, kap1]
    terms = [
        dict(region=regions[c], amp=amp[c], E=cols[c][0], H=cols[c][1], kap=kaps[c], s=ss[c], zref=zref[c])
        for c in range(6)
    ]
    terms.append(dict(region=0, amp=1.0, E=E_i, H=H_i, kap=kap0, s=s_inc, zref=0.0))
    return _field_function(terms, h)
