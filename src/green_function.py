"""Spectral dyadic Green's function of an air | substrate | air stack.

Derivation, equations and references: see GREEN_FUNCTION.md. Summary:

* Time convention exp(+j*omega*t); fields exp(-j(kx*x + ky*y)).
* An electric surface current J (tangential, A/m) sits in the plane z = 0, which
  is the air/substrate interface (substrate occupies -h < z < 0).
* For each spectral pair (kx, ky) the field splits into TE^z and TM^z parts with
  in-plane unit vectors

      p = (kx, ky)/kt       (TM direction, along k_t)
      q = z x p = (-ky, kx)/kt   (TE direction)

  and each part obeys a transmission-line (TL) equation along z.
* The tangential electric field on the current plane is

      E_t(kx, ky) = G(kx, ky) . J(kx, ky),
      G = -( Z_TM p p^T + Z_TE q q^T )                      [ohm]

  where Z_TE, Z_TM are the parallel combinations of the TL input impedances seen
  looking up (semi-infinite air) and down (substrate of thickness h, then
  semi-infinite air). The minus sign is physical: a sheet in free space radiates
  E = -(eta0/2) J.
* G is a 2x2 complex-symmetric matrix. At kt = 0 the TE/TM directions are
  undefined but Z_TE = Z_TM there, so G = -Z * identity (see ``spectral_dyad``).

G is per unit *Floquet coefficient* of the current: the 1/(Px*Py) factor of the
Floquet series belongs to the later MoM assembly, not to G.

Phase 4 adds tangential E/H at arbitrary z: the field radiated by the sheet
(``sheet_fields``) and the field of the unpatterned stack for a plane wave incident
from either side (``stack_fields``); see S_PARAMETERS.md. Exponentials are always
written with non-positive real exponent (|exp(-2jbh)| <= 1 etc.), so nothing overflows.

NOT here: MoM matrix, S-parameter extraction (mom_matrix.py, sparameters.py).

Phase 12A adds an optional PEC ground plane at z = -h (``backing="pec"``; default "air",
unchanged): air (z>0) | pattern at z=0 | substrate (-h<z<0) | PEC ground at z=-h. Derivation,
closed forms and their numerical verification against the textbook ``tan(bh)`` shorted-line
formula: PHASE12A_GROUND_PLANE.md and GREEN_FUNCTION.md Sec. 4/12A. Summary: grounding replaces
the bottom-face load reflection coefficient (the ``rho`` of ``_mode_constants``, physically "seen
from inside the substrate, looking down at whatever terminates it") with the short-circuit value
-1, everywhere that reflection coefficient is used (the equivalent impedance seen by the sheet,
the internal substrate standing wave, and the field below z=-h, which becomes exactly zero
because the (1+rho) factor that already multiplies the "z < -h" field vanishes at rho=-1 -- no
separate "no field below the ground" code path is needed). The top-of-stack reflection Gamma
(needed for the unpatterned-stack excitation) does NOT reduce to that same substitution, because
its existing closed form assumes a symmetric slab (air both above and below, so the bottom load
reflection is exactly minus the top interface's own reflection); Phase 12A generalizes it to an
arbitrary bottom load reflection (Sec. 12A below) and recovers the existing formula exactly when
the load is air.
"""

from dataclasses import dataclass

import numpy as np

from .constants import C0, EPS0, MU0
from .frequency import angular_frequency, validate_frequency
from .materials import Substrate
from .periodic_bc import RectangularLattice, normal_wavenumber

# |k_t| <= _KT_ZERO_RTOL * k0 is treated as k_t = 0 (the p/q directions are then
# numerically meaningless). There Z_TM - Z_TE = O((kt/k0)^2) ~ 1e-24 relative, so
# the isotropic form is exact to double precision.
_KT_ZERO_RTOL = 1e-12

# kz values that are exactly 0 are replaced by _KZ_FLOOR*k0 in the field formulas (see _fields_kz).
_KZ_FLOOR = 1e-12


def te_tm_unit_vectors(kx, ky, kt_floor=0.0):
    """In-plane unit vectors (p, q) for the TM and TE parts of each (kx, ky).

    p = (kx, ky)/kt, q = z_hat x p = (-ky, kx)/kt. Where kt <= ``kt_floor`` the
    direction is undefined and the arbitrary orthonormal pair p = x_hat,
    q = y_hat is returned (physically irrelevant, since Z_TE = Z_TM at kt = 0).

    Returns arrays ``p`` and ``q`` with a trailing axis of length 2 (x, y).
    """
    kx = np.asarray(kx, dtype=float)
    ky = np.asarray(ky, dtype=float)
    kx, ky = np.broadcast_arrays(kx, ky)
    kt = np.hypot(kx, ky)
    zero = kt <= kt_floor
    kt_safe = np.where(zero, 1.0, kt)
    px = np.where(zero, 1.0, kx / kt_safe)
    py = np.where(zero, 0.0, ky / kt_safe)
    p = np.stack([px, py], axis=-1)
    q = np.stack([-py, px], axis=-1)
    return p, q


_BACKINGS = ("air", "pec")


@dataclass(frozen=True)
class LayeredGreenFunction:
    """Spectral Green's function for air (z>0) | substrate (-h<z<0) | backing (z<-h).

    The current plane is z = 0. All inputs SI; ``frequency`` in Hz. ``backing`` (Phase 12A):
    ``"air"`` (default; every pre-Phase-12A call site, unchanged) terminates the substrate in
    semi-infinite air, as before; ``"pec"`` terminates it in a short circuit (a ground plane) at
    z = -h. Every existing constructor call (positional or keyword, 2 arguments) is untouched.
    """

    substrate: Substrate
    frequency: float
    backing: str = "air"

    def __post_init__(self):
        object.__setattr__(self, "frequency", validate_frequency(self.frequency))
        if self.backing not in _BACKINGS:
            raise ValueError(f"backing must be one of {_BACKINGS}, got {self.backing!r}")

    # ---- derived quantities -------------------------------------------------
    @property
    def omega(self):
        return angular_frequency(self.frequency)

    @property
    def k0(self):
        return self.omega / C0

    @property
    def k_substrate(self):
        """Complex substrate wavenumber, Im <= 0."""
        return self.substrate.wavenumber(self.frequency)

    # ---- kz in each region --------------------------------------------------
    def normal_wavenumbers(self, kx, ky):
        """(kz_air, kz_substrate), both on the branch Im(kz) <= 0."""
        return (
            normal_wavenumber(self.k0, kx, ky),
            normal_wavenumber(self.k_substrate, kx, ky),
        )

    # ---- TE / TM equivalent impedances at the current plane -----------------
    def equivalent_impedances(self, kx, ky):
        """(Z_TE, Z_TM) [ohm]: parallel of the up- and down-looking TL impedances.

        With a = kz_air, b = kz_substrate, eps = eps_c, h = thickness and

            sigma = (1 - exp(-2j*b*h)) / (2*b)      (-> j*h as b -> 0),

        the closed forms (derived in GREEN_FUNCTION.md, Sec. 4) are, for ``backing="air"``
        (default; unchanged since Phase 2-4):

            Z_TE = omega*mu0 * (1 + sigma*(a - b)) / (2a + sigma*(a - b)^2)
            Z_TM = a/(omega*eps0) * (a - b*sigma*(a - b/eps))
                                  / (2a + eps*sigma*(a - b/eps)^2)

        These are algebraically identical to the textbook tan(b*h) input-impedance
        expressions but stay finite at a = 0 or b = 0 and cannot overflow, because
        |exp(-2j*b*h)| <= 1 on the Im(b) <= 0 branch.

        For ``backing="pec"`` (Phase 12A, GREEN_FUNCTION.md Sec. 12A) the bottom-face load is a
        short circuit (reflection coefficient -1) instead of semi-infinite air, giving

            Z_TE = omega*mu0 * sigma / (1 + sigma*(a - b))
            Z_TM = A*B*b*sigma / (omega*eps0 * (A - b*sigma*(A - B))),   A = a, B = b/eps

        derived the same way (substitute the short's Gamma_L = -1 into
        Z_down = Z_sub*(1 + Gamma_L*u)/(1 - Gamma_L*u), form the parallel combination with
        Z_air, then apply the same sigma substitution); verified numerically against the direct
        Z_eq = 1/(1/Z_air + 1/(j*Z_sub*tan(b*h))) form to ~1e-14 (PHASE12A_GROUND_PLANE.md).
        """
        a, b = self.normal_wavenumbers(kx, ky)
        return self._impedances_from_kz(a, b)

    def _impedances_from_kz(self, a, b):
        """The closed forms of ``equivalent_impedances`` for given kz values (a, b)."""
        eps = self.substrate.eps_complex
        h = self.substrate.thickness

        with np.errstate(divide="ignore", invalid="ignore"):
            sigma = np.where(b == 0.0, 1j * h, -np.expm1(-2j * b * h) / (2.0 * np.where(b == 0.0, 1.0, b)))
        if self.backing == "air":
            d = a - b
            z_te = self.omega * MU0 * (1.0 + sigma * d) / (2.0 * a + sigma * d**2)
            d_tm = a - b / eps
            z_tm = (a / (self.omega * EPS0)) * (a - b * sigma * d_tm) / (2.0 * a + eps * sigma * d_tm**2)
        else:  # "pec": bottom-face load reflection Gamma_L = -1 (short circuit at z = -h)
            z_te = self.omega * MU0 * sigma / (1.0 + sigma * (a - b))
            aa, bb = a, b / eps
            z_tm = (aa * bb * b * sigma) / (self.omega * EPS0 * (aa - b * sigma * (aa - bb)))
        return z_te, z_tm

    # ---- the dyad -------------------------------------------------------------
    def spectral_dyad(self, kx, ky):
        """G(kx, ky): array of shape ``broadcast(kx, ky).shape + (2, 2)`` [ohm].

        Ẽ_t = G · J̃ with G = -(Z_TM p pᵀ + Z_TE q qᵀ). Written out with
        kt² = kx² + ky²:

            G_xx = -(Z_TM kx² + Z_TE ky²)/kt²
            G_yy = -(Z_TM ky² + Z_TE kx²)/kt²
            G_xy = G_yx = -(Z_TM - Z_TE) kx ky/kt²

        At kt = 0 (|kt| <= 1e-12 k0), G = -(Z_TE + Z_TM)/2 · identity, which is
        the exact kt -> 0 limit because Z_TE = Z_TM there.
        """
        kx = np.asarray(kx, dtype=float)
        ky = np.asarray(ky, dtype=float)
        kx, ky = np.broadcast_arrays(kx, ky)
        z_te, z_tm = self.equivalent_impedances(kx, ky)

        kt2 = kx**2 + ky**2
        zero = kt2 <= (_KT_ZERO_RTOL * self.k0) ** 2
        kt2_safe = np.where(zero, 1.0, kt2)
        z_avg = 0.5 * (z_te + z_tm)

        gxx = np.where(zero, -z_avg, -(z_tm * kx**2 + z_te * ky**2) / kt2_safe)
        gyy = np.where(zero, -z_avg, -(z_tm * ky**2 + z_te * kx**2) / kt2_safe)
        gxy = np.where(zero, 0.0, -(z_tm - z_te) * kx * ky / kt2_safe)

        g = np.empty(kx.shape + (2, 2), dtype=complex)
        g[..., 0, 0] = gxx
        g[..., 1, 1] = gyy
        g[..., 0, 1] = gxy
        g[..., 1, 0] = gxy
        return g

    # ---- Floquet harmonics ----------------------------------------------------
    def floquet_dyads(self, lattice: RectangularLattice, kx0, ky0, max_order):
        """Green's dyad for every Floquet harmonic |m|, |n| <= max_order.

        Returns ``(kx, ky, G)`` with ``kx, ky`` of shape (2N+1, 2N+1) (index
        [m+N, n+N], ``indexing='ij'``) and ``G`` of shape (2N+1, 2N+1, 2, 2).
        """
        m, n = lattice.mode_grid(max_order)
        kx, ky = lattice.transverse_wavenumbers(kx0, ky0, m, n)
        return kx, ky, self.spectral_dyad(kx, ky)

    # =====================================================================
    # Tangential E and H at arbitrary z  (Phase 4; derivation: S_PARAMETERS.md, Sec. 3)
    # =====================================================================
    #
    # Per (kx, ky) and TE/TM mode, E_t = V(z) e, H_t = I(z) z x e (GREEN_FUNCTION.md Sec. 2);
    # with a, b = kz(air), kz(substrate), Y = 1/Z the characteristic admittances
    # (TE: a/(w mu0), b/(w mu0);  TM: w eps0/a, w eps0 eps_c/b), and
    #
    #     rho = (Z_air - Z_sub)/(Z_air + Z_sub)   [TE: (b-a)/(b+a);  TM: (a-B)/(a+B), B = b/eps_c]
    #     u   = exp(-2j b h)
    #
    # the substrate solution is  V = W [ exp(+jbz) + rho exp(-jb(2h+z)) ]  (both terms bounded),
    # I = Y_sub W [ rho exp(-jb(2h+z)) - exp(+jbz) ].
    # Regions: z > 0 (air above), -h <= z < 0 (substrate), z < -h (air below).

    def _fields_kz(self, kx, ky):
        """(a, b) with an exact zero (Rayleigh point kz = 0) replaced by 1e-12*k0.

        The field formulas divide by kz (I = V/Z_c); the *value* of every field is continuous
        through kz = 0, so the replacement changes nothing beyond ~1e-12 relative to the
        neighbourhood, while avoiding 0/0.
        """
        a, b = self.normal_wavenumbers(kx, ky)
        floor = _KZ_FLOOR * self.k0
        return np.where(a == 0.0, floor, a), np.where(b == 0.0, floor, b)

    def _mode_constants(self, a, b):
        """Per-mode (Y_air, Y_sub, rho) for 'te' and 'tm' as dicts."""
        w = self.omega
        eps = self.substrate.eps_complex
        bb = b / eps
        return {
            "te": (a / (w * MU0), b / (w * MU0), (b - a) / (b + a)),
            "tm": (w * EPS0 / a, w * EPS0 * eps / b, (a - bb) / (a + bb)),
        }

    def _regions(self, z):
        """Clipped local coordinates and masks for the three regions."""
        h = self.substrate.thickness
        return (
            np.maximum(z, 0.0),  # z_up  >= 0
            np.clip(z, -h, 0.0),  # z_sub in [-h, 0]
            np.minimum(z, -h) + h,  # z_dn  <= 0  (distance below the bottom face)
        )

    def _bottom_reflection(self, rho_natural):
        """The bottom-face (z=-h) load reflection coefficient actually used: the natural
        air-load value (unchanged) for ``backing="air"``, or -1 (short circuit) for "pec". Every
        field formula below is written purely in terms of this one reflection coefficient, so a
        PEC ground needs no separate "z < -h is empty" code path: with rho = -1, the (1+rho)
        factor that already multiplies every "z < -h" field is exactly zero (GREEN_FUNCTION.md
        Sec. 12A), and the substrate standing wave built from it satisfies V(-h) = 0 exactly.
        """
        return rho_natural if self.backing == "air" else -1.0

    def _sheet_mode_transfer(self, a, b, z, sheet_side):
        """{'te': (V, I), 'tm': (V, I)} of the sheet field per unit source J_e."""
        h = self.substrate.thickness
        z_te, z_tm = self._impedances_from_kz(a, b)
        consts = self._mode_constants(a, b)
        z_up, z_sub, z_dn = self._regions(z)
        u = np.exp(-2j * b * h)
        above = (z > 0.0) | ((z == 0.0) & (sheet_side == "above"))
        below = z < -h
        out = {}
        for mode, zeq in (("te", z_te), ("tm", z_tm)):
            y_a, y_s, rho_natural = consts[mode]
            rho = self._bottom_reflection(rho_natural)
            v0 = -zeq
            den = 1.0 + rho * u
            v_up = v0 * np.exp(-1j * a * z_up)
            i_up = y_a * v_up
            e1 = np.exp(1j * b * z_sub)
            e2 = np.exp(-1j * b * (2.0 * h + z_sub))
            v_sub = v0 * (e1 + rho * e2) / den
            i_sub = v0 * y_s * (rho * e2 - e1) / den
            v_dn = v0 * (1.0 + rho) * np.exp(-1j * b * h) / den * np.exp(1j * a * z_dn)
            i_dn = -y_a * v_dn
            out[mode] = (
                np.where(above, v_up, np.where(below, v_dn, v_sub)),
                np.where(above, i_up, np.where(below, i_dn, i_sub)),
            )
        return out

    def _stack_mode_transfer(self, a, b, z, port):
        """{'te': (V, I), 'tm': (V, I)} of the unpatterned stack for unit incident mode voltage.

        Port 'top': incident wave travels -z with V = 1 at z = 0.
        Port 'bottom' (``backing="air"`` only; the caller rejects it otherwise): incident wave
        travels +z with V = 1 at z = -h (mirror image, z -> -h - z, I -> -I, of the top-port
        solution; the unpatterned *ungrounded* stack is symmetric).

        Top-of-stack reflection Gamma, generalized to an arbitrary bottom-face load reflection
        ``rho`` (GREEN_FUNCTION.md Sec. 12A; ``r`` is the fixed, backing-independent air/substrate
        interface reflection):

            Gamma       = (r + rho*u) / (1 + rho*u*r)
            1 + Gamma   = (1 + r)*(1 + rho*u) / (1 + rho*u*r)

        which is exactly the existing ``r*(1-u)/(1-r^2*u)`` / ``(1+r)*(1-r*u)/(1-r^2*u)`` formula
        when rho = -r (the natural, symmetric ungrounded case: the bottom load is air too), and a
        new, independently-derived formula when rho = -1 (PEC ground).
        """
        h = self.substrate.thickness
        if port == "bottom":
            z = -h - z
        consts = self._mode_constants(a, b)
        z_up, z_sub, z_dn = self._regions(z)
        u = np.exp(-2j * b * h)
        above = z >= 0.0
        below = z < -h
        out = {}
        for mode in ("te", "tm"):
            y_a, y_s, rho_natural = consts[mode]
            r = -rho_natural
            rho = self._bottom_reflection(rho_natural)
            gamma = (r + rho * u) / (1.0 + rho * u * r)
            one_plus_gamma = (1.0 + r) * (1.0 + rho * u) / (1.0 + rho * u * r)
            w_sub = one_plus_gamma / (1.0 + rho * u)
            inc = np.exp(1j * a * z_up)
            ref = gamma * np.exp(-1j * a * z_up)
            v_up, i_up = inc + ref, y_a * (ref - inc)
            e1 = np.exp(1j * b * z_sub)
            e2 = np.exp(-1j * b * (2.0 * h + z_sub))
            v_sub, i_sub = w_sub * (e1 + rho * e2), w_sub * y_s * (rho * e2 - e1)
            v_dn = w_sub * (1.0 + rho) * np.exp(-1j * b * h) * np.exp(1j * a * z_dn)
            i_dn = -y_a * v_dn
            v = np.where(above, v_up, np.where(below, v_dn, v_sub))
            i = np.where(above, i_up, np.where(below, i_dn, i_sub))
            out[mode] = (v, -i) if port == "bottom" else (v, i)
        return out

    def _apply_modes(self, kx, ky, vec, modes):
        """Combine per-mode (V, I) into tangential vectors for an in-plane input ``vec``.

        E_t = V_TM (vec.p) p + V_TE (vec.q) q,   H_t = -I_TE (vec.q) p + I_TM (vec.p) q,
        p = k_t/|k_t|, q = z x p. At kt = 0 (|kt| <= 1e-12 k0) p, q are arbitrary and TE = TM,
        so the two modes are replaced by their mean: E_t = V vec, H_t = I z x vec exactly.
        """
        p, q = te_tm_unit_vectors(kx, ky, kt_floor=_KT_ZERO_RTOL * self.k0)
        (v_te, i_te), (v_tm, i_tm) = modes["te"], modes["tm"]
        near0 = np.hypot(kx, ky) <= _KT_ZERO_RTOL * self.k0
        v_avg, i_avg = 0.5 * (v_te + v_tm), 0.5 * (i_te + i_tm)
        v_te, v_tm = np.where(near0, v_avg, v_te), np.where(near0, v_avg, v_tm)
        i_te, i_tm = np.where(near0, i_avg, i_te), np.where(near0, i_avg, i_tm)
        vec = np.broadcast_to(np.asarray(vec, dtype=complex), p.shape)
        vp = np.sum(vec * p, axis=-1)[..., None]
        vq = np.sum(vec * q, axis=-1)[..., None]
        e_t = v_tm[..., None] * vp * p + v_te[..., None] * vq * q
        h_t = -i_te[..., None] * vq * p + i_tm[..., None] * vp * q
        return e_t, h_t

    def sheet_fields(self, kx, ky, current, z, sheet_side="above"):
        """Tangential (E_t, H_t) at height ``z`` radiated by the current sheet.

        ``current`` is the Floquet coefficient J~(kx, ky) (shape (..., 2), A/m); returns two
        arrays of shape ``broadcast(kx, ky, z).shape + (2,)`` (V/m and A/m), the harmonic's
        amplitude at (x, y) = 0 (multiply by exp(-j(kx x + ky y)) for the field).

        E_t is continuous everywhere; H_t jumps by z x [H(0+) - H(0-)] = J on the sheet, so at
        exactly z = 0 ``sheet_side`` ('above': limit z -> 0+, 'below': z -> 0-) selects the side.
        At z = 0 E_t equals ``spectral_dyad(kx, ky) @ current``. ``backing="pec"`` (Phase 12A):
        both fields are exactly zero for z < -h (nothing exists below the ground) and E_t(-h) = 0.
        """
        if sheet_side not in ("above", "below"):
            raise ValueError(f"sheet_side must be 'above' or 'below', got {sheet_side!r}")
        kx, ky, z = (np.asarray(v, dtype=float) for v in np.broadcast_arrays(kx, ky, z))
        a, b = self._fields_kz(kx, ky)
        modes = self._sheet_mode_transfer(a, b, z, sheet_side)
        return self._apply_modes(kx, ky, current, modes)

    def stack_fields(self, kx, ky, e_incident, z, port="top"):
        """Tangential (E_t, H_t) of the *unpatterned* stack at height ``z`` for a unit plane wave.

        ``e_incident`` (shape (..., 2)) is the incident tangential E at the port's reference
        interface: port 'top' (wave travelling -z, from z > 0) referenced at z = 0; port 'bottom'
        (wave travelling +z, from z < -h) referenced at z = -h. The incident wave must propagate
        (|k_t| < k0). Fields include the incident wave and the stack's reflected/transmitted waves.

        ``backing="pec"``: there is nothing below the ground, so port 'bottom' does not exist and
        raises ``ValueError`` (Phase 12A; a real, physical restriction, not a missing feature).
        """
        if port not in ("top", "bottom"):
            raise ValueError(f"port must be 'top' or 'bottom', got {port!r}")
        if port == "bottom" and self.backing == "pec":
            raise ValueError(
                "no bottom port exists for a grounded stack (backing='pec'): the PEC ground at "
                "z = -h blocks any wave from below"
            )
        kx, ky, z = (np.asarray(v, dtype=float) for v in np.broadcast_arrays(kx, ky, z))
        if np.any(np.hypot(kx, ky) >= self.k0):
            raise ValueError("the incident wave must propagate in air (|k_t| < k0)")
        a, b = self._fields_kz(kx, ky)
        modes = self._stack_mode_transfer(a, b, z, port)
        return self._apply_modes(kx, ky, e_incident, modes)
