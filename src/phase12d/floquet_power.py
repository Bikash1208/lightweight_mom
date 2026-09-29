"""Task B: the aperture-level design metric -- reflected power in EVERY propagating Floquet
order of a (super)cell's own lattice, and the phase of the intended (anomalous-reflection) order,
computed directly from the solved surface current. This decides Phase 13 (brief Sec. "Why" point
2), not Phase 9's dipole-moment proxy.

Method. The solved rooftop current's Fourier coefficient at any Floquet harmonic (kx, ky) is
``basis_functions.RooftopBasis.current_coefficients`` (exact closed form, already used by every
S-parameter extraction in this project). Feeding it to
``green_function.LayeredGreenFunction.sheet_fields`` gives that harmonic's own tangential (E, H) in
the air above the sheet -- a genuine single-plane-wave field for that (kx, ky), because the
Floquet harmonics are mutually orthogonal over one lattice cell (Parseval: integrating the exact
cross term of two DIFFERENT harmonics over the cell is exactly zero). The time-average power this
harmonic alone carries through one cell is therefore

    P_mn = (A/2) * Re{ E_mn,x * conj(H_mn,y) - E_mn,y * conj(H_mn,x) },   A = cell area,

exactly as for a single plane wave, and the total reflected power is exactly ``sum_mn P_mn``
(no cross terms, by the same orthogonality) -- this is the "document the method" the brief asks
for. Order (0,0) additionally carries the unpatterned-stack (incident + Airy-reflected) field, so
it needs the same up/down split ``sparameters.decompose_air_fields`` already uses to isolate the
REFLECTED (up-going) part before applying the formula above; every other order has no incident
contribution (only the specular order is a port mode, CONVENTIONS.md) so its ``sheet_fields``
result at z > 0 is already purely outgoing.

Normalizing by the incident power (``A/(2*eta0)`` for a unit-amplitude normally-incident plane
wave) gives ``fraction_mn = P_mn / P_incident``; at (0,0) this is exactly ``|S11|**2`` by
construction (validated: ``tests/test_phase12d_floquet_power.py``), and ``sum_mn fraction_mn`` over
every order (propagating or not; evanescent orders contribute ~0 real power automatically, since a
reactive near-field's E x conj(H) is purely imaginary in a lossless medium) equals the lossy
grounded cell's own total reflected-power fraction.

Validation required by the brief:
- Uniform array (Sec. below, ``validate_uniform_array``): treating a homogeneous nxn supercell
  with an ASSUMED period n times its true period must put 100% of the power in the orders that are
  exact multiples of n (which must reproduce the true isolated |S11|**2 exactly) and ~0 in every
  other, spurious order the larger assumed period introduces.
- A trivial case with a known answer (``validate_full_pec_sheet``): a patch covering the WHOLE
  period has a spatially uniform (constant) current, whose Fourier transform is exactly zero at
  every harmonic except (0,0) -- so every higher order's power must be exactly (not just
  approximately) zero, and |S11| = 1 (full reflection, no loss possible below a solid PEC sheet).
"""

from dataclasses import dataclass

import numpy as np

from ..constants import ETA0
from ..sparameters import decompose_air_fields


@dataclass(frozen=True)
class FloquetPowerResult:
    m: np.ndarray                  # (2N+1, 2N+1) int
    n: np.ndarray                  # (2N+1, 2N+1) int
    kx: np.ndarray
    ky: np.ndarray
    propagating: np.ndarray        # bool
    fraction: np.ndarray           # P_mn / P_incident, real
    e_reflected_x: np.ndarray      # (2N+1, 2N+1) complex: reflected E_x Floquet amplitude per
                                    # order (the up-going field for order (0,0); phase of this is
                                    # the "phase of the intended order" the brief asks for)

    @property
    def total_reflected_fraction(self):
        return float(np.sum(self.fraction))

    @property
    def total_propagating_fraction(self):
        return float(np.sum(self.fraction[self.propagating]))

    def _index(self, m, n):
        idx = np.argwhere((self.m == m) & (self.n == n))
        if idx.size == 0:
            raise ValueError(f"order ({m}, {n}) not in this result's max_order range")
        return tuple(idx[0])

    def fraction_at(self, m, n):
        return float(self.fraction[self._index(m, n)])

    def phase_deg_at(self, m, n):
        """Phase (deg, principal value) of the reflected E_x Floquet amplitude of order (m, n)."""
        return float(np.degrees(np.angle(self.e_reflected_x[self._index(m, n)])))


def floquet_order_powers(basis, currents, green, lattice, max_order, incident_polarization="x"):
    """``FloquetPowerResult`` for every Floquet order |m|, |n| <= ``max_order`` of ``lattice``
    (the SUPERCELL's own lattice -- pass ``supercell.lattice``, whose period is n times the
    sub-cell period, so these orders are the aperture's own anomalous-reflection diffraction
    orders, not the sub-cell's)."""
    m, n = lattice.mode_grid(max_order)
    kx, ky = lattice.transverse_wavenumbers(0.0, 0.0, m, n)
    kz = green.normal_wavenumbers(kx, ky)[0]
    propagating = kz.imag == 0.0

    j = basis.current_coefficients(currents, kx, ky)   # (2N+1, 2N+1, 2)
    e, h = green.sheet_fields(kx, ky, j, 0.0, sheet_side="above")

    i0 = int(max_order)   # index of order (0, 0) in the (2*max_order+1) grid
    pol_vec = np.array([1.0, 0.0], dtype=complex) if incident_polarization == "x" else np.array([0.0, 1.0], dtype=complex)
    e_stack, h_stack = green.stack_fields(0.0, 0.0, pol_vec, 0.0, port="top")
    up00, _ = decompose_air_fields(e[i0, i0] + e_stack, h[i0, i0] + h_stack)
    h_up00 = np.array([-up00[1], up00[0]], dtype=complex) / ETA0   # H_up = z x E_up / eta0 (kt=0, exact)

    e = e.copy()
    h = h.copy()
    e[i0, i0] = up00
    h[i0, i0] = h_up00

    cross_z = (e[..., 0] * np.conj(h[..., 1]) - e[..., 1] * np.conj(h[..., 0])).real
    area = lattice.cell_area
    p_mn = 0.5 * area * cross_z
    p_incident = 0.5 * area / ETA0
    fraction = p_mn / p_incident
    return FloquetPowerResult(m, n, kx, ky, propagating, fraction, e[..., 0].copy())


def intended_order_for_gradient(n_cells, gradient_axis="x"):
    """The blazed-grating "intended" (anomalous-reflection) order for an N-cell linear phase-
    gradient supercell spanning one 2*pi cycle along ``gradient_axis``: order (+1, 0) or (0, +1)
    of the SUPERCELL's own lattice (whose period along that axis is N times the sub-cell period) --
    the standard generalized-Snell's-law blazed-grating result (one full 2*pi phase cycle per
    supercell period steers all of one polarization's reflected power into the first order)."""
    return (1, 0) if gradient_axis == "x" else (0, 1)
