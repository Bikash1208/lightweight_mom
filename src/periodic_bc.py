"""Rectangular lattice and Floquet-mode bookkeeping.

Convention (time factor exp(+j*omega*t)): a plane wave varies as
exp(-j*k.r), so a field quasi-periodic under a lattice translation obeys

    E(x + Px, y) = E(x, y) * exp(-j*kx0*Px)
    E(x, y + Py) = E(x, y) * exp(-j*ky0*Py)

where (kx0, ky0) is the transverse wavevector of the incident wave. The
Floquet harmonics are exp(-j*(kxm*x + kyn*y)) with

    kxm = kx0 + 2*pi*m/Px,    kyn = ky0 + 2*pi*n/Py,    m, n integers.

Only lattice bookkeeping lives here. No Green's function is implemented.
"""

import math
from dataclasses import dataclass

import numpy as np

from ._validate import finite_real


@dataclass(frozen=True)
class RectangularLattice:
    """Rectangular periodic lattice, periods along x and y [m]."""

    period_x: float
    period_y: float

    def __post_init__(self):
        for name in ("period_x", "period_y"):
            v = finite_real(name, getattr(self, name))
            if v <= 0.0:
                raise ValueError(f"{name} must be > 0 m, got {v!r}")
            object.__setattr__(self, name, v)

    @property
    def cell_area(self):
        return self.period_x * self.period_y

    @property
    def reciprocal_x(self):
        """2*pi/Px [rad/m]."""
        return 2.0 * math.pi / self.period_x

    @property
    def reciprocal_y(self):
        """2*pi/Py [rad/m]."""
        return 2.0 * math.pi / self.period_y

    def bloch_phase(self, kx0, ky0):
        """Phase factors (exp(-j*kx0*Px), exp(-j*ky0*Py)) picked up per cell."""
        return (
            np.exp(-1j * kx0 * self.period_x),
            np.exp(-1j * ky0 * self.period_y),
        )

    def transverse_wavenumbers(self, kx0, ky0, m, n):
        """Floquet transverse wavenumbers (kxm, kyn) for integer index arrays."""
        m = np.asarray(m)
        n = np.asarray(n)
        if not (np.issubdtype(m.dtype, np.integer) and np.issubdtype(n.dtype, np.integer)):
            raise ValueError("Floquet indices m, n must be integers")
        return kx0 + self.reciprocal_x * m, ky0 + self.reciprocal_y * n

    def mode_grid(self, max_order):
        """Integer index grids (m, n), each of shape (2N+1, 2N+1), indexing='ij'.

        m, n range over -N..N with N = ``max_order``.
        """
        if not isinstance(max_order, (int, np.integer)) or max_order < 0:
            raise ValueError(f"max_order must be an integer >= 0, got {max_order!r}")
        idx = np.arange(-max_order, max_order + 1)
        return np.meshgrid(idx, idx, indexing="ij")


def normal_wavenumber(k, kx, ky):
    """kz = sqrt(k^2 - kx^2 - ky^2) on the physical branch.

    ``k`` may be complex (lossy medium, Im(k) <= 0); ``kx``, ``ky`` are real.
    Branch choice, required by exp(+j*omega*t) and waves exp(-j*kz*|z|):

        Im(kz) <= 0     (decaying / outgoing),   Re(kz) >= 0 when propagating.

    * propagating, lossless: kz = +sqrt(k^2 - kt^2)  (real, > 0)
    * evanescent, lossless : kz = -j*sqrt(kt^2 - k^2)
    * lossy                : principal root already has Im < 0, Re > 0
    """
    k = np.asarray(k, dtype=complex)
    kt2 = np.asarray(kx, dtype=float) ** 2 + np.asarray(ky, dtype=float) ** 2
    if np.any(k.imag > 0.0):
        raise ValueError("Im(k) must be <= 0 under exp(+j*omega*t) (passive medium)")
    kz = np.sqrt(k**2 - kt2 + 0j)
    kz = np.where(kz.imag > 0.0, -kz, kz)
    return kz if kz.ndim else complex(kz)
