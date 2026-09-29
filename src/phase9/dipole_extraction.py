"""Sec. 20: the "local field metric" Phase 9B compares between the isolated and embedded cases --
here, the induced-current dipole moment of one sub-cell's footprint, reconstructed from the
already-solved rooftop coefficients via `basis_functions.RooftopBasis.evaluate` (Phase 2,
unmodified). This is a standard EM quantity (the zeroth-order/DC spatial moment of the induced
surface current), not an invented one: it is exactly the quantity that determines a periodic
cell's own specular (0, 0) Floquet scattering coefficient (`basis_functions.RooftopBasis.
current_coefficients` computes the same integral at kx=ky=0 for a *whole*-cell, single-patch
basis; this module applies the identical integral restricted to one sub-cell's footprint of a
larger composite basis).
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DipoleMoment:
    px: complex
    py: complex


def local_dipole_moment(
    basis, coeffs, region_center_x, region_center_y, half_width_x, half_width_y, kx0=0.0, ky0=0.0, n_samples=48
):
    """Integral of J(x, y) over the axis-aligned rectangle
    [cx-hw, cx+hw] x [cy-hw, cy+hw], via `basis.evaluate` (reused, unmodified) on an
    ``n_samples`` x ``n_samples`` grid (a Riemann-sum approximation, refined by
    `test_dipole_extraction_matches_current_coefficients_at_kx_ky_zero`'s cross-check against the
    basis's own closed-form Fourier-transform integral)."""
    lin_x = (
        region_center_x + np.linspace(-half_width_x, half_width_x, n_samples, endpoint=False) + half_width_x / n_samples
    )
    lin_y = (
        region_center_y + np.linspace(-half_width_y, half_width_y, n_samples, endpoint=False) + half_width_y / n_samples
    )
    x, y = np.meshgrid(lin_x, lin_y, indexing="xy")
    j = basis.evaluate(x, y, coeffs, kx0, ky0)
    d_area = (2 * half_width_x / n_samples) * (2 * half_width_y / n_samples)
    px = complex(np.sum(j[..., 0]) * d_area)
    py = complex(np.sum(j[..., 1]) * d_area)
    return DipoleMoment(px, py)
