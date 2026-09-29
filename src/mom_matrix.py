"""Galerkin MoM assembly: impedance matrix Z and excitation vector V.

Solve Z I = V for the rooftop coefficients I [A/m]. Exact equations,
normalization and sign conventions: MOM_FORMULATION.md. Summary:

    Z_ij = -(1/A) * sum_{m,n} conj(Phi_i(k_mn)) * G_{d_i d_j}(k_mn) * Phi_j(k_mn)
    V_i  =  conj(Phi_i(k_00)) * E_stack,d_i

with A = Px*Py, Phi_b the transform of rooftop b (basis_functions.py), G the
layered Green's dyad (green_function.py), d_b in {x, y} the rooftop direction,
k_mn = (kx0 + 2 pi m/Px, ky0 + 2 pi n/Py), |m|,|n| <= floquet_order.
"""

import numpy as np

from .basis_functions import RooftopBasis
from .green_function import LayeredGreenFunction

# Upper bound on the size (elements) of the temporary (n_unknowns x n_modes_chunk)
# complex arrays, to keep memory ~ tens of MB regardless of the Floquet order.
_CHUNK_ELEMENTS = 2_000_000


def assemble_impedance_matrix(
    basis: RooftopBasis, green: LayeredGreenFunction, kx0, ky0, floquet_order
):
    """Z (n_unknowns x n_unknowns, complex), units ohm*m^2 (see MOM_FORMULATION.md)."""
    kx, ky, g = green.floquet_dyads(basis.lattice, kx0, ky0, floquet_order)
    kx, ky = kx.ravel(), ky.ravel()
    gxx, gxy, gyy = g[..., 0, 0].ravel(), g[..., 0, 1].ravel(), g[..., 1, 1].ravel()

    nb, nxb = basis.n_unknowns, basis.n_x
    z = np.zeros((nb, nb), dtype=complex)
    chunk = max(1, _CHUNK_ELEMENTS // nb)
    for s in range(0, kx.size, chunk):
        sl = slice(s, s + chunk)
        phi = basis.fourier_transform(kx[sl], ky[sl])          # (nb, kc)
        phi_x, phi_y = phi[:nxb], phi[nxb:]
        cx, cy = np.conj(phi_x), np.conj(phi_y)
        z[:nxb, :nxb] += (cx * gxx[sl]) @ phi_x.T
        z[:nxb, nxb:] += (cx * gxy[sl]) @ phi_y.T
        z[nxb:, :nxb] += (cy * gxy[sl]) @ phi_x.T
        z[nxb:, nxb:] += (cy * gyy[sl]) @ phi_y.T
    return -z / basis.lattice.cell_area


def assemble_excitation_vector(basis: RooftopBasis, e_stack_t, kx0, ky0):
    """V_i = conj(Phi_i(kx0, ky0)) * E_stack,d_i  (units V*m)."""
    e_stack_t = np.asarray(e_stack_t, dtype=complex)
    if e_stack_t.shape != (2,):
        raise ValueError(f"e_stack_t must have shape (2,), got {e_stack_t.shape}")
    phi = basis.fourier_transform(kx0, ky0)
    return np.conj(phi) * e_stack_t[basis.direction]
