"""Rooftop basis functions on a uniform rectangular mesh of one periodic unit cell.

Deliberately minimal: x- and y-directed rooftops only. Derivation and
normalization: MOM_FORMULATION.md.

Mesh. The cell [-Px/2, Px/2] x [-Py/2, Py/2] is divided into ``nx`` x ``ny``
rectangles of size dx = Px/nx, dy = Py/ny. A rectangle is *metal* when its
centroid is metal (``UnitCell.is_metal``, which already returns the complement
for ``slot``: slots are represented by metal-complement currents). Features whose
edges lie on mesh lines are represented exactly; others are staircased.

Rooftops. One unknown per mesh edge whose two neighbouring rectangles are both
metal (so the current component normal to a metal/non-metal boundary is zero by
construction: the PEC boundary condition on J):

* x-directed rooftop on the edge x = x_e between rectangles (i-1, j) and (i, j):
      f(x, y) = tri((x - x_e)/dx) * pulse((y - y_c)/dy)  * x_hat,
  tri(t) = max(0, 1 - |t|), pulse(t) = 1 for -1/2 <= t < 1/2. Peak value 1 [A/m per
  unit coefficient] at the edge, support 2dx x dy.
* y-directed rooftop: same with x and y exchanged.

Periodicity. Edge indices wrap: the edge at x = -Px/2 is the same unknown as the
one at x = +Px/2 (likewise in y). Its Bloch-periodic extension continues the
current into neighbouring cells (needed for wires and for slot complements,
whose metal touches the cell boundary). The unknown is referenced to the copy
inside the cell, at x = -Px/2 (or y = -Py/2).

Transform convention (CONVENTIONS.md): F(kx, ky) = integral of f exp(+j k.r) dA over
the *unwrapped* rooftop, which equals the transform of its whole Bloch extension
over one cell. Closed form (np.sinc(z) = sin(pi z)/(pi z)):

    F_x = exp(j(kx xc + ky yc)) * dx sinc^2(kx dx/2pi) * dy sinc(ky dy/2pi)
    F_y = exp(j(kx xc + ky yc)) * dx sinc(kx dx/2pi)   * dy sinc^2(ky dy/2pi)
"""

from dataclasses import dataclass

import numpy as np

from .geometry import UnitCell


def _check_integer(name, v, minimum):
    if isinstance(v, bool) or not isinstance(v, (int, np.integer)):
        raise ValueError(f"{name} must be an integer, got {v!r}")
    if v < minimum:
        raise ValueError(f"{name} must be >= {minimum}, got {v}")


@dataclass(frozen=True)
class MeshSettings:
    """Discretization controls (numerical, not physical).

    Rank condition: rooftop positions repeat with the mesh, so Floquet harmonics m and
    m + cells_x (n and n + cells_y) are indistinguishable at the rooftop centres. The
    truncated harmonic set must therefore contain every residue class,
    2*floquet_order + 1 >= max(cells_x, cells_y); otherwise Z is rank-deficient and the
    solution non-unique. (This is only the *minimum*: accuracy needs a several times
    larger order, see MOM_FORMULATION.md.)
    """

    cells_x: int
    cells_y: int
    floquet_order: int  # Floquet harmonics |m|, |n| <= floquet_order in the Z assembly

    def __post_init__(self):
        _check_integer("cells_x", self.cells_x, 2)
        _check_integer("cells_y", self.cells_y, 2)
        _check_integer("floquet_order", self.floquet_order, 0)
        if 2 * self.floquet_order + 1 < max(self.cells_x, self.cells_y):
            raise ValueError(
                f"floquet_order={self.floquet_order} is too small for the mesh "
                f"({self.cells_x} x {self.cells_y} cells): need 2*floquet_order + 1 >= "
                f"{max(self.cells_x, self.cells_y)} or Z is rank-deficient"
            )


class RooftopBasis:
    """x/y rooftop basis on the metal of ``unit_cell``; unknown ordering: all
    x-directed (row-major in edge index [i, j]) then all y-directed."""

    def __init__(self, unit_cell: UnitCell, nx: int, ny: int):
        _check_integer("nx", nx, 2)
        _check_integer("ny", ny, 2)
        self.unit_cell = unit_cell
        self.lattice = unit_cell.lattice
        self.nx, self.ny = int(nx), int(ny)
        px, py = self.lattice.period_x, self.lattice.period_y
        self.dx, self.dy = px / self.nx, py / self.ny

        cell_x = -px / 2 + (np.arange(self.nx) + 0.5) * self.dx
        cell_y = -py / 2 + (np.arange(self.ny) + 0.5) * self.dy
        gx, gy = np.meshgrid(cell_x, cell_y, indexing="ij")
        self.metal_mask = np.asarray(unit_cell.is_metal(gx, gy), dtype=bool)  # [i, j]

        # x-edge (i, j): between rectangles (i-1, j) and (i, j); y-edge: (i, j-1) and (i, j)
        active_x = self.metal_mask & np.roll(self.metal_mask, 1, axis=0)
        active_y = self.metal_mask & np.roll(self.metal_mask, 1, axis=1)
        ix, jx = np.nonzero(active_x)
        iy, jy = np.nonzero(active_y)
        self.n_x, self.n_y = ix.size, iy.size
        if self.n_x + self.n_y == 0:
            raise ValueError(f"mesh has no interior metal edges (metal too small for this mesh): nx={nx}, ny={ny}")
        self.edge_index = np.concatenate([np.stack([ix, jx], 1), np.stack([iy, jy], 1)])
        self.direction = np.concatenate([np.zeros(self.n_x, int), np.ones(self.n_y, int)])
        self.center_x = np.concatenate([-px / 2 + ix * self.dx, cell_x[iy]])
        self.center_y = np.concatenate([cell_y[jx], -py / 2 + jy * self.dy])

    @property
    def n_unknowns(self):
        return self.n_x + self.n_y

    # ------------------------------------------------------------------ transform
    def fourier_transform(self, kx, ky):
        """Phi_b(kx, ky) = integral f_b exp(+j k.r) dA. Shape (n_unknowns,) + shape(kx)."""
        kx, ky = np.broadcast_arrays(np.asarray(kx, float), np.asarray(ky, float))
        shape = kx.shape
        kxf, kyf = kx.ravel(), ky.ravel()
        sx = np.sinc(kxf * self.dx / (2 * np.pi))
        sy = np.sinc(kyf * self.dy / (2 * np.pi))
        amp = np.empty((self.n_unknowns, kxf.size))
        amp[: self.n_x] = (self.dx * self.dy) * (sx**2 * sy)[None, :]
        amp[self.n_x :] = (self.dx * self.dy) * (sx * sy**2)[None, :]
        phase = np.exp(1j * (np.outer(self.center_x, kxf) + np.outer(self.center_y, kyf)))
        return (amp * phase).reshape((self.n_unknowns,) + shape)

    def current_coefficients(self, coeffs, kx, ky):
        """Floquet coefficients J~(kx, ky) = (1/A) sum_b I_b Phi_b, shape shape(kx)+(2,)."""
        coeffs = np.asarray(coeffs, dtype=complex)
        if coeffs.shape != (self.n_unknowns,):
            raise ValueError(f"expected {self.n_unknowns} coefficients, got shape {coeffs.shape}")
        phi = self.fourier_transform(kx, ky)
        area = self.lattice.cell_area
        jx = np.tensordot(coeffs[: self.n_x], phi[: self.n_x], axes=1) / area
        jy = np.tensordot(coeffs[self.n_x :], phi[self.n_x :], axes=1) / area
        return np.stack([jx, jy], axis=-1)

    # --------------------------------------------------------------- pointwise
    def evaluate(self, x, y, coeffs, kx0=0.0, ky0=0.0):
        """Bloch-periodic current J(x, y) = sum_b I_b f_b^Bloch, shape shape(x)+(2,).

        f^Bloch sums all lattice copies with phase exp(-j(kx0 p Px + ky0 q Py)) so
        J(x + Px, y) = J(x, y) exp(-j kx0 Px) for any point in the plane.
        """
        coeffs = np.asarray(coeffs, dtype=complex)
        if coeffs.shape != (self.n_unknowns,):
            raise ValueError(f"expected {self.n_unknowns} coefficients, got shape {coeffs.shape}")
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        shape = x.shape
        px, py = self.lattice.period_x, self.lattice.period_y
        # fold points into the reference cell; remember the lattice shift for the phase
        p = np.floor(x.ravel() / px + 0.5)
        q = np.floor(y.ravel() / py + 0.5)
        xf, yf = x.ravel() - p * px, y.ravel() - q * py
        bloch = np.exp(-1j * (kx0 * p * px + ky0 * q * py))

        jx = np.zeros(xf.size, complex)
        jy = np.zeros(xf.size, complex)
        used = np.flatnonzero(coeffs != 0.0)  # zero coefficients contribute nothing
        cx, cy = self.center_x[used, None], self.center_y[used, None]
        c_used, dir_used = coeffs[used], self.direction[used]
        is_x = (dir_used == 0)[:, None]
        for sp in (-1, 0, 1):  # basis copies that can overlap the reference cell
            for sq in (-1, 0, 1):
                w = np.exp(-1j * (kx0 * sp * px + ky0 * sq * py))
                u = (xf[None, :] - (cx + sp * px)) / self.dx
                v = (yf[None, :] - (cy + sq * py)) / self.dy
                tri_u = np.clip(1.0 - np.abs(u), 0.0, None)
                tri_v = np.clip(1.0 - np.abs(v), 0.0, None)
                pul_u = ((u >= -0.5) & (u < 0.5)).astype(float)  # pulse over one mesh width
                pul_v = ((v >= -0.5) & (v < 0.5)).astype(float)
                val = np.where(is_x, tri_u * pul_v, pul_u * tri_v)
                contrib = w * (c_used[:, None] * val)
                jx += contrib[dir_used == 0].sum(axis=0)
                jy += contrib[dir_used == 1].sum(axis=0)
        return (np.stack([jx, jy], axis=-1) * bloch[:, None]).reshape(shape + (2,))
