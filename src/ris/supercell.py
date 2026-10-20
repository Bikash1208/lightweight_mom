"""Supercell MoM with per-site diode states and Floquet-order extraction (Phase 20 R3 item 4: locality test). MODEL_ONLY.

A cell of `src/ris/cells.py` has a state-independent metal pattern; only the diode loads change with the state. A supercell of `nx x ny` such
cells with a state per site therefore has ONE impedance matrix Z0 (the metal pattern) and one sparse load contribution per coding. `solve_orders`
assembles Z0 once for an incidence direction and returns, for every coding, the reflection amplitude of every propagating Floquet order
(m, n) of the SUPERCELL lattice as a power-wave 2x2 matrix in the (TE, TM) basis of that order's own plane (TM = tangential E along the order's
transverse direction p_mn, TE = z x p_mn; the (0,0) order is exactly `compute_s_parameters_oblique`'s `s11`).

Order extraction: with the solved rooftop coefficients I, the Floquet current harmonic of order (m, n) is
`basis.current_coefficients(I, kx_mn, ky_mn)` (the transform is evaluated at that order's wavenumber), its radiated field above the sheet is
`green.sheet_fields`, and `decompose_oblique_fields` at (theta_mn, phi_mn) gives the up-going amplitude; the power-wave normalisation is
sqrt(Y_out / Y_in) with Y the TE/TM wave admittances at the respective angles (as in `sparameters._oblique_scattering`, whose (0,0) result this
reproduces; tested). Evanescent orders are skipped.

Memory: Z0 is dense (n^2 complex); each coding adds its loads in place at the few touched entries, solves on a copy, and restores the entries.
"""

import math
from dataclasses import replace

import numpy as np

from ..basis_functions import MeshSettings, RooftopBasis
from ..composite_geometry import AxisRect
from ..green_function import LayeredGreenFunction
from ..loads import _add_region, _add_rotated
from ..mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from ..periodic_bc import RectangularLattice
from ..polygon_geometry import Poly, PolyComposite
from ..sparameters import (
    _check_planes,
    _mode_admittances,
    decompose_oblique_fields,
    oblique_basis_vectors,
)


# ----------------------------------------------------------------------------------------------------------------------- geometry
def _translate_item(item, dx, dy):
    if isinstance(item, Poly):
        return Poly(tuple((x + dx, y + dy) for x, y in item.vertices))
    if isinstance(item, AxisRect):
        return AxisRect(item.cx + dx, item.cy + dy, item.length, item.width)
    raise TypeError(f"cannot translate {type(item).__name__}")


def build_supercell(cell, state_grid, kit):
    """`UnitCell` of `ny x nx` copies of `cell` (a `Cp2Cell`/`Lp*Cell`) with diode state `state_grid[j][i]` at site (i, j); (0, 0) is the
    lower-left site. The supercell origin is its centre; sub-cell (i, j) is centred at ((i - (nx-1)/2) P, (j - (ny-1)/2) P)."""
    ny, nx = len(state_grid), len(state_grid[0])
    p = cell.period
    template = cell.unit_cell(state_grid[0][0], kit)
    metal = template.geometry.metal
    apertures = getattr(template.geometry, "apertures", ())
    new_metal, new_apertures, loads = [], [], []
    for j in range(ny):
        for i in range(nx):
            dx, dy = (i - (nx - 1) / 2) * p, (j - (ny - 1) / 2) * p
            new_metal += [_translate_item(m, dx, dy) for m in metal]
            new_apertures += [_translate_item(a, dx, dy) for a in apertures]
            loads += [replace(ld, x=ld.x + dx, y=ld.y + dy) for ld in cell.unit_cell(state_grid[j][i], kit).loads]
    geometry = PolyComposite(tuple(new_metal), tuple(new_apertures))
    return replace(template, lattice=RectangularLattice(nx * p, ny * p), geometry=geometry, loads=tuple(loads))


# ----------------------------------------------------------------------------------------------------------------------- sparse loads
class _Accumulator:
    """Dict-backed stand-in for the dense `out` matrix of `loads._add_region` / `_add_rotated` (they only do `out[i, j] += v` and
    `out[np.ix_(a, b)] += block`)."""

    def __init__(self):
        self.d = {}

    def __getitem__(self, key):
        r, c = key
        if isinstance(r, np.ndarray) or isinstance(c, np.ndarray):
            rr, cc = np.broadcast_arrays(np.asarray(r), np.asarray(c))
            out = np.zeros(rr.shape, dtype=complex)
            for idx in np.ndindex(rr.shape):
                out[idx] = self.d.get((int(rr[idx]), int(cc[idx])), 0.0)
            return out
        return self.d.get((int(r), int(c)), 0.0)

    def __setitem__(self, key, value):
        r, c = key
        if isinstance(r, np.ndarray) or isinstance(c, np.ndarray):
            rr, cc = np.broadcast_arrays(np.asarray(r), np.asarray(c))
            val = np.broadcast_to(np.asarray(value), rr.shape)
            for idx in np.ndindex(rr.shape):
                self.d[(int(rr[idx]), int(cc[idx]))] = complex(val[idx])
        else:
            self.d[(int(r), int(c))] = complex(value)


def sparse_load_matrix(basis, loads, omega):
    """Same entries as `loads.load_matrix` as a dict {(row, col): value} (no n x n dense array)."""
    acc = _Accumulator()
    for ld in loads:
        if ld.is_rotated:
            _add_rotated(acc, basis, ld, ld.impedance(omega) * ld.sheet_factor)
            continue
        region = ld.region()
        _add_region(acc, basis, ld.direction, ld.impedance(omega) * ld.sheet_factor, region)
        if ld.open_circuit or not ld.transverse_conducts:
            _add_region(acc, basis, "y" if ld.direction == "x" else "x", complex(1e12), region)
    return acc.d


# ----------------------------------------------------------------------------------------------------------------------- solve
def solve_orders(unit_cell, frequency, mesh, planes, theta, phi, load_sets, max_order=3):
    """Per load set (tuple of `LumpedLoad`): `{(m, n): S}` for every propagating order |m|, |n| <= max_order of the cell's own lattice;
    S is the 2x2 power-wave matrix `[out (TE, TM), in (TE, TM)]` (the (0,0) entry equals `compute_s_parameters_oblique(...).s11`).
    Z0 is assembled once; the loads are added per set (`unit_cell.loads` is ignored)."""
    if not 0.0 <= theta < math.pi / 2:
        raise ValueError(f"theta must satisfy 0 <= theta < pi/2, got {theta}")
    backing = getattr(unit_cell, "backing", "air")
    if backing != "pec":
        raise ValueError("solve_orders handles the grounded (one-port) stack only")
    green = LayeredGreenFunction(unit_cell.substrate, frequency, backing=backing)
    _check_planes(green, planes)
    k0 = green.k0
    kx0, ky0 = k0 * math.sin(theta) * math.cos(phi), k0 * math.sin(theta) * math.sin(phi)
    basis = RooftopBasis(unit_cell, mesh.cells_x, mesh.cells_y)
    z0 = assemble_impedance_matrix(basis, green, kx0, ky0, mesh.floquet_order)
    q, p = oblique_basis_vectors(phi)
    e_inc = (q, p)
    v = np.empty((basis.n_unknowns, 2), dtype=complex)
    for mode in range(2):
        e_stack = green.stack_fields(kx0, ky0, e_inc[mode].astype(complex), 0.0, "top")[0]
        v[:, mode] = assemble_excitation_vector(basis, e_stack, kx0, ky0)
    z_planes = np.array([planes.z_top])
    y_in = np.array(_mode_admittances(frequency, theta))  # (TE, TM)
    lat = unit_cell.lattice
    orders = []
    for m in range(-max_order, max_order + 1):
        for n in range(-max_order, max_order + 1):
            kx, ky = lat.transverse_wavenumbers(kx0, ky0, np.array(m), np.array(n))
            kt = math.hypot(float(kx), float(ky))
            if kt < k0 * (1.0 - 1e-9):
                orders.append((m, n, float(kx), float(ky), kt))
    results = []
    for load_set in load_sets:
        entries = sparse_load_matrix(basis, tuple(load_set), green.omega) if load_set else {}
        saved = {k: z0[k] for k in entries}
        for k, val in entries.items():
            z0[k] += val
        try:
            currents = np.linalg.solve(z0, v)
        finally:
            for k, val in saved.items():
                z0[k] = val
        if not np.all(np.isfinite(currents)):
            raise np.linalg.LinAlgError("solution contains non-finite values")
        a_modes = np.zeros(2, dtype=complex)  # incident (down-going) amplitude at the plane, per mode
        for mode in range(2):
            e_stack, h_stack = green.stack_fields(kx0, ky0, e_inc[mode].astype(complex), z_planes, "top")
            e_sheet, h_sheet = green.sheet_fields(
                kx0, ky0, basis.current_coefficients(currents[:, mode], kx0, ky0), z_planes
            )
            _, down = decompose_oblique_fields(e_stack + e_sheet, h_stack + h_sheet, phi, *y_in)
            a_modes[mode] = down[0, mode]
        out = {}
        for m, n, kx, ky, kt in orders:
            theta_o = math.asin(min(kt / k0, 1.0))
            phi_o = math.atan2(ky, kx) if kt > 0.0 else phi
            y_te, y_tm = _mode_admittances(frequency, theta_o)
            y_out = np.array([y_te, y_tm])
            s = np.zeros((2, 2), dtype=complex)
            for mode in range(2):
                e_sheet, h_sheet = green.sheet_fields(
                    kx, ky, basis.current_coefficients(currents[:, mode], kx, ky), z_planes
                )
                if (m, n) == (0, 0):
                    e_stack, h_stack = green.stack_fields(kx0, ky0, e_inc[mode].astype(complex), z_planes, "top")
                    e_sheet, h_sheet = e_stack + e_sheet, h_stack + h_sheet
                up, _ = decompose_oblique_fields(e_sheet, h_sheet, phi_o, y_te, y_tm)
                s[:, mode] = up[0] / a_modes[mode] * np.sqrt(y_out / y_in[mode])
            out[(m, n)] = s
        results.append(out)
    return results


def default_mesh(n_per_cell, nx, ny):
    """`MeshSettings` for an `nx x ny` supercell of cells meshed with `n_per_cell` cells per period: same pitch, Floquet cutoff
    `n_total // 2 + 1` (the isolated-cell convention applied to the supercell mesh; the rank condition needs 2F + 1 >= cells)."""
    return MeshSettings(n_per_cell * nx, n_per_cell * ny, n_per_cell * max(nx, ny) // 2 + 1)
