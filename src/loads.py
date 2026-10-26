"""Lumped RLC loads across metal gaps (Phase 20 R1 item 2). Derivation and normalisation: MOM_FORMULATION.md Sec. 14.

A `LumpedLoad` is a small axis-aligned rectangle of the unit cell that the mesh treats as CONDUCTOR (rooftop basis functions cross
it) but whose conduction along `direction` carries a series impedance `Z_L(f) = R + j w L + 1/(j w C)` (`C=None`: no
capacitor term; `open_circuit=True`: `Z_L = OPEN_OHM`). Transverse to `direction` the region stays a perfect conductor (`transverse_conducts=True`, the default) unless
`transverse_conducts=False` or `open_circuit=True`, which remove the transverse conduction (an open circuit is then exactly the empty gap).
The load is the sheet impedance `Z_s = Z_L * width / length` (length along `direction`, width across it), entering the
Galerkin system as  (Z + Z_load) I = V,  Z_load,ij = Z_s * integral_region f_i . f_j dA  over the basis functions of `direction`.

No loads (the default, `UnitCell.loads == ()`) means none of this code runs: every pre-existing result is bit-identical.
"""

import math
import warnings
from dataclasses import dataclass

import numpy as np

from ._validate import finite_real

OPEN_OHM = 1.0e12  # "open circuit" series impedance [ohm]; the solver sees a very large but finite Z_L
_DIRECTIONS = ("x", "y")
_GAUSS2 = (np.array([-1.0, 1.0]) / np.sqrt(3.0), np.array([1.0, 1.0]))  # 2-point Gauss-Legendre: exact for quadratics


@dataclass(frozen=True)
class LumpedLoad:
    x: float  # centre of the region [m]
    y: float
    length: float  # extent ALONG `direction` [m]
    width: float  # extent ACROSS `direction` [m]
    direction: "str | float"  # "x" or "y": the conduction direction the impedance acts on; a number = angle [deg] from +x (Phase 20 R3: rotated load)
    R: float = 0.0  # series resistance [ohm] >= 0
    L: float = 0.0  # series inductance [H] >= 0
    C: "float | None" = None  # series capacitance [F] > 0, or None for no capacitor term
    open_circuit: bool = False  # Z_L = OPEN_OHM AND no conduction transverse to `direction`: the region is an empty gap
    transverse_conducts: bool = True  # False: no current transverse to `direction` in the region (a lumped element only conducts along its axis)
    z_of_f: "object" = None  # optional callable f[Hz] -> Z_L [ohm] replacing R, L, C (Phase 20 R2: diode models with a parallel branch)

    def __post_init__(self):
        for name in ("x", "y", "length", "width", "R", "L"):
            object.__setattr__(self, name, finite_real(name, getattr(self, name)))
        if self.C is not None:
            object.__setattr__(self, "C", finite_real("C", self.C))
            if self.C <= 0.0:
                raise ValueError(f"load C must be > 0 F (or None), got {self.C}")
        if isinstance(self.direction, int | float) and not isinstance(self.direction, bool):
            ang = float(finite_real("direction", self.direction)) % 180.0
            # an angle of 0 / 90 deg IS the axis-aligned load (same code path, bit-identical)
            object.__setattr__(self, "direction", "x" if ang == 0.0 else "y" if ang == 90.0 else ang)
        if self.direction not in _DIRECTIONS and not isinstance(self.direction, float):
            raise ValueError(f"load direction must be one of {_DIRECTIONS} or an angle in degrees, got {self.direction!r}")
        if self.length <= 0.0 or self.width <= 0.0:
            raise ValueError("load length and width must be > 0")
        if self.R < 0.0 or self.L < 0.0:
            raise ValueError("load R and L must be >= 0 (passive)")

    @property
    def is_rotated(self):
        """True for a conduction axis that is not x or y (Phase 20 R3): the region is a rectangle rotated by `direction` degrees."""
        return isinstance(self.direction, float)

    @property
    def angle(self):
        """Conduction direction as an angle [rad] from +x."""
        return {"x": 0.0, "y": 0.5 * math.pi}.get(self.direction, math.radians(self.direction) if self.is_rotated else 0.0)

    def contains(self, x, y):
        """Points inside the load region (rotated rectangle for a rotated load)."""
        if not self.is_rotated:
            x0, x1, y0, y1 = self.region()
            x, y = np.asarray(x, float), np.asarray(y, float)
            return (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1)
        c, s_ = math.cos(self.angle), math.sin(self.angle)
        dx, dy = np.asarray(x, float) - self.x, np.asarray(y, float) - self.y
        return (np.abs(c * dx + s_ * dy) <= self.length / 2) & (np.abs(-s_ * dx + c * dy) <= self.width / 2)

    def corners(self):
        """Four corners (x, y) of the region, counter-clockwise."""
        c, s_ = math.cos(self.angle), math.sin(self.angle)
        pts = []
        for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
            u, v = a * self.length / 2, b * self.width / 2
            pts.append((self.x + c * u - s_ * v, self.y + s_ * u + c * v))
        return pts

    @property
    def sheet_factor(self):
        """Z_s / Z_L = width / length."""
        return self.width / self.length

    def impedance(self, omega):
        if self.open_circuit:
            return complex(OPEN_OHM)
        if self.z_of_f is not None:
            return complex(self.z_of_f(omega / (2.0 * np.pi)))
        z = complex(self.R, omega * self.L)
        if self.C is not None:
            z += 1.0 / (1j * omega * self.C)
        return z

    def region(self):
        """`(x0, x1, y0, y1)` of the load rectangle (the bounding box of a rotated one)."""
        if self.is_rotated:
            xs, ys = zip(*self.corners(), strict=True)
            return min(xs), max(xs), min(ys), max(ys)
        if self.direction == "x":
            hx, hy = self.length / 2, self.width / 2
        else:
            hx, hy = self.width / 2, self.length / 2
        return self.x - hx, self.x + hx, self.y - hy, self.y + hy


def _tent_product_integral(c1, c2, h, a, b):
    """Integral over [a, b] of tent(x - c1; h) * tent(x - c2; h), tent(u; h) = max(0, 1 - |u|/h). Exact (piecewise quadratic)."""
    pts = sorted({a, b, *(c + s * h for c in (c1, c2) for s in (-1.0, 0.0, 1.0))})
    pts = [p for p in pts if a <= p <= b]
    total = 0.0
    for lo, hi in zip(pts[:-1], pts[1:], strict=True):
        if hi - lo <= 0.0:
            continue
        mid, half = 0.5 * (lo + hi), 0.5 * (hi - lo)
        xs = mid + half * _GAUSS2[0]
        f = np.maximum(0.0, 1.0 - np.abs(xs - c1) / h) * np.maximum(0.0, 1.0 - np.abs(xs - c2) / h)
        total += half * float(np.sum(f * _GAUSS2[1]))
    return total


def _interval_overlap(c, width, a, b):
    return max(0.0, min(c + width / 2, b) - max(c - width / 2, a))


def check_load_geometry(loads, unit_cell, nx=None, ny=None):
    """Raise if a load does not fit in the cell (with a one-cell margin for the basis-function images); warn if its edges are
    off the mesh grid."""
    px, py = unit_cell.lattice.period_x, unit_cell.lattice.period_y
    for k, ld in enumerate(loads):
        x0, x1, y0, y1 = ld.region()
        mx = px / nx if nx else 0.0
        my = py / ny if ny else 0.0
        if (
            x0 < -px / 2 + mx - 1e-12
            or x1 > px / 2 - mx + 1e-12
            or y0 < -py / 2 + my - 1e-12
            or y1 > py / 2 - my + 1e-12
        ):
            raise ValueError(f"load {k} must lie at least one mesh cell inside the unit cell (basis-function images)")
        # the region must be conductor in the pattern: otherwise no rooftop crosses it and the load has no effect at all
        xs = np.linspace(x0, x1, 9)[1:-1]
        ys = np.linspace(y0, y1, 9)[1:-1]
        gx, gy = np.meshgrid(xs, ys, indexing="ij")
        if ld.is_rotated:
            keep = ld.contains(gx, gy)
            gx, gy = gx[keep], gy[keep]
        if not np.all(unit_cell.is_metal(gx, gy)):
            warnings.warn(
                f"load {k}: part of its region is not metal in the unit-cell pattern, so no basis function crosses it there "
                "and the load has no (or a partial) effect: include the region in the geometry as conductor",
                UserWarning,
                stacklevel=3,
            )
        if nx and ny and not ld.is_rotated:
            for v, p, n in ((x0, px, nx), (x1, px, nx), (y0, py, ny), (y1, py, ny)):
                u = (v + p / 2) / (p / n)
                if abs(u - round(u)) > 1e-6:
                    warnings.warn(
                        f"load {k}: region edge is off the mesh grid ({u:.3f} cells); results depend on the mesh alignment",
                        UserWarning,
                        stacklevel=3,
                    )
                    break


def _add_region(out, basis, direction, zs, region):
    """Add `zs * integral_region f_i . f_j` for the basis functions of `direction` into `out`."""
    x0, x1, y0, y1 = region
    dx, dy = basis.dx, basis.dy
    if direction == "x":
        idx = np.arange(basis.n_x)
        tent_c, tent_h, tent_lo, tent_hi = basis.center_x[idx], dx, x0, x1
        row_c, row_w, row_lo, row_hi = basis.center_y[idx], dy, y0, y1
    else:
        idx = np.arange(basis.n_x, basis.n_unknowns)
        tent_c, tent_h, tent_lo, tent_hi = basis.center_y[idx], dy, y0, y1
        row_c, row_w, row_lo, row_hi = basis.center_x[idx], dx, x0, x1
    touch = (
        (tent_c + tent_h > tent_lo)
        & (tent_c - tent_h < tent_hi)
        & (row_c + row_w / 2 > row_lo)
        & (row_c - row_w / 2 < row_hi)
    )
    sel = np.flatnonzero(touch)
    for a in sel:
        for b in sel:
            if abs(row_c[a] - row_c[b]) > 1e-9 * row_w:
                continue  # different rows: the rect factors do not overlap
            if abs(tent_c[a] - tent_c[b]) >= 2 * tent_h:
                continue
            t = _tent_product_integral(tent_c[a], tent_c[b], tent_h, tent_lo, tent_hi)
            if t == 0.0:
                continue
            r = _interval_overlap(row_c[a], row_w, row_lo, row_hi)
            out[idx[a], idx[b]] += zs * t * r


def rotated_load_points(ld, h):
    """Composite 2-point Gauss rule on the rotated load rectangle with sub-cells of side <= h/4: returns `(px, py, weights)`."""
    ns = max(1, math.ceil(ld.length / (h / 4.0)))
    nt = max(1, math.ceil(ld.width / (h / 4.0)))
    g = _GAUSS2[0]
    ls = (-ld.length / 2 + (np.arange(ns)[:, None] + 0.5 + 0.5 * g[None, :]) * ld.length / ns).ravel()
    ts = (-ld.width / 2 + (np.arange(nt)[:, None] + 0.5 + 0.5 * g[None, :]) * ld.width / nt).ravel()
    sl, st = np.meshgrid(ls, ts, indexing="ij")
    c, s_ = math.cos(ld.angle), math.sin(ld.angle)
    px = ld.x + c * sl - s_ * st
    py = ld.y + s_ * sl + c * st
    w = np.full(sl.shape, (ld.length / (2 * ns)) * (ld.width / (2 * nt)))
    return px.ravel(), py.ravel(), w.ravel()


def _add_rotated(out, basis, ld, zs):
    """Add the sheet-impedance tensor  zs u u^T  (+ OPEN_OHM v v^T when transverse conduction is blocked) integrated over the rotated
    rectangle: Z_ij = sum_k w_k [zs (u.f_i)(u.f_j) + zt (v.f_i)(v.f_j)](p_k). The rooftops are the SAME functions as `_add_region`'s
    (unit-height tent x unit box, x- and y-directed); the quadrature is a fine composite Gauss rule (convergent, not exact: the region's
    edges cut the mesh cells arbitrarily)."""
    dx, dy = basis.dx, basis.dy
    px, py, w = rotated_load_points(ld, min(dx, dy))
    x0, x1, y0, y1 = ld.region()
    nx_ = basis.n_x
    idx = np.arange(basis.n_unknowns)
    cxs, cys = basis.center_x, basis.center_y
    is_x = idx < nx_
    # support (bounding box) of each basis function, to select the ones that can touch the region
    half_x = np.where(is_x, dx, dx / 2)
    half_y = np.where(is_x, dy / 2, dy)
    touch = (cxs + half_x > x0) & (cxs - half_x < x1) & (cys + half_y > y0) & (cys - half_y < y1)
    sel = np.flatnonzero(touch)
    if sel.size == 0:
        return
    fx = np.zeros((sel.size, px.size))
    fy = np.zeros((sel.size, px.size))
    for r, i in enumerate(sel):
        if is_x[i]:
            fx[r] = np.maximum(0.0, 1.0 - np.abs(px - cxs[i]) / dx) * (np.abs(py - cys[i]) < dy / 2)
        else:
            fy[r] = np.maximum(0.0, 1.0 - np.abs(py - cys[i]) / dy) * (np.abs(px - cxs[i]) < dx / 2)
    c, s_ = math.cos(ld.angle), math.sin(ld.angle)
    fu = fx * c + fy * s_
    block = zs * (fu * w) @ fu.T
    if ld.open_circuit or not ld.transverse_conducts:
        fv = -fx * s_ + fy * c
        block = block + complex(OPEN_OHM) * (fv * w) @ fv.T
    out[np.ix_(sel, sel)] += block


def load_matrix(basis, loads, omega):
    """`Z_load` (n_unknowns x n_unknowns, complex) for `loads` on the rooftop `basis` at angular frequency `omega`."""
    n = basis.n_unknowns
    out = np.zeros((n, n), dtype=complex)
    for ld in loads:
        if ld.is_rotated:
            _add_rotated(out, basis, ld, ld.impedance(omega) * ld.sheet_factor)
            continue
        region = ld.region()
        _add_region(out, basis, ld.direction, ld.impedance(omega) * ld.sheet_factor, region)
        if ld.open_circuit or not ld.transverse_conducts:
            _add_region(out, basis, "y" if ld.direction == "x" else "x", complex(OPEN_OHM), region)
    return out


def absorbed_power(basis, loads, omega, currents):
    """Time-average power dissipated in the loads per unit cell: 0.5 Re(I^H Z_load I) [W] for rooftop coefficients `currents` (peak phasors)."""
    zl = load_matrix(basis, loads, omega)
    i = np.asarray(currents, dtype=complex)
    return 0.5 * float(np.real(np.vdot(i, zl @ i)))
