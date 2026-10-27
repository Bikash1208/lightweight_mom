"""Phase 28 (HFSS-first path): exact polygon description of the G cell (`RationalSectorGeometry`, 6.4 mm, 'bowtie' mode) for the AEDT scripts.

Everything is derived from the repo's `RationalSectorCell` (hub, arms, diagonal pieces from the geometry fields; the eight neck rectangles and the four ON-sets of the states from
`neck_rect` / `states()`), in cell-local mm. `check_against_mask()` rasterises the polygons and compares them with `in_feature` of the MoM geometry."""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from phase28_small_board import g_cell  # noqa: E402

from src.ris.cells import ON, madp_kit  # noqa: E402


def _rot90(pts, k):
    out = np.asarray(pts, dtype=float)
    for _ in range(k % 4):
        out = np.stack([-out[:, 1], out[:, 0]], axis=1)
    return out


def cell_polygons_mm():
    """`{piece_name: [(x, y), ...]}` of the metal pieces (hub, 4 arms, 4 diagonal pieces); the necks are the diode sheets (see `necks_mm`)."""
    g = g_cell().geometry()
    r, t, h, a, d = g.r_in, g.gap_t, g.hub_half, g.arm_outer, g.diag_outer
    hub = [(-h, -h), (h, -h), (h, h), (-h, h)]
    arm = [(r, -(r / 2 - t)), (a, -(a / 2 - t)), (a, a / 2 - t), (r, r / 2 - t)]
    diag = [(r, r / 2 + t), (d, d / 2 + t), (d, d), (d / 2 + t, d), (r / 2 + t, r), (r, r)]
    out = {"hub": np.asarray(hub) * 1e3}
    for k in range(4):
        out[f"arm{k}"] = _rot90(np.asarray(arm) * 1e3, k)
    for k in range(4):
        out[f"diag{k}"] = _rot90(np.asarray(diag) * 1e3, k)
    return {n: [(float(x), float(y)) for x, y in p] for n, p in out.items()}


def necks_mm():
    """The eight diode sites in the order of `RationalSectorCell._sites()`: `(rect corners (x0,x1,y0,y1) mm, current line ((xa,ya),(xb,yb)) mm, size length x width mm)`."""
    cell = g_cell()
    geo = cell.geometry()
    out = []
    for kind, k, xc, yc, length, width, axis in cell._sites():
        x0, x1, y0, y1 = (v * 1e3 for v in geo.neck_rect("arm_neck" if kind == "arm_neck" else "diag", k))
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        line = ((x0, cy), (x1, cy)) if axis == "x" else ((cx, y0), (cx, y1))
        out.append(
            {
                "site": len(out),
                "kind": kind,
                "k": k,
                "rect_mm": (x0, x1, y0, y1),
                "current_line_mm": line,
                "length_mm": float(length * 1e3),
                "width_mm": float(width * 1e3),
                "axis": axis,
            }
        )
    return out


def state_on_sites():
    """`[[sites ON in state 0], ...]` (site index into `necks_mm()`)."""
    cell = g_cell()
    return [[i for i, s in enumerate(st) if s == ON] for st in cell.states()]


def diode_values():
    kit = madp_kit()
    return {"r_on_ohm": float(kit.on.R_s), "c_off_pf": float(kit.off.C_T * 1e12), "l_s_nh": float(kit.on.L_s * 1e9)}


def _inside(poly, x, y):
    p = np.asarray(poly)
    ok = np.zeros(x.shape, dtype=bool)
    n = len(p)
    for i in range(n):
        x1, y1 = p[i]
        x2, y2 = p[(i + 1) % n]
        cond = (y1 > y) != (y2 > y)
        with np.errstate(divide="ignore", invalid="ignore"):
            xi = (x2 - x1) * (y - y1) / (y2 - y1) + x1
        ok ^= cond & (x < xi)
    return ok


def check_against_mask(n=801):
    """Fraction of grid points (0.0075 mm step, away from edges) where polygons + necks and the MoM `in_feature` mask disagree. Returns `(n_disagree, n_total, n_near_edge)`."""
    g = g_cell().geometry()
    c = g.clip * 1e3
    xs = np.linspace(-c, c, n)
    xg, yg = np.meshgrid(xs, xs)
    mask = g.in_feature(xg * 1e-3, yg * 1e-3)
    poly = np.zeros(xg.shape, dtype=bool)
    for p in cell_polygons_mm().values():
        poly |= _inside(p, xg, yg)
    for nk in necks_mm():
        x0, x1, y0, y1 = nk["rect_mm"]
        poly |= (xg >= x0) & (xg <= x1) & (yg >= y0) & (yg <= y1)
    bad = mask != poly
    # points within one grid step of an edge may legitimately differ
    step = xs[1] - xs[0]
    near = np.zeros(xg.shape, dtype=bool)
    for p in list(cell_polygons_mm().values()):
        pa = np.asarray(p)
        for i in range(len(pa)):
            a, b = pa[i], pa[(i + 1) % len(pa)]
            t = np.clip(((xg - a[0]) * (b[0] - a[0]) + (yg - a[1]) * (b[1] - a[1])) / max(np.sum((b - a) ** 2), 1e-18), 0, 1)
            dist = np.hypot(xg - (a[0] + t * (b[0] - a[0])), yg - (a[1] + t * (b[1] - a[1])))
            near |= dist < 1.5 * step
    for nk in necks_mm():
        x0, x1, y0, y1 = nk["rect_mm"]
        near |= (np.abs(xg - x0) < 1.5 * step) | (np.abs(xg - x1) < 1.5 * step) | (np.abs(yg - y0) < 1.5 * step) | (np.abs(yg - y1) < 1.5 * step)
    return int(np.sum(bad & ~near)), int(bad.size), int(np.sum(bad & near))
