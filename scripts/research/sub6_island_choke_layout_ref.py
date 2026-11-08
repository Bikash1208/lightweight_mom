"""SUB6 island cell: bias-choke LAYOUT (meander + radial stub) on the bottom copper, fit check under one switched cell and against its periodic neighbours, and the FR4-thickness table.
MODEL_ONLY.  Frame of the island package: slot along x, island via at (0, -0.20), disc via at (0, +2.35) (pad r 0.30, drill 0.30).

    python sub6_island_choke_layout.py table      -> results/sub6_island_choke/{choke_table.json, CHOKE.md}
"""
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sub6_island_choke as ch

ROOT = ch.ROOT
RES = ROOT / "results" / "sub6_island_choke"
P = 11.6
PAD_R = 0.30
VIA_I, VIA_D = (0.0, -0.20), (0.0, 2.35)
LEAD, LINK = 0.50, 0.50  # mm: straight lead from the via centre; straight link between the last meander run and the stub apex
LS_MIN, LS_MAX = 0.8, 4.0
DIRS = {"+x": (1, 0), "-x": (-1, 0), "+y": (0, 1), "-y": (0, -1)}


def rect(c0, c1, w):
    """axis-aligned rectangle around the segment c0-c1 (square ends extended by w/2)."""
    (x0, y0), (x1, y1) = c0, c1
    h = w / 2
    return [(min(x0, x1) - h, min(y0, y1) - h), (max(x0, x1) + h, min(y0, y1) - h), (max(x0, x1) + h, max(y0, y1) + h), (min(x0, x1) - h, max(y0, y1) + h)]


def circle_poly(c, r, n=24):
    return [(c[0] + r * math.cos(2 * math.pi * k / n), c[1] + r * math.sin(2 * math.pi * k / n)) for k in range(n)]


def sector_poly(apex, d, r2, ang_deg, n=16):
    a0 = math.atan2(d[1], d[0])
    half = math.radians(ang_deg) / 2
    pts = [apex] + [(apex[0] + r2 * math.cos(a0 - half + 2 * half * k / n), apex[1] + r2 * math.sin(a0 - half + 2 * half * k / n)) for k in range(n + 1)]
    return pts


def build(via, dname, flip, w, s, l_mm, r2, ang, n_runs):
    """-> dict(path=[centreline points], parts=[polygons in order], apex, runs, ls) or None when Ls is out of range."""
    d = np.array(DIRS[dname], float)
    nrm = np.array([-d[1], d[0]]) * flip
    pm = w + s
    ls = (l_mm - LEAD - LINK - (n_runs - 1) * pm) / n_runs
    if not (LS_MIN <= ls <= LS_MAX):
        return None
    p = np.array(via, float)
    path = [tuple(p)]
    p = p + d * LEAD
    path.append(tuple(p))
    parts = [circle_poly(via, PAD_R), rect(path[0], path[1], w)]
    side = 1.0
    for k in range(n_runs):
        q = p + nrm * side * ls
        parts.append(rect(tuple(p), tuple(q), w))
        path.append(tuple(q))
        p = q
        if k < n_runs - 1:
            q = p + d * pm
            parts.append(rect(tuple(p), tuple(q), w))
            path.append(tuple(q))
            p = q
        side = -side
    q = p + d * LINK
    parts.append(rect(tuple(p), tuple(q), w))
    path.append(tuple(q))
    parts.append(sector_poly(tuple(q), tuple(d), r2, ang))
    return dict(via=via, dir=dname, flip=flip, path=path, parts=parts, apex=tuple(q), n_runs=n_runs, ls=ls, pitch=pm, w=w, r2=r2, angle=ang, spacing=s)


def _seg_dist(p, a, b):
    ab, ap = np.array(b) - np.array(a), np.array(p) - np.array(a)
    t = np.clip(ap @ ab / (ab @ ab + 1e-30), 0, 1)
    return float(np.linalg.norm(ap - t * ab))


def _segs_cross(a, b, c, d):
    def o(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    return o(a, b, c) * o(a, b, d) < 0 and o(c, d, a) * o(c, d, b) < 0


def _inside(pt, poly):
    x, y = pt
    c = False
    for i in range(len(poly)):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % len(poly)]
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            c = not c
    return c


def _bbox_gap(A, B):
    ax, ay = zip(*A)
    bx, by = zip(*B)
    dx = max(min(bx) - max(ax), min(ax) - max(bx), 0.0)
    dy = max(min(by) - max(ay), min(ay) - max(by), 0.0)
    return math.hypot(dx, dy)


def poly_dist(A, B):
    g = _bbox_gap(A, B)
    if g > 0.5:
        return g  # far apart: the bounding-box gap is a lower bound, good enough for the >= 0.15 mm test
    nA, nB = len(A), len(B)
    for i in range(nA):
        for j in range(nB):
            if _segs_cross(A[i], A[(i + 1) % nA], B[j], B[(j + 1) % nB]):
                return 0.0
    if _inside(A[0], B) or _inside(B[0], A):
        return 0.0
    d = min(_seg_dist(a, B[j], B[(j + 1) % nB]) for a in A for j in range(nB))
    return min(d, min(_seg_dist(b, A[i], A[(i + 1) % nA]) for b in B for i in range(nA)))


def shift(parts, dx, dy):
    return [[(x + dx, y + dy) for x, y in poly] for poly in parts]


def min_gap(chokes, images=True, self_check=True):
    """smallest gap (mm) between the metal of different chokes (and their periodic images), plus self-spacing between non-adjacent parts of one choke."""
    gaps = {"between_chokes_in_cell": 1e9, "to_periodic_images": 1e9, "self": 1e9}
    for (i, a), (j, b) in itertools.combinations(enumerate(chokes), 2):
        gaps["between_chokes_in_cell"] = min(gaps["between_chokes_in_cell"], min(poly_dist(pa, pb) for pa in a["parts"] for pb in b["parts"]))
    if images:
        for a in chokes:
            for b in chokes:
                for ix, iy in itertools.product((-1, 0, 1), repeat=2):
                    if (ix, iy) == (0, 0):
                        continue
                    gaps["to_periodic_images"] = min(gaps["to_periodic_images"], min(poly_dist(pa, pb) for pa in a["parts"] for pb in shift(b["parts"], ix * P, iy * P)))
    if self_check:
        for a in chokes:
            pr = a["parts"]
            for i in range(len(pr)):
                for j in range(i + 2, len(pr)):
                    if i == 0:
                        continue  # the pad touches the lead/first run
                    gaps["self"] = min(gaps["self"], poly_dist(pr[i], pr[j]))
    return gaps


def spacing_rule(h):
    return max(ch.G_MIN, round(3.0 * h, 3))  # meander gap >= 3 h (coupling below about -20 dB), never below the 0.15 mm limit


def candidates(des, h):
    """all layout combinations for one design -> list of (gap dict, chokes) that satisfy trace/gap >= 0.15 mm in the cell and with the periodic neighbours."""
    s = spacing_rule(h)
    single = {}
    for name, via in (("I", VIA_I), ("D", VIA_D)):
        single[name] = []
        for dn, fl, n in itertools.product(DIRS, (1, -1), range(1, 9)):
            c = build(via, dn, fl, des["w_mm"], s, des["l_mm"], des["r2_mm"], des["angle_deg"], n)
            if c is None:
                continue
            g = min_gap([c])
            if min(g["to_periodic_images"], g["self"]) >= ch.G_MIN - 1e-9:
                single[name].append(c)
    out = []
    for ca in single["I"]:
        for cb in single["D"]:
            g = min_gap([ca, cb])
            if min(g.values()) >= ch.G_MIN - 1e-9:
                out.append((g, [ca, cb]))
    return out


def isolated(des, h):
    """fits ignoring the periodic neighbours (one switched cell alone)."""
    s = spacing_rule(h)
    out = []
    for dA, dB, fA, fB in itertools.product(DIRS, DIRS, (1, -1), (1, -1)):
        for nA in range(1, 9):
            ca = build(VIA_I, dA, fA, des["w_mm"], s, des["l_mm"], des["r2_mm"], des["angle_deg"], nA)
            if ca is None or min_gap([ca], images=False)["self"] < ch.G_MIN:
                continue
            for nB in range(1, 9):
                cb = build(VIA_D, dB, fB, des["w_mm"], s, des["l_mm"], des["r2_mm"], des["angle_deg"], nB)
                if cb is None or min_gap([cb], images=False)["self"] < ch.G_MIN:
                    continue
                g = min_gap([ca, cb], images=False, self_check=False)
                if g["between_chokes_in_cell"] >= ch.G_MIN:
                    out.append((g, [ca, cb]))
                    return out
    return out


def best_layout(des, h):
    ok = candidates(des, h)
    ok_cell = ok if ok else isolated(des, h)
    def key(c):
        return (-(min(c[0]["to_periodic_images"], 2.0)), sum(x["n_runs"] for x in c[1]))
    pick = sorted(ok, key=key)[0] if ok else None
    pick_cell = sorted(ok_cell, key=lambda c: sum(x["n_runs"] for x in c[1]))[0] if ok_cell else None
    return pick, pick_cell


def table():
    RES.mkdir(parents=True, exist_ok=True)
    rows = []
    for h in ch.H_LIST:
        best = None
        for ang in (90, 120):
            des = ch.design(h, ang)
            if des is None:
                continue
            pick, pick_cell = best_layout(des, h)
            rec = dict(des, spacing_mm=spacing_rule(h), fits_all_neighbours=bool(pick), fits_isolated_cell=bool(pick_cell),
                       layout=None if pick is None else dict(gaps=pick[0], chokes=[dict(via=c["via"], dir=c["dir"], flip=c["flip"], n_runs=c["n_runs"], ls_mm=c["ls"], apex=c["apex"]) for c in pick[1]]))
            rec["_chokes"] = None if pick is None else pick[1]
            rows.append(rec)
            print(h, ang, "Z0 %.0f |Zin|min %.0f Yload %.0f uS fits(all)=%s fits(isolated)=%s" % (rec["ms"]["z0"], rec["zin_min_ohm"], rec["Yload_max_uS"], rec["fits_all_neighbours"], rec["fits_isolated_cell"]), flush=True)
    for r in rows:
        r.pop("_chokes")
    (RES / "choke_table.json").write_text(json.dumps(rows, indent=1, default=lambda o: [o.real, o.imag] if isinstance(o, complex) else float(o)))
    return rows


if __name__ == "__main__":
    {"table": table}[sys.argv[1]]()
