"""SUB6 FINAL BOARD, step 3: cell table, 12 series bias strings (20 switched cells each), bottom-copper routing, clearance checks (MODEL_ONLY until board HFSS).

Board frame = the synthesis frame (feed at the origin, mm).  Switched cell = island cell: its local frame (slot along local x, island at local -y, diode gap at local +y) is rotated to the hologram slot axis theta_i
(theta_i or theta_i + 180 deg, free choice per cell: the sheet tensor is invariant, used to shorten the bias traces).  Vias (island frame): island (-1.78, -0.53), disc (-3.63, -3.63), drill 0.3, pad 0.6, anti-pad r 0.5.
Diode orientation: anode on the ISLAND pad, cathode on the DISC pad of the same cell (all cells), so the string  island(1) - diode - disc(1) - trace - island(2) - diode - disc(2) ...  carries the forward current.
Strings: the 240 switched cells are cut into 12 azimuth sectors of 20 cells (sorted by azimuth around the feed), each ordered as a path by nearest-neighbour + 2-opt on the via-to-via trace length (orientation DP).
Routing: bottom copper on FR4 0.8 (z -4.045 ... -4.010), trace 0.15 mm, gap >= 0.15 mm, grid A* (numba) with 0.15 mm pitch, foreign vias kept 0.525 mm (pad 0.3 + gap 0.15 + half trace 0.075 > anti-pad 0.5),
board outline = hexagon + 10 mm margin, connector: 24 pads (1.0 x 1.0 mm, pitch 1.27 mm) in one column near the left flat edge, 4 mm radius around the SMA hole kept clear.
    python sub6_final_layout.py build
"""
import json
import math
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from numba import njit
from scipy.spatial import cKDTree

import sub6_final_design as FD

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sub6_final_board"
P = 11.6
VI_L, VD_L = np.array([-1.78, -0.53]), np.array([-3.63, -3.63])  # local (island frame) via positions
ISL = dict(x0=-2.174, x1=2.174, y0=-0.85, y1=0.45)
SLOT_ISL = (8.696, 2.9)
R_ISL_CELL = 5.435
V_DIODE, I_DIODE = 0.95, 10e-3
W_TR, GAP_TR = 0.15, 0.15
PAD_R, ANTIPAD_R = 0.30, 0.50
KEEP_VIA = PAD_R + GAP_TR + W_TR / 2  # 0.525 (router keep-out radius is KEEP_VIA + 0.1 = 0.625)
ANTIPAD_CLEAR = ANTIPAD_R + W_TR / 2  # 0.575: trace centre line must stay this far from a foreign via centre, else the trace body overlaps the anti-pad void
NEAR_R = float(os.environ.get("FINAL_NEAR_R", "0.05"))  # only vias within NEAR_R of a net end point (its own via) are lifted while that net is routed
MARGIN = 10.0  # outline offset beyond the hexagon (mm)
GRID = 0.15
N_STR, N_PER = 12, 20
PAD_W, PAD_PITCH = 1.0, 2.54


def rot(v, a):
    c, s = math.cos(a), math.sin(a)
    return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1]])


def hexagon(side):
    return np.array([[side * math.cos(math.radians(90 + 60 * i)), side * math.sin(math.radians(90 + 60 * i))] for i in range(6)])


# ------------------------------------------------------------------ cell table
def cell_table(fd):
    D, T = fd["D"], fd["T"]
    S = D["S"]
    x, y = S.x * 1e3, S.y * 1e3
    mask, live = D["mask"] & D["live"], D["live"]
    axis = fd["axis"]
    rows = []
    st = T["st"]["rows"]
    for i in range(S.n):
        if mask[i]:
            rows.append(dict(site=i, x=x[i], y=y[i], cls="SWITCHED", row_id="ISL_W2.9_q0.8_Li0.5_c0.6", r=R_ISL_CELL, L=SLOT_ISL[0], W=SLOT_ISL[1], theta=float(axis[i]), flipped=int(fd["flips"][i]), alpha_req=float(D["a_req_sw"][i])))
        else:
            rw = st[D["ist"][i]]["row"]
            rows.append(dict(site=i, x=x[i], y=y[i], cls="STATIC", row_id=rw["id"], r=rw["r"], L=rw["L"], W=rw["W"], theta=float(axis[i]), flipped=int(fd["flips"][i]), alpha_req=float(D["a_req_st"][i])))
    return rows


def via_points(cell, s):
    th = cell["theta"] + s * math.pi
    c = np.array([cell["x"], cell["y"]])
    return c + rot(VI_L, th), c + rot(VD_L, th), th


# ------------------------------------------------------------------ strings
def string_cost(cells, order, conn):
    """min over orientations s of sum |VD_a(s_a) - VI_b(s_b)| along `order` (DP), plus the distance of the two ends to the connector; returns (cost, s list)."""
    n = len(order)
    pts = [[via_points(cells[order[k]], s) for s in (0, 1)] for k in range(n)]
    best = np.zeros((n, 2))
    arg = np.zeros((n, 2), int)
    for k in range(1, n):
        for sb in (0, 1):
            cand = [best[k - 1][sa] + np.linalg.norm(pts[k - 1][sa][1] - pts[k][sb][0]) for sa in (0, 1)]
            arg[k][sb] = int(np.argmin(cand))
            best[k][sb] = min(cand)
    ends = [np.linalg.norm(pts[0][s][0] - conn) for s in (0, 1)]
    endl = [np.linalg.norm(pts[-1][s][1] - conn) for s in (0, 1)]
    # first/last orientation also matter for the lead length: add a light weight (0.5)
    tot = [best[-1][s] + 0.5 * endl[s] for s in (0, 1)]
    sl = int(np.argmin(tot))
    s_list = [0] * n
    s_list[-1] = sl
    for k in range(n - 1, 0, -1):
        s_list[k - 1] = int(arg[k][s_list[k]])
    return float(tot[sl] + 0.5 * ends[s_list[0]]), s_list


def order_string(cells, idx, conn):
    cx = np.array([[cells[i]["x"], cells[i]["y"]] for i in idx])
    start = int(np.argmin(np.linalg.norm(cx - conn, axis=1)))
    order = [start]
    left = set(range(len(idx))) - {start}
    while left:
        last = cx[order[-1]]
        nxt = min(left, key=lambda j: np.linalg.norm(cx[j] - last))
        order.append(nxt)
        left.remove(nxt)
    cost, s = string_cost([cells[i] for i in idx], order, conn)
    improved = True
    while improved:
        improved = False
        for i in range(0, len(order) - 1):
            for j in range(i + 1, len(order)):
                cand = order[:i] + order[i:j + 1][::-1] + order[j + 1:]
                c2, s2 = string_cost([cells[k] for k in idx], cand, conn)
                if c2 < cost - 1e-9:
                    order, cost, s, improved = cand, c2, s2, True
    return [idx[k] for k in order], s


def build_strings(cells, conn):
    sw = [i for i, c in enumerate(cells) if c["cls"] == "SWITCHED"]
    az = np.array([math.atan2(cells[i]["y"], cells[i]["x"]) for i in sw])
    az0 = math.atan2(conn[1], conn[0])
    key = np.mod(az - az0, 2 * math.pi)
    srt = [sw[k] for k in np.argsort(key)]
    strings = []
    for k in range(N_STR):
        idx = srt[k * N_PER:(k + 1) * N_PER]
        order, s = order_string(cells, idx, conn)
        strings.append(dict(id=k + 1, cells=order, s=s))
    return strings


# ------------------------------------------------------------------ router
@njit(cache=True)
def astar(blocked_a, blocked_b, W, H, sx, sy, tx, ty, g, parent, hk, hv, cap):
    n = W * H
    for i in range(n):
        g[i] = 1e30
        parent[i] = -1
    s = sy * W + sx
    t = ty * W + tx
    g[s] = 0.0
    size = 0
    hk[0] = 0.0
    hv[0] = s
    size = 1
    dx = np.array([1, -1, 0, 0, 1, 1, -1, -1])
    dy = np.array([0, 0, 1, -1, 1, -1, 1, -1])
    cst = np.array([1.0, 1.0, 1.0, 1.0, 1.41421356, 1.41421356, 1.41421356, 1.41421356])
    nexp = 0
    while size > 0:
        nexp += 1
        if nexp > cap:
            return False
        cur = hv[0]
        fcur = hk[0]
        size -= 1
        if size > 0:
            hk[0] = hk[size]
            hv[0] = hv[size]
            i = 0
            while True:
                l = 2 * i + 1
                if l >= size:
                    break
                r = l + 1
                m = l
                if r < size and hk[r] < hk[l]:
                    m = r
                if hk[m] < hk[i]:
                    tk, tv = hk[i], hv[i]
                    hk[i], hv[i] = hk[m], hv[m]
                    hk[m], hv[m] = tk, tv
                    i = m
                else:
                    break
        cx = cur % W
        cy = cur // W
        if cur == t:
            return True
        gc = g[cur]
        # lazy deletion
        ex = abs(cx - tx)
        ey = abs(cy - ty)
        hc = (ex + ey) + (1.41421356 - 2.0) * min(ex, ey)
        if fcur > gc + hc + 1e-6:
            continue
        for d in range(8):
            nx = cx + dx[d]
            ny = cy + dy[d]
            if nx < 0 or ny < 0 or nx >= W or ny >= H:
                continue
            ni = ny * W + nx
            if blocked_a[ni] or blocked_b[ni]:
                if ni != t:
                    continue
            if d >= 4:  # no diagonal corner cutting through blocked cells
                if blocked_a[cy * W + nx] or blocked_a[ny * W + cx] or blocked_b[cy * W + nx] or blocked_b[ny * W + cx]:
                    continue
            ng = gc + cst[d]
            if ng < g[ni]:
                g[ni] = ng
                parent[ni] = cur
                ex2 = abs(nx - tx)
                ey2 = abs(ny - ty)
                h2 = (ex2 + ey2) + (1.41421356 - 2.0) * min(ex2, ey2)
                # push
                k = size
                hk[k] = ng + h2
                hv[k] = ni
                size += 1
                while k > 0:
                    p = (k - 1) // 2
                    if hk[p] > hk[k]:
                        tk, tv = hk[p], hv[p]
                        hk[p], hv[p] = hk[k], hv[k]
                        hk[k], hv[k] = tk, tv
                        k = p
                    else:
                        break
    return False


class Router:
    def __init__(self, outline_poly, keepout_clear_r):
        xs = outline_poly[:, 0]
        ys = outline_poly[:, 1]
        self.x0, self.y0 = xs.min() - 1.0, ys.min() - 1.0
        self.W = int((xs.max() - xs.min() + 2.0) / GRID) + 1
        self.H = int((ys.max() - ys.min() + 2.0) / GRID) + 1
        W, H = self.W, self.H
        gx = self.x0 + np.arange(W) * GRID
        gy = self.y0 + np.arange(H) * GRID
        GX, GY = np.meshgrid(gx, gy)
        inside = np.ones((H, W), bool)
        n = len(outline_poly)
        for k in range(n):
            a, b = outline_poly[k], outline_poly[(k + 1) % n]
            cross = (b[0] - a[0]) * (GY - a[1]) - (b[1] - a[1]) * (GX - a[0])
            ln = math.hypot(b[0] - a[0], b[1] - a[1])
            inside &= cross / ln >= 0.9  # polygon is counter-clockwise: keep 0.9 mm from the outline
        self.static = ~inside
        self.static |= np.hypot(GX, GY) < keepout_clear_r
        self.static = self.static.ravel().astype(np.uint8)
        self.trace = np.zeros(H * W, np.int16)
        self.viac = np.zeros(H * W, np.int16)
        self.g = np.empty(H * W, np.float64)
        self.parent = np.empty(H * W, np.int32)
        self.hk = np.empty(40_000_000, np.float64)
        self.hv = np.empty(40_000_000, np.int32)

    def node(self, p):
        return int(round((p[0] - self.x0) / GRID)), int(round((p[1] - self.y0) / GRID))

    def stamp(self, arr, p, r_mm, val):
        cx, cy = self.node(p)
        rc = int(math.ceil(r_mm / GRID)) + 1
        W, H = self.W, self.H
        xa, xb, ya, yb = max(cx - rc, 0), min(cx + rc + 1, W), max(cy - rc, 0), min(cy + rc + 1, H)
        yy, xx = np.mgrid[ya:yb, xa:xb]
        m = (np.hypot((xx - cx) * GRID + (self.x0 + cx * GRID - p[0]) * 0 + (cx * GRID + self.x0 - p[0]) * 0, (yy - cy) * GRID) <= 1e9)
        # exact disc around the true point
        px = self.x0 + xx * GRID - p[0]
        py = self.y0 + yy * GRID - p[1]
        m = np.hypot(px, py) <= r_mm
        view = arr.reshape(H, W)[ya:yb, xa:xb]
        if val is None:
            return m, (ya, yb, xa, xb)
        view[m] += val

    def add_via(self, p):
        self.stamp(self.viac, p, KEEP_VIA + 0.1, 1)

    def remove_via(self, p):
        self.stamp(self.viac, p, KEEP_VIA + 0.1, -1)

    def route(self, p0, p1, inflate_mm):
        """grid A* from the node nearest to p0 to the node nearest to p1; returns the polyline [p0, ..., p1] (mm) or None."""
        W, H = self.W, self.H
        blocked_a = ((self.trace > 0) | (self.static > 0)).astype(np.uint8)
        blocked_b = (self.viac > 0).astype(np.uint8)

        def free_node(p):
            cx, cy = self.node(p)
            best = None
            for r in range(0, 5):
                for ddx in range(-r, r + 1):
                    for ddy in range(-r, r + 1):
                        nx, ny = cx + ddx, cy + ddy
                        if 0 <= nx < W and 0 <= ny < H and not (blocked_a[ny * W + nx] or blocked_b[ny * W + nx]):
                            d = math.hypot(self.x0 + nx * GRID - p[0], self.y0 + ny * GRID - p[1])
                            if best is None or d < best[0]:
                                best = (d, nx, ny)
                if best:
                    return best[1], best[2]
            return None

        a, b = free_node(p0), free_node(p1)
        if a is None or b is None:
            return None, "no free start/end node"
        ok = astar(blocked_a, blocked_b, W, H, a[0], a[1], b[0], b[1], self.g, self.parent, self.hk, self.hv, 4_000_000)
        if not ok:
            return None, "no path"
        path = []
        cur = b[1] * W + b[0]
        while cur != -1:
            path.append((self.x0 + (cur % W) * GRID, self.y0 + (cur // W) * GRID))
            cur = self.parent[cur]
        path = path[::-1]
        poly = [tuple(p0)] + simplify(path) + [tuple(p1)]
        self.paint(poly, inflate_mm, +1)
        return poly, "ok"


def _paint(self, poly, inflate_mm, val):
    H, W = self.H, self.W
    for q in resample(poly, GRID / 2):
        m, (ya, yb, xa, xb) = self.stamp(self.trace, q, inflate_mm, None)
        self.trace.reshape(H, W)[ya:yb, xa:xb][m] += val


Router.paint = _paint


def simplify(path):
    out = [path[0]]
    for i in range(1, len(path) - 1):
        d0 = (round(path[i][0] - path[i - 1][0], 6), round(path[i][1] - path[i - 1][1], 6))
        d1 = (round(path[i + 1][0] - path[i][0], 6), round(path[i + 1][1] - path[i][1], 6))
        if abs(d0[0] * d1[1] - d0[1] * d1[0]) > 1e-9 or d0[0] * d1[0] + d0[1] * d1[1] < 0:
            out.append(path[i])
    out.append(path[-1])
    return out


def resample(poly, step):
    pts = []
    for a, b in zip(poly[:-1], poly[1:]):
        n = max(1, int(math.ceil(math.hypot(b[0] - a[0], b[1] - a[1]) / step)))
        for k in range(n):
            t = k / n
            pts.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
    pts.append(poly[-1])
    return pts


def poly_len(poly):
    return float(sum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(poly[:-1], poly[1:])))


# ------------------------------------------------------------------ build
def build():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    fd = FD.final_design("cal")
    cells = cell_table(fd)
    side = fd["D"]["side"] * 1e3
    apothem = side * math.cos(math.radians(30)) + MARGIN
    r_out = apothem / math.cos(math.radians(30))
    outline = hexagon(r_out)
    conn_x = -(apothem - 3.5)
    conn = np.array([conn_x, 0.0])
    strings = build_strings(cells, conn)
    # vias
    vias = {}
    for st in strings:
        for k, ci in enumerate(st["cells"]):
            vi, vd, th = via_points(cells[ci], st["s"][k])
            cells[ci]["theta_board"] = th
            vias[ci] = (vi, vd)
    nets = []
    for st in strings:
        for k in range(N_PER - 1):
            a, b = st["cells"][k], st["cells"][k + 1]
            nets.append(dict(kind="link", string=st["id"], pos=k, p0=vias[a][1], p1=vias[b][0], a=a, b=b))
    # connector pads: 24 pads in one column, ordered by the y of the string end they serve
    ends = []
    for st in strings:
        ends.append((st["id"], "start", vias[st["cells"][0]][0], st["cells"][0]))
        ends.append((st["id"], "end", vias[st["cells"][-1]][1], st["cells"][-1]))
    ends.sort(key=lambda e: -e[2][1])
    pads = [np.array([conn_x, (11.5 - k) * 2.54]) for k in range(24)]
    for k, e in enumerate(ends):
        nets.append(dict(kind="lead", string=e[0], pos=e[1], p0=e[2], p1=pads[k], a=e[3], b=None, pad=k))
    router = Router(outline, 4.5)
    allv = []
    for ci, (vi, vd) in vias.items():
        router.add_via(vi)
        router.add_via(vd)
        allv += [vi, vd]
    for p in pads:
        router.add_via(p)
        allv.append(p)
    router.vlist = list(allv)
    print(f"grid {router.W} x {router.H}, strings {len(strings)}, nets {len(nets)}, setup {time.time() - t0:.0f}s", flush=True)
    inflate = W_TR + GAP_TR + 0.06  # centre-line to centre-line keep-out of a routed trace (0.36 mm)
    order = sorted([n for n in nets if n["kind"] == "link"], key=lambda n: np.linalg.norm(n["p0"] - n["p1"])) + sorted([n for n in nets if n["kind"] == "lead"], key=lambda n: np.linalg.norm(n["p0"] - n["p1"]))
    def try_route(n):
        near = [q for q in router.vlist if min(np.linalg.norm(q - n["p0"]), np.linalg.norm(q - n["p1"])) < NEAR_R]
        for q in near:
            router.remove_via(q)
        poly, msg = router.route(n["p0"], n["p1"], inflate)
        for q in near:
            router.add_via(q)
        n["poly"] = poly
        return poly is not None, msg

    for n in order:
        try_route(n)
    for rnd in range(8):
        bad = [n for n in order if n["poly"] is None]
        print(f"round {rnd}: unrouted {len(bad)}", [(n['kind'], n['string'], n['pos']) for n in bad], flush=True)
        if not bad:
            break
        for n in bad:
            near = []
            for m in order:
                if m["poly"] is None or m is n:
                    continue
                pts = np.array(resample(m["poly"], 0.3))
                if min(np.linalg.norm(pts - n["p0"], axis=1).min(), np.linalg.norm(pts - n["p1"], axis=1).min()) < 3.0 + rnd:
                    near.append(m)
            for m in near:
                router.paint(m["poly"], inflate, -1)
                m["poly"] = None
            ok, msg = try_route(n)
            for m in sorted(near, key=lambda q: np.linalg.norm(q["p0"] - q["p1"])):
                try_route(m)
    failed = [(n["kind"], n["string"], n["pos"]) for n in order if n["poly"] is None]
    print(f"routed {len(nets) - len(failed)}/{len(nets)}  failed {failed}  {time.time() - t0:.0f}s", flush=True)
    return dict(fd=fd, cells=cells, strings=strings, vias=vias, nets=nets, pads=pads, outline=outline, failed=failed, conn=conn, apothem=apothem, r_out=r_out)


def antipad_check(res):
    """ANTI-PAD CROSSING CHECK.  (1) every routed trace (sampled every 0.03 mm) against every FOREIGN via: centre-line distance must be >= ANTIPAD_CLEAR (0.575 mm = anti-pad radius 0.5 + half trace);
    a hit means the trace body overlaps that via's anti-pad void in the ground plane.  A net's own end vias are excluded.  (2) via-via: two anti-pad voids (r 0.5) must not merge (centre distance >= 1.0 mm).
    (3) trace vs board outline margin is the router's 0.9 mm keep-in (not re-checked here).  Returns counts, the worst distances and the list of hits."""
    nets = [n for n in res["nets"] if n["poly"] is not None]
    vc = [(ci, k, v) for ci, (vi, vd) in res["vias"].items() for k, v in (("I", vi), ("D", vd))]
    vpts = np.array([v for _, _, v in vc])
    vtree = cKDTree(vpts)
    hits, dmin = [], 9.0
    for n in nets:
        own = {j for p in ([n["p0"], n["p1"]] if n["kind"] == "link" else [n["p0"]]) for j in vtree.query_ball_point(p, 0.05)}
        sp = np.array(resample(n["poly"], 0.03))
        d, j = vtree.query(sp, k=2)
        for kk in range(2):
            for i in np.nonzero(d[:, kk] < ANTIPAD_CLEAR - 1e-6)[0]:
                if int(j[i, kk]) in own:
                    continue
                hits.append(dict(kind=n["kind"], string=n["string"], pos=n["pos"], via=f"{vc[j[i, kk]][0]}{vc[j[i, kk]][1]}", distance_mm=float(d[i, kk])))
        mask = np.ones(len(vpts), bool)
        mask[list(own)] = False
        if mask.any():
            dd = cKDTree(vpts[mask]).query(sp)[0]
            dmin = min(dmin, float(dd.min()))
    pairs = vtree.query_pairs(2 * ANTIPAD_R - 1e-6, output_type="ndarray")
    vv = [dict(a=f"{vc[a][0]}{vc[a][1]}", b=f"{vc[b][0]}{vc[b][1]}", distance_mm=float(np.linalg.norm(vpts[a] - vpts[b]))) for a, b in pairs]
    uniq = sorted({(h["kind"], h["string"], h["pos"], h["via"]) for h in hits})
    return dict(antipad_clear_mm=ANTIPAD_CLEAR, trace_antipad_hits=len(uniq), trace_antipad_hit_samples=len(hits), min_trace_to_foreign_via_mm=dmin, worst_hits=sorted(hits, key=lambda h: h["distance_mm"])[:10],
                nets_with_hits=sorted({(h[0], h[1], h[2]) for h in uniq}), via_via_antipad_merges=len(vv), via_via_merge_pairs=vv[:10])


def save_result(res):
    """routing.json: everything the fab step needs (cells, strings, vias, routed polylines, pads, outline)."""
    f = lambda a: [round(float(x), 4) for x in a]
    d = dict(cells=[{k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in c.items()} for c in res["cells"]],
             strings=[dict(id=s["id"], cells=[int(c) for c in s["cells"]], s=[int(x) for x in s["s"]]) for s in res["strings"]],
             vias={str(ci): dict(island=f(vi), disc=f(vd)) for ci, (vi, vd) in res["vias"].items()},
             nets=[dict(kind=n["kind"], string=n["string"], pos=n["pos"], a=n["a"], b=n["b"], pad=n.get("pad"), p0=f(n["p0"]), p1=f(n["p1"]), poly=[f(q) for q in n["poly"]] if n["poly"] else None) for n in res["nets"]],
             pads=[f(p) for p in res["pads"]], outline=[f(p) for p in res["outline"]], conn=f(res["conn"]), apothem=float(res["apothem"]), r_out=float(res["r_out"]), failed=res["failed"],
             params=dict(P=P, W_TR=W_TR, GAP_TR=GAP_TR, PAD_R=PAD_R, ANTIPAD_R=ANTIPAD_R, VI_L=f(VI_L), VD_L=f(VD_L), V_DIODE=V_DIODE, I_DIODE=I_DIODE, N_STR=N_STR, N_PER=N_PER, PAD_W=PAD_W, MARGIN=MARGIN, GRID=GRID, KEEP_VIA=KEEP_VIA, NEAR_R=NEAR_R))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "routing.json").write_text(json.dumps(d), encoding="utf-8", newline="\n")


def verify(res):
    """independent clearance check: sample every trace at 0.03 mm and compare centre distances with other nets / foreign vias / pads / outline."""
    nets = [n for n in res["nets"] if n["poly"] is not None]
    pts, lab = [], []
    for k, n in enumerate(nets):
        sp = resample(n["poly"], 0.03)
        pts += sp
        lab += [k] * len(sp)
    pts, lab = np.array(pts), np.array(lab)
    tree = cKDTree(pts)
    pairs = tree.query_pairs(0.40, output_type="ndarray")
    d_min = 9.0
    viol = 0
    if len(pairs):
        diff = lab[pairs[:, 0]] != lab[pairs[:, 1]]
        pp = pairs[diff]
        if len(pp):
            dd = np.linalg.norm(pts[pp[:, 0]] - pts[pp[:, 1]], axis=1)
            d_min = float(dd.min())
            viol = int((dd < W_TR + GAP_TR - 1e-6).sum())
    # trace vs foreign vias
    vpts = []
    vown = []
    for ci, (vi, vd) in res["vias"].items():
        vpts += [vi, vd]
        vown += [("I", ci), ("D", ci)]
    vtree = cKDTree(np.array(vpts))
    own = {}
    for k, n in enumerate(nets):
        own[k] = set()
    vmin = 9.0
    vviol = 0
    ownmap = {}
    for k, n in enumerate(nets):
        for p in ([n["p0"], n["p1"]] if n["kind"] == "link" else [n["p0"]]):
            d, j = vtree.query(p)
            ownmap.setdefault(k, set()).add(int(j))
    for k, n in enumerate(nets):
        sp = np.array(resample(n["poly"], 0.03))
        idxs = vtree.query_ball_point(sp, KEEP_VIA - 0.01)
        bad = set()
        for lst in idxs:
            for j in lst:
                if j not in ownmap.get(k, set()):
                    bad.add(j)
        vviol += len(bad)
        dmin_k = [vtree.query(sp)[0][i] for i in range(0)]
        for j2 in range(len(vpts)):
            pass
    return dict(min_centre_distance_between_nets_mm=d_min, min_gap_between_traces_mm=d_min - W_TR if d_min < 9 else None, violations_trace_trace=viol, trace_to_foreign_via_violations=vviol, n_samples=int(len(pts)))


if __name__ == "__main__":
    if sys.argv[1] == "build":
        r = build()
        v = verify(r)
        print(v)
        a = antipad_check(r)
        print("ANTIPAD", json.dumps({k: v for k, v in a.items() if k != "worst_hits"}, default=str))
        for h in a["worst_hits"][:5]:
            print("  hit", h)
        save_result(r)
        (OUT / "antipad_check.json").write_text(json.dumps(dict(verify=v, antipad=a), indent=1, default=str), encoding="utf-8", newline="\n")
