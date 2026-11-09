"""SUB6 FINAL BOARD, 13-pad bias re-route: 12 drive pads S01..S12 + 1 common RETURN pad.  MODEL_ONLY / NOT_HFSS_VALIDATED.
    SUB6_MESH=96 SUB6_VARIANT=13pad python sub6_final_13pad.py
Everything RF is taken unchanged from results/sub6_final_board/mesh96/routing.json (cells, vias, strings, order).  The 228 link traces (disc via(i) -> island via(i+1)) are kept
as they are; the 24 lead traces are replaced by 12 start leads (island via of position 1 -> pad S_k) and one RETURN tree (disc via of position 20 of every string -> RETURN pad),
RETURN branches 0.5 mm wide where the routing allows (ladder 0.5/0.4/0.3/0.2/0.15).  Pad column: 13 pads 1.0 x 1.0, pitch 2.54, RETURN = middle pad (k = 6, y = 0).
Writes <OUT>/fab/{traces.json, cell_table.csv, drill.csv, drill.drl, string_table.csv, layout_summary.json, sub6_final_board.dxf, connectivity.json, drc.json, rf_unchanged.json}.
Gates: CONNECTIVITY_OK/FAIL, DRC_OK/FAIL, RF_UNCHANGED/RF_CHANGED.
"""
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

os.environ["SUB6_VARIANT"] = "13pad"
os.environ.setdefault("SUB6_MESH", "96")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.ndimage import distance_transform_edt
from scipy.spatial import cKDTree

import sub6_final_layout as L

SRC = L.ROOT / "results" / "sub6_final_board" / "mesh96"
OUT = L.OUT
FAB = OUT / "fab"
R3 = lambda v: round(float(v) + 0.0, 3)
LADDER = [0.5, 0.4, 0.3, 0.2, 0.15]
W_LEAD = 0.15
ANTIPAD_MIN = 0.172
N_PAD = 13
RET_K = 6
PITCH = L.PAD_PITCH
MARGIN_PAINT = 0.06


def load():
    d = json.load(open(SRC / "routing.json"))
    vias = {int(k): (np.array(v["island"]), np.array(v["disc"])) for k, v in d["vias"].items()}
    nets = d["nets"]
    return d, vias, nets


def route_masks(router, ba, bb, p0, p1, cap=14_000_000):
    W, H = router.W, router.H

    def free_node(p):
        cx, cy = router.node(p)
        best = None
        for r in range(0, 6):
            for ddx in range(-r, r + 1):
                for ddy in range(-r, r + 1):
                    nx, ny = cx + ddx, cy + ddy
                    if 0 <= nx < W and 0 <= ny < H and not (ba[ny * W + nx] or bb[ny * W + nx]):
                        d = math.hypot(router.x0 + nx * L.GRID - p[0], router.y0 + ny * L.GRID - p[1])
                        if best is None or d < best[0]:
                            best = (d, nx, ny)
            if best:
                return best[1], best[2]
        return None

    a, b = free_node(p0), free_node(p1)
    if a is None or b is None:
        return None
    if not L.astar(ba, bb, W, H, a[0], a[1], b[0], b[1], router.g, router.parent, router.hk, router.hv, cap):
        return None
    path = []
    cur = b[1] * W + b[0]
    while cur != -1:
        path.append((router.x0 + (cur % W) * L.GRID, router.y0 + (cur // W) * L.GRID))
        cur = router.parent[cur]
    return [tuple(map(float, p0))] + [tuple(map(float, q)) for q in L.simplify(path[::-1])] + [tuple(map(float, p1))]


class Ctx:
    def __init__(self, d, vias, pads):
        self.d, self.vias, self.pads = d, vias, pads
        self.outline = np.array(d["outline"])
        self.router = L.Router(self.outline, 4.5)
        self.vlist = []
        for ci, (vi, vd) in vias.items():
            for v in (vi, vd):
                self.router.add_via(v)
                self.vlist.append(v)
        for p in pads:
            self.router.add_via(p, L.PAD_KEEP)
            self.vlist.append(p)
        self.router.vlist = list(self.vlist)

    def paint(self, poly, w, val):
        self.router.paint(poly, w / 2 + 0.075 + L.GAP_TR + MARGIN_PAINT, val)

    def route(self, p0, p1, w, own_pts=()):
        """route centre line p0 -> p1 with trace width w; vias within NEAR_R of p0/p1 (own vias / own pad) are lifted."""
        r = self.router
        near = [q for q in r.vlist if min(np.linalg.norm(q - np.asarray(p0)), np.linalg.norm(q - np.asarray(p1))) < L.NEAR_R]
        for q in near:
            r.remove_via(q)
        ba = ((r.trace > 0) | (r.static > 0)).astype(np.uint8)
        bb = (r.viac > 0).astype(np.uint8)
        if w > L.W_TR + 1e-9:
            e = (w - L.W_TR) / 2 + 0.01
            bl = (ba | bb).reshape(r.H, r.W).astype(bool)
            dist = distance_transform_edt(~bl) * L.GRID
            dil = (dist <= e) & ~bl
            # via keep-outs get the full dilation; trace/static too (same map)
            ba = (ba.reshape(r.H, r.W) | dil).astype(np.uint8).ravel()
        poly = route_masks(r, ba, bb, p0, p1)
        for q in near:
            r.add_via(q)
        return poly


def build_pads(conn_x):
    return [np.array([conn_x, (RET_K - k) * PITCH]) for k in range(N_PAD)]


def main():
    import pickle
    import random
    t0 = time.time()
    FAB.mkdir(parents=True, exist_ok=True)
    d, vias, nets0 = load()
    conn_x = d["conn"][0]
    pads = build_pads(conn_x)
    ctx = Ctx(d, vias, pads)
    links = [n for n in nets0 if n["kind"] == "link"]
    for n in links:
        ctx.paint([tuple(p) for p in n["poly"]], L.W_TR, +1)
    strings = d["strings"]
    ends = {s["id"]: vias[s["cells"][-1]][1] for s in strings}
    starts = {s["id"]: vias[s["cells"][0]][0] for s in strings}
    print(f"setup {time.time() - t0:.0f}s links {len(links)} grid {ctx.router.W}x{ctx.router.H}", flush=True)

    def route_ladder(p0, p1, ladder=LADDER):
        for w in ladder:
            poly = ctx.route(p0, p1, w)
            if poly is not None:
                return w, poly
        return None

    # ---------------- 1. RETURN ring in the outline margin; reaches the RETURN pad from the west side (corridor between pads and outline)
    ret_pad = pads[RET_K]
    a_l = d["apothem"] - 1.9
    Rl = a_l / math.cos(math.radians(30))
    hexv = [(0.0, Rl), (a_l, Rl / 2), (a_l, -Rl / 2), (0.0, -Rl), (-a_l, -Rl / 2), (-a_l, Rl / 2)]
    north = [(-a_l, 0.0), (-a_l, Rl / 2), hexv[0], hexv[1], hexv[2]]
    south = [(-a_l, 0.0), (-a_l, -Rl / 2), hexv[3], hexv[2]]
    branches = []
    for arm, name in ((north, "N"), (south, "S")):
        prev = ret_pad
        for wp in arm:
            got = route_ladder(prev, np.array(wp))
            if got is None:
                raise SystemExit(f"ring segment {name} {wp} unrouted")
            w, poly = got
            branches.append(dict(string=0, w=w, poly=poly, kind=f"ring{name}"))
            prev = np.array(poly[-1])
        print(f"  ring arm {name} done", flush=True)
    for b in branches:
        ctx.paint(b["poly"], b["w"], +1)
    leads, fails = {}, []
    # ---------------- 2. RETURN branches (Prim, ladder 0.5 ... 0.15) from every string end to the nearest tree point
    tree_pts = [ret_pad] + [np.array(q) for b in branches for q in L.resample(b["poly"], 1.0)]
    remaining = dict(ends)
    while remaining:
        tp = np.array(tree_pts)
        sid = min(remaining, key=lambda s: np.linalg.norm(tp - remaining[s], axis=1).min())
        e = remaining[sid]
        # candidates: the closest point of every tree polyline (ring arms, earlier branches), nearest 8
        cl = []
        for b in branches:
            pts = np.array(L.resample(b["poly"], 1.0))
            i = int(np.argmin(np.linalg.norm(pts - e, axis=1)))
            cl.append((float(np.linalg.norm(pts[i] - e)), pts[i]))
        cand = [c for _, c in sorted(cl, key=lambda x: x[0])[:8]]
        got = None
        for w in LADDER:
            best = None
            for c in cand:
                poly = ctx.route(e, c, w)
                if poly is not None:
                    ln = L.poly_len(poly)
                    if best is None or ln < best[0]:
                        best = (ln, poly)
            if best:
                got = (w, best[1])
                break
        if got is None:
            fails.append(("branch", sid))
            del remaining[sid]
            print(f"  RETURN branch string {sid} UNROUTED", flush=True)
            continue
        w, poly = got
        branches.append(dict(string=sid, w=w, poly=poly, kind="branch"))
        ctx.paint(poly, w, +1)
        tree_pts += [np.array(q) for q in L.resample(poly, 1.0)]
        del remaining[sid]
        print(f"  RETURN branch string {sid}: width {w}, length {L.poly_len(poly):.1f} mm", flush=True)
    # ---------------- 3. S leads (0.15 mm): angular fan assignment (highest angle seen from the connector -> top pad), routed outside-in
    spads = [k for k in range(N_PAD) if k != RET_K]
    ang = sorted(starts, key=lambda s: -math.atan2(starts[s][1], starts[s][0] - conn_x))
    assign = {sid: spads[i] for i, sid in enumerate(ang)}
    pad_order = [0, 12, 1, 11, 2, 10, 3, 9, 4, 8, 5, 7]
    inv = {k: sid for sid, k in assign.items()}
    for k in pad_order:
        sid = inv[k]
        poly = ctx.route(starts[sid], pads[k], W_LEAD)
        if poly is None:
            fails.append((sid, k))
            continue
        leads[sid] = poly
        ctx.paint(poly, W_LEAD, +1)
        print(f"  lead string {sid} -> pad index {k}: {L.poly_len(poly):.1f} mm", flush=True)
    print("lead failures", fails, f"{time.time() - t0:.0f}s", flush=True)
    if any(f[0] != "branch" for f in fails):
        pickle.dump(dict(branches=branches, leads=leads, assign=assign, fails=fails), open(OUT / "partial13.pkl", "wb"))
        return None
    print(f"done {time.time() - t0:.0f}s fails {fails}", flush=True)
    res = dict(d=d, vias=vias, links=links, pads=pads, branches=branches, leads=leads, assign=assign, fails=fails, starts=starts, ends=ends)
    pickle.dump(res, open(OUT / "route13.pkl", "wb"))
    return res


if __name__ == "__main__":
    main()
