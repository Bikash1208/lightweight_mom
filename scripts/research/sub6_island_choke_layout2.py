"""SUB6 island cell: butterfly bias-choke LAYOUT for the new via positions (MODEL_ONLY).  Takes over sub6_island_choke.py / sub6_island_choke_layout.py (commit 32f7c8a, unchanged copies:
sub6_island_choke.py, sub6_island_choke_layout_ref.py).

Design (user, 2026-11-09): FR4 0.8 mm, choke line 0.15 mm (Z0 121.4 ohm, lambda_g/4 = 7.58 mm, meander length l = 7.571 mm from the optimiser), BUTTERFLY = two 120-degree radial stubs, r2 = 3.148 mm, apex at the end
of the choke line, one on each side of the line (axes at +-90 deg to the line direction), bias trace (via pad -> choke line) 0.3 mm wide, 0.5 mm long.  Vias at the |J|-minimum positions (island frame, slot along x):
island (-1.78, -0.53), disc (-3.63, -3.63).  Periodic model: nothing may cross the 11.6 mm cell boundary; trace / gap >= 0.15 mm inside the cell, between the two chokes, to the pad of the other via and to the
periodic images of both chokes.  Anti-pads (r 0.5 in the ground) keep clear of each other and of the cell boundary.
    python sub6_island_choke_layout2.py fit      -> results/sub6_island_choke/layout_fit.json
"""
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sub6_island_choke as ch
import sub6_island_choke_layout_ref as L0

ROOT = ch.ROOT
RES = ROOT / "results" / "sub6_island_choke"
P = 11.6
H_FR4 = 0.8
ANG = 120.0
W_CH, W_BIAS = 0.15, 0.30
LEAD, LINK = 0.50, 0.50
GAP_MIN = 0.15
PAD_R, ANTIPAD_R = 0.30, 0.50
VIA_I, VIA_D = (-1.78, -0.53), (-3.63, -3.63)
DIRS = L0.DIRS
LINKS = (2.9, 3.1, 3.3, 3.6, 4.0, 4.5)  # straight approach to the butterfly apex: the 120-degree sectors reach 0.866 r2 = 2.73 mm back along the line
SPACINGS = (2.4, 1.6, 1.2, 0.8, 0.6, 0.45, 0.3, 0.2, 0.15)


def design():
    d = ch.design(H_FR4, ANG, W_CH)
    # butterfly = two identical sectors in parallel: Z_stub / 2
    fs = np.linspace(ch.BAND[0], ch.BAND[1], 7)
    zin = []
    for f in fs:
        z, zs = ch.z_in(f, H_FR4, W_CH, d["l_mm"], d["r2_mm"], ANG)
        p = ch.ms_params(W_CH, H_FR4, f)
        g = p["alpha_np_m"] + 1j * p["beta0_rad_m"]
        t = np.tanh(g * d["l_mm"] * 1e-3)
        zs2 = zs / 2.0
        zin.append(p["z0"] * (zs2 + p["z0"] * t) / (p["z0"] + zs2 * t))
    d["zin_butterfly_min_ohm"] = float(min(abs(z) for z in zin))
    d["zin_butterfly_f0"] = complex(zin[3])
    return d


def sector(apex, axis_angle, r2, ang_deg, n=16):
    half = math.radians(ang_deg) / 2
    return [apex] + [(apex[0] + r2 * math.cos(axis_angle - half + 2 * half * k / n), apex[1] + r2 * math.sin(axis_angle - half + 2 * half * k / n)) for k in range(n + 1)]


def build(via, dname, flip, s, l_mm, r2, n_runs, link=LINK):
    """parts (in order): [pad, bias trace, (meander runs and joints...), link, stub+, stub-].  Bias trace (0.3 mm) from the via centre, then the 0.15 mm choke line of length l_mm
    (lead excluded) with n_runs runs, then the butterfly at the end."""
    d = np.array(DIRS[dname], float)
    nrm = np.array([-d[1], d[0]]) * flip
    pm = W_CH + s
    ls = (l_mm - link - (n_runs - 1) * pm) / n_runs
    if not (L0.LS_MIN <= ls <= L0.LS_MAX):
        return None
    p = np.array(via, float)
    path = [tuple(p)]
    p = p + d * LEAD
    path.append(tuple(p))
    parts = [L0.circle_poly(via, PAD_R), L0.rect(path[0], path[1], W_BIAS)]
    side = 1.0
    for k in range(n_runs):
        q = p + nrm * side * ls
        parts.append(L0.rect(tuple(p), tuple(q), W_CH))
        path.append(tuple(q))
        p = q
        if k < n_runs - 1:
            q = p + d * pm
            parts.append(L0.rect(tuple(p), tuple(q), W_CH))
            path.append(tuple(q))
            p = q
        side = -side
    q = p + d * link
    parts.append(L0.rect(tuple(p), tuple(q), W_CH))
    path.append(tuple(q))
    a0 = math.atan2(d[1], d[0])
    parts.append(sector(tuple(q), a0 + math.pi / 2, r2, ANG))
    parts.append(sector(tuple(q), a0 - math.pi / 2, r2, ANG))
    return dict(via=via, dir=dname, flip=flip, path=path, parts=parts, apex=tuple(q), n_runs=n_runs, ls=ls, pitch=pm, spacing=s, link=link, w=W_CH, r2=r2, angle=ANG, axis_angle=a0)


def _tag_pairs(c):
    """index pairs exempt from the self-gap test: neighbours in the chain and the two stubs with the link / each other."""
    n = len(c["parts"])
    ex = {(i, i + 1) for i in range(n - 3)}
    ex |= {(n - 3, n - 2), (n - 3, n - 1), (n - 2, n - 1)}
    ex |= {(0, 1), (0, 2)}  # pad touches the bias trace and the first run
    ex |= {(1, 2)}
    return ex


def self_gap(c):
    pr, ex = c["parts"], _tag_pairs(c)
    g = 1e9
    for i in range(len(pr)):
        for j in range(i + 1, len(pr)):
            if (i, j) in ex or (i, j) == (0, 3 if len(pr) > 4 else 99):
                continue
            g = min(g, L0.poly_dist(pr[i], pr[j]))
    return g


def gaps(chokes):
    out = dict(between=1e9, images=1e9, self=1e9, outside_cell=0.0)
    for a, b in itertools.combinations(chokes, 2):
        out["between"] = min(out["between"], min(L0.poly_dist(pa, pb) for pa in a["parts"] for pb in b["parts"]))
    for a in chokes:
        for b in chokes:
            for ix, iy in itertools.product((-1, 0, 1), repeat=2):
                if (ix, iy) == (0, 0):
                    continue
                out["images"] = min(out["images"], min(L0.poly_dist(pa, pb) for pa in a["parts"] for pb in L0.shift(b["parts"], ix * P, iy * P)))
        out["self"] = min(out["self"], self_gap(a))
        xs = [x for poly in a["parts"] for x, _ in poly]
        ys = [y for poly in a["parts"] for _, y in poly]
        out["outside_cell"] = max(out["outside_cell"], max(max(xs) - P / 2, -P / 2 - min(xs), max(ys) - P / 2, -P / 2 - min(ys), 0.0))
    return out


def antipads_ok():
    d = math.hypot(VIA_I[0] - VIA_D[0], VIA_I[1] - VIA_D[1])
    edge = min(P / 2 - abs(v) for via in (VIA_I, VIA_D) for v in via) - ANTIPAD_R
    return dict(antipad_to_antipad_mm=d - 2 * ANTIPAD_R, antipad_to_cell_edge_mm=edge, ok=bool(d - 2 * ANTIPAD_R >= GAP_MIN and edge >= GAP_MIN))


def fit():
    RES.mkdir(parents=True, exist_ok=True)
    des = design()
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in des.items() if k in ("r1_mm", "r2_mm", "l_mm", "zin_min_ohm", "zin_butterfly_min_ohm", "Yload_max_uS", "via_L_nH")}, "Z0 %.1f" % des["ms"]["z0"], flush=True)
    apd = antipads_ok()
    print("anti-pads", apd, flush=True)
    best = None
    tried = 0
    for s in SPACINGS:
        sing = {}
        for name, via in (("I", VIA_I), ("D", VIA_D)):
            sing[name] = []
            for dn, fl, n, lk in itertools.product(DIRS, (1, -1), range(1, 9), LINKS):
                c = build(via, dn, fl, s, des["l_mm"], des["r2_mm"], n, lk)
                if c is None:
                    continue
                g = gaps([c])
                tried += 1
                if g["outside_cell"] <= 1e-9 and min(g["self"], g["images"]) >= GAP_MIN - 1e-9:
                    sing[name].append((g, c))
        print(f"s = {s}: single-choke candidates inside the cell, island {len(sing['I'])}, disc {len(sing['D'])}", flush=True)
        ok = []
        for (ga, a), (gb, b) in itertools.product(sing["I"], sing["D"]):
            g = gaps([a, b])
            if g["outside_cell"] <= 1e-9 and min(g["between"], g["images"], g["self"]) >= GAP_MIN - 1e-9:
                ok.append((g, [a, b]))
        if ok:
            ok.sort(key=lambda t: -min(t[0]["between"], t[0]["images"], t[0]["self"]))
            best = (s, ok[0], len(ok))
            break
    res = dict(design={k: v for k, v in des.items() if k not in ("ms",)}, ms=des["ms"], antipads=apd, vias=dict(island=VIA_I, disc=VIA_D), widths=dict(choke=W_CH, bias=W_BIAS), n_tried=tried,
               fits=best is not None)
    if best:
        s, (g, chokes), n_ok = best
        res.update(spacing_mm=s, n_layouts_ok=n_ok, gaps=g, chokes=[dict(via=c["via"], dir=c["dir"], flip=c["flip"], n_runs=c["n_runs"], ls_mm=c["ls"], apex=c["apex"], axis_angle_rad=c["axis_angle"], link=c["link"], spacing=c["spacing"]) for c in chokes])
        print("FIT", s, g, [(c["dir"], c["flip"], c["n_runs"], round(c["ls"], 2), c["link"]) for c in chokes], flush=True)
    else:
        print("NO FIT with trace/gap >= 0.15 mm inside the 11.6 mm cell", flush=True)
    (RES / "layout_fit.json").write_text(json.dumps(res, indent=1, default=lambda o: [o.real, o.imag] if isinstance(o, complex) else float(o)))
    return res


if __name__ == "__main__":
    {"fit": fit}[sys.argv[1]]()
