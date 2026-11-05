"""Segmented comparator, real model: 4 probes + SP4T, quarter-aperture holograms, RT5880 chosen cell (alpha 1.84 Np/m), MODEL_ONLY.
State k uses the wedge W_k of the hexagon (90 deg, centred on phi_k) fed by its own probe at the wedge centroid (or a fraction f of it); only the cells of W_k radiate (ideal isolation of the other wedges: optimistic upper bound),
static cells only (no diodes), same cell options and power bookkeeping as the coded board with the radial geometry measured from the probe.  Polarisation per state as in set B (RHCP/LHCP/RHCP/LHCP). The SP4T loss is a parameter:
gain = geometric gain - IL; launch fraction 0.627 as for the coded board."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_arch_study as st
import sidequest_coded_improve as I
import sidequest_coded_k4 as K


def clip_halfplane(poly, a, b, c):
    """Keep a x + b y + c >= 0 (Sutherland-Hodgman, convex)."""
    out = []
    for i, p in enumerate(poly):
        q = poly[(i + 1) % len(poly)]
        fp, fq = a * p[0] + b * p[1] + c, a * q[0] + b * q[1] + c
        if fp >= 0:
            out.append(tuple(p))
        if fp * fq < 0:
            t = fp / (fp - fq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def wedge_poly(phi, half=np.radians(45.0)):
    hexa = [tuple(v) for v in st.hexagon()]
    # half-planes through the origin: angle within [phi - half, phi + half]
    for sgn, ang in ((1, phi - half), (-1, phi + half)):
        a, b = -np.sin(ang) * sgn, np.cos(ang) * sgn  # n . r >= 0 with n perpendicular to the edge direction
        hexa = clip_halfplane(hexa, a, b, 0.0)
    return np.array(hexa)


def centroid(poly):
    x, y = poly[:, 0], poly[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    a = cr.sum() / 2
    return np.array([((x + x1) * cr).sum() / (6 * a), ((y + y1) * cr).sum() / (6 * a)]), abs(a)


class QBoard(K.KBoard):
    def __init__(self, phi_deg, pol, feed_frac=1.0, half_deg=45.0):
        phi = np.radians(phi_deg)
        poly = wedge_poly(phi, np.radians(half_deg))
        cen, area = centroid(poly)
        self.feed = cen * feed_frac
        self.area = area
        _, x, y, _ = st.aperture("hex")
        ang = np.arctan2(y, x)
        d = (ang - phi + np.pi) % (2 * np.pi) - np.pi
        inside = (np.abs(d) <= np.radians(half_deg) + 1e-9)
        x, y = x[inside], y[inside]
        rx, ry = x - self.feed[0], y - self.feed[1]
        keep = np.hypot(rx, ry) >= st.RHO_KEEP
        self.x, self.y = x[keep], y[keep]
        self.n = self.x.size
        self.phis, self.pols, self.K = np.array([phi]), [pol], 1
        rx, ry = self.x - self.feed[0], self.y - self.feed[1]
        self.rho, self.psi = np.hypot(rx, ry), np.arctan2(ry, rx)
        self.rhat = np.array([np.cos(self.psi), np.sin(self.psi)])
        self.redge = np.array([float(np.ravel(st.polygon_ray_distance(poly, tuple(self.feed), a))[0]) for a in self.psi])
        self.az = ((self.psi + np.pi) / (2 * np.pi) * 12).astype(int) % 12
        self.arr = [np.exp(1j * st.K0 * np.sin(K.THETA) * (self.x * np.cos(phi) + self.y * np.sin(phi)))]
        self.profile = I.Board.profile.__get__(self)


def span(bd):
    """Angular extent [rad] of the wedge as seen from its probe."""
    d = (bd.psi - np.angle(np.exp(1j * bd.psi).mean()) + np.pi) % (2 * np.pi) - np.pi
    return float(d.max() - d.min())


if __name__ == "__main__":
    STAT, PAIRS = K.library()
    vs = np.linspace(STAT.min(), STAT.max(), 9)
    out = {}
    pols = ["RHCP", "LHCP", "RHCP", "LHCP"]
    for frac in (0.3, 0.5, 0.75, 1.0, 1.3, 1.5):
        rows = []
        try:
            QBoard(0, "RHCP", feed_frac=frac).profile(np.zeros(1), np.zeros(1)) if False else None
        except Exception:  # noqa: BLE001
            pass
        for k, phi in enumerate((0, 90, 180, 270)):
            bd = QBoard(phi, pols[k], feed_frac=frac)
            best = None
            for vset in (np.array([-1.0, 1.0]), vs):
                for g in (1.0, 2.0):
                    try:
                        r, _ = K.run(bd, K.options(1, vset, []), [0.0], lam=0.0, g=g, iters=3)
                    except ValueError:
                        continue
                    if best is None or r["states"][0]["gain"] > best[0]["states"][0]["gain"]:
                        best = (r, g, len(vset))
            if best is None:
                rows = None
                break
            r = best[0]
            s = r["states"][0]
            dq = 10 * np.log10(4 * np.pi * bd.area / K.LAM0**2)  # quadrant directivity limit
            th = span(bd)
            rows.append({"phi": phi, "pol": pols[k], "cells": int(bd.n), "area_mm2": bd.area * 1e6, "feed_xy_mm": (bd.feed * 1e3).round(1).tolist(), "gain_geom": s["gain"], "SLL": s["SLL"], "XPD": s["XPD"], "AR": s["AR"],
                         "ap_eff_quadrant": 10 ** ((s["gain"] - dq) / 10), "edge": r["edge"], "SW_to_beam": r["SW_to_beam"], "span_deg": np.degrees(th),
                         "levels": best[2], "g": best[1]})
        if rows is None:
            print(frac, "not feasible (empty azimuth sector)")
            continue
        out[f"feed_{frac}"] = rows
        print(frac, [(round(x["gain_geom"], 2), round(x["SLL"], 1), round(x["SW_to_beam"], 2), round(x["edge"], 2)) for x in rows], flush=True)
    Path("results/sidequest_hybrid/quadrant_board.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
