"""Phase 26: tile-amplitude equalisation of board B (Butler-ready): flatten the aperture amplitude across the shared edges.

    python scripts/research/phase26_tile_equalise.py <out.json> [--quick]

The tile-period (259.8 mm) amplitude modulation of the four hexagon apertures puts grating lobes ~6 deg from the beam.
Candidate designs change only the TARGET DENSITY of the march (cell choices stay table-limited): `b` boosts the
density with the distance from the own feed, g *= 1 + b (rho/rho_edge)^2 (the wave is depleted outward, so the
aperture amplitude falls toward every edge); `w, g0` taper only toward the FREE (outer) edges: g *= g0 + (1-g0) smoothstep(d_free / w),
no taper toward the shared edges. Model only (main patterns, no pickup: the pickup is a separate -23 dB term).
Metric: strongest RHCP local maximum 4-20 deg from the peak ("tile lobe") relative to the peak, for equal phase and the four
Butler inputs (outputs = ports 2, 3, 1, 4; steps -45, +135, -135, +45), anchored gain as in `phase26_fourhex_design.py`.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import maximum_filter

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase26_fourhex_design as fd  # noqa: E402
import phase26_tensor_cp_hexagon as cp  # noqa: E402
from src.hologram.tensor_leaky_cp import SupercellTable  # noqa: E402
from src.pattern.vector_aperture_pattern import copol_component  # noqa: E402

APOTHEM = 0.150 * np.cos(np.radians(30.0))
ORDER = [2, 3, 1, 4]
STEPS = {"equal": 0.0, "1L": -45.0, "2R": 135.0, "2L": -135.0, "1R": 45.0}
D_X = 0.129904


def free_edge_normals_deg(j):
    shared = []
    for q, f in fd.FEEDS.items():
        if q == j:
            continue
        v = np.array(f) - np.array(fd.FEEDS[j])
        if abs(np.hypot(*v) - 259.8) < 1.0:
            shared.append(np.degrees(np.arctan2(v[1], v[0])) % 360.0)
    return [n for n in (0.0, 60.0, 120.0, 180.0, 240.0, 300.0) if all(abs(((n - s + 180.0) % 360.0) - 180.0) > 1.0 for s in shared)]


def factory(b, w_mm, g0):
    def make(j):
        free = free_edge_normals_deg(j)

        def fn(psi, rho, r_edge):
            x = rho[None, :] / r_edge[:, None]
            g = 1.0 + b * x**2
            if w_mm > 0:
                px, py = rho[None, :] * np.cos(psi)[:, None], rho[None, :] * np.sin(psi)[:, None]
                d = np.min([APOTHEM - (px * np.cos(np.radians(n)) + py * np.sin(np.radians(n))) for n in free], axis=0)
                s = np.clip(d / (w_mm * 1e-3), 0.0, 1.0)
                g = g * (g0 + (1.0 - g0) * s * s * (3.0 - 2.0 * s))
            return g

        return fn

    return make


def evaluate(bd, green, theta_c, phi_c, theta_r, phi_r):
    """Per-hexagon fields on a coarse full-hemisphere grid (power) and a fine regional grid (lobes)."""
    out = {}
    for j in fd.FEEDS:
        xg, yg = bd[j]["x_mm"] * 1e-3, bd[j]["y_mm"] * 1e-3
        d = bd[j]["design"]
        ec = cp.pattern(d, xg, yg, green, theta_c, phi_c)
        er = cp.pattern(d, xg, yg, green, theta_r, phi_r)
        eta = bd[j]["march"]["eta_conv"]
        sc = np.sqrt(fd.LAUNCH * eta / fd.prad(*ec, theta_c, phi_c))
        out[j] = (er[0] * sc, er[1] * sc, eta)
    return out


def lobes(et, ep, theta, phi, nominal):
    rh = np.abs(copol_component(et, ep, "rhcp")) ** 2 * 4 * np.pi * 10 ** (fd.OFFSET_DB / 10)
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    t0, p0 = np.radians(nominal[0]), np.radians(nominal[1])
    cosd = np.cos(np.radians(tg)) * np.cos(t0) + np.sin(np.radians(tg)) * np.sin(t0) * np.cos(np.radians(pg) - p0)
    near = np.degrees(np.arccos(np.clip(cosd, -1, 1))) <= 12.0
    i = np.unravel_index(np.argmax(np.where(near, rh, -1)), rh.shape)
    cd = np.cos(np.radians(tg)) * np.cos(np.radians(tg[i])) + np.sin(np.radians(tg)) * np.sin(np.radians(tg[i])) * np.cos(np.radians(pg) - np.radians(pg[i]))
    dpk = np.degrees(np.arccos(np.clip(cd, -1, 1)))
    mx = (rh == maximum_filter(rh, size=(3, 3), mode="nearest")) & (rh > 0)
    out = {}
    for tag, lo, hi in (("tile_lobe_4_20deg", 4.0, 20.0), ("beyond_20deg", 20.0, 90.0)):
        m = mx & (dpk >= lo) & (dpk < hi)
        out[tag] = float(10 * np.log10(rh[m].max() / rh[i])) if m.any() else -99.0
    return {"peak_theta_deg": float(tg[i]), "peak_phi_deg": float(pg[i]), "gain_dbi": float(10 * np.log10(rh[i])), **out}


def butler(fields, theta, phi, cases=STEPS):
    res = {}
    u0, v0 = np.sin(np.radians(30)) * np.cos(np.radians(225)), np.sin(np.radians(30)) * np.sin(np.radians(225))
    for name, step in cases.items():
        et = sum(0.5 * np.exp(1j * np.radians(n * step)) * fields[p][0] for n, p in enumerate(ORDER))
        ep = sum(0.5 * np.exp(1j * np.radians(n * step)) * fields[p][1] for n, p in enumerate(ORDER))
        u1 = u0 - np.radians(step) / (cp.K0 * D_X)
        nominal = (float(np.degrees(np.arcsin(np.hypot(u1, v0)))), float(np.degrees(np.arctan2(v0, u1)) % 360.0))
        res[name] = lobes(et, ep, theta, phi, nominal)
    return res


def main(out, quick=False):
    table = SupercellTable.from_files(cp.TABLE_GLOB)
    green = cp.LayeredGreenFunction(cp.SUB, cp.F, backing="pec")
    theta_c, phi_c = np.arange(0.0, 90.0 + 1e-9, 2.0), np.arange(0.0, 360.0, 4.0)
    theta_r, phi_r = np.arange(5.0, 60.0 + 1e-9, 1.0), np.arange(180.0, 271.0, 1.0)
    cands = [("baseline taper_c 0.5", dict(taper_c=0.5, b=0.0, w=0.0, g0=1.0))]
    grid = [(0.0, 0.0), (0.0, 40.0), (0.0, 80.0), (1.0, 0.0), (2.0, 0.0), (4.0, 0.0), (2.0, 40.0), (4.0, 40.0), (8.0, 40.0), (4.0, 80.0)]
    if quick:
        grid = grid[:2]
    for b, w in grid:
        cands.append((f"b {b:g} free-edge ramp {w:g} mm", dict(taper_c=0.0, b=b, w=w, g0=0.3)))
    res = []
    for name, p in cands:
        bd = fd.design_board("B", table, taper_c=p["taper_c"], eta_design=0.9, density_factory=factory(p["b"], p["w"], p["g0"]))
        fl = evaluate(bd, green, theta_c, phi_c, theta_r, phi_r)
        rec = {"name": name, "params": p, "eta_conv": float(np.mean([fl[j][2] for j in fd.FEEDS])), "butler": butler(fl, theta_r, phi_r),
               "alone_gain_dbi": {str(j): lobes(fl[j][0], fl[j][1], theta_r, phi_r, (30.0, 225.0))["gain_dbi"] for j in fd.FEEDS}}
        res.append(rec)
        b_ = rec["butler"]
        print(f"{name:34s} eta {rec['eta_conv']:.3f} | eq {b_['equal']['gain_dbi']:.1f} dBi lobe {b_['equal']['tile_lobe_4_20deg']:.1f} | "
              f"1L {b_['1L']['tile_lobe_4_20deg']:.1f} 1R {b_['1R']['tile_lobe_4_20deg']:.1f} | 2R {b_['2R']['tile_lobe_4_20deg']:.1f} 2L {b_['2L']['tile_lobe_4_20deg']:.1f}", flush=True)
        Path(out).write_text(json.dumps({"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "SOLVER_SUPERCELL"], "candidates": res}, indent=1))




def board_taper_factory(p, ax_mm=330.0, ay_mm=270.0, floor=0.0):
    """Board-wide smooth taper in GLOBAL coordinates: g = max(cos(pi/2 r)^p, floor), r = sqrt((x/ax)^2 + (y/ay)^2)."""
    def make(j):
        f = np.array(fd.FEEDS[j]) * 1e-3

        def fn(psi, rho, r_edge):
            x = (f[0] + rho[None, :] * np.cos(psi)[:, None]) * 1e3
            y = (f[1] + rho[None, :] * np.sin(psi)[:, None]) * 1e3
            r = np.clip(np.sqrt((x / ax_mm) ** 2 + (y / ay_mm) ** 2), 0.0, 1.0)
            return np.maximum(np.cos(0.5 * np.pi * r) ** p, floor)

        return fn

    return make


def ideal_reference(table, green, theta_r, phi_r, bd):
    """Uniform and cos-tapered IDEAL circular aperture on the same cells: the board outline's own lobe level."""
    kt = cp.K0 * np.sin(np.radians(30.0))
    u = np.array([np.cos(np.radians(225.0)), np.sin(np.radians(225.0))])
    xx = np.concatenate([bd[j]["x_mm"] for j in fd.FEEDS])
    yy = np.concatenate([bd[j]["y_mm"] for j in fd.FEEDS])
    ph = np.exp(-1j * kt * (xx * 1e-3 * u[0] + yy * 1e-3 * u[1]))
    out = {}
    for name, a in (("uniform", np.ones_like(xx)), ("cos_taper_ellipse_330x270", np.cos(0.5 * np.pi * np.clip(np.sqrt((xx / 330.0) ** 2 + (yy / 270.0) ** 2), 0, 1)))):
        e = cp.pattern({"jx": a * ph, "jy": -1j * a * ph}, xx * 1e-3, yy * 1e-3, green, theta_r, phi_r)
        out[name] = lobes(e[0], e[1], theta_r, phi_r, (30.0, 225.0))["tile_lobe_4_20deg"]
    return out


def curated(out):
    table = SupercellTable.from_files(cp.TABLE_GLOB)
    green = cp.LayeredGreenFunction(cp.SUB, cp.F, backing="pec")
    tc, pc = np.arange(0.0, 90.0 + 1e-9, 2.0), np.arange(0.0, 360.0, 4.0)
    tr, pr = np.arange(5.0, 60.0 + 1e-9, 1.0), np.arange(180.0, 271.0, 1.0)
    cands = [("baseline (eta_d 0.9, taper 0.5)", dict(taper_c=0.5, eta=0.9)),
             ("uniform target density (taper 0)", dict(taper_c=0.0, eta=0.9)),
             ("outward boost b=4", dict(taper_c=0.0, eta=0.9, b=4.0)),
             ("free-edge ramp 60 mm, g0 0.3", dict(taper_c=0.0, eta=0.9, w=60.0)),
             ("low density eta_d 0.4", dict(taper_c=0.0, eta=0.4)),
             ("low density eta_d 0.3", dict(taper_c=0.0, eta=0.3)),
             ("board cos taper p2 (330x270), eta_d 0.9", dict(taper_c=0.0, eta=0.9, board=(2, 330.0, 270.0))),
             ("board cos taper p2 (330x270), eta_d 0.3", dict(taper_c=0.0, eta=0.3, board=(2, 330.0, 270.0))),
             ("quiet zone 50 mm (no radiation near the pin), taper 0", dict(taper_c=0.0, eta=0.9, quiet=50.0)),
             ("quiet zone 50 mm, taper 0.5", dict(taper_c=0.5, eta=0.9, quiet=50.0))]
    res, ref = [], None
    for name, p in cands:
        if "board" in p:
            fac = board_taper_factory(*p["board"])
        else:
            fac = factory(p.get("b", 0.0), p.get("w", 0.0), 0.3)
        bd = fd.design_board("B", table, taper_c=p["taper_c"], eta_design=p["eta"], density_factory=fac)
        if ref is None:
            ref = ideal_reference(table, green, tr, pr, bd)
        fl, rem = {}, []
        for j in fd.FEEDS:
            d = dict(bd[j]["design"])
            s = (d["rho"] * 1e3 >= p.get("quiet", 0.0)).astype(float)
            frac = float(np.sum(d["amp"] ** 2 * (1 - s)) / np.sum(d["amp"] ** 2))
            rem.append(frac)
            d["jx"], d["jy"] = d["jx"] * s, d["jy"] * s
            x, y = bd[j]["x_mm"] * 1e-3, bd[j]["y_mm"] * 1e-3
            ec, er = cp.pattern(d, x, y, green, tc, pc), cp.pattern(d, x, y, green, tr, pr)
            eta_eff = bd[j]["march"]["eta_conv"] * (1.0 - frac)
            sc = np.sqrt(fd.LAUNCH * eta_eff / fd.prad(*ec, tc, pc))
            fl[j] = (er[0] * sc, er[1] * sc, eta_eff)
        rec = {"name": name, "params": p, "eta_radiated": float(np.mean([fl[j][2] for j in fd.FEEDS])), "removed_by_quiet_zone": float(np.mean(rem)),
               "butler": butler(fl, tr, pr)}
        res.append(rec)
        b_ = rec["butler"]
        print(f"{name:55s} eta {rec['eta_radiated']:.3f} | eq {b_['equal']['gain_dbi']:.1f}/{b_['equal']['tile_lobe_4_20deg']:.1f} | "
              f"1L {b_['1L']['gain_dbi']:.1f}/{b_['1L']['tile_lobe_4_20deg']:.1f} 1R {b_['1R']['gain_dbi']:.1f}/{b_['1R']['tile_lobe_4_20deg']:.1f} | "
              f"2R {b_['2R']['gain_dbi']:.1f}/{b_['2R']['tile_lobe_4_20deg']:.1f} 2L {b_['2L']['gain_dbi']:.1f}/{b_['2L']['tile_lobe_4_20deg']:.1f}", flush=True)
    # ideal uniform-amplitude, ideal-phase Butler reference on the baseline phases (amplitude removed): steering quantisation floor
    bd = fd.design_board("B", table)
    fl = {}
    for j in fd.FEEDS:
        d = dict(bd[j]["design"])
        s = np.where(d["amp"] > 0, 1.0, 0.0) / np.maximum(d["amp"], 1e-12)
        d["jx"], d["jy"] = d["jx"] * s, d["jy"] * s
        e = cp.pattern(d, bd[j]["x_mm"] * 1e-3, bd[j]["y_mm"] * 1e-3, green, tr, pr)
        fl[j] = (e[0], e[1], 1.0)
    uni = butler(fl, tr, pr)
    Path(out).write_text(json.dumps({"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "SOLVER_SUPERCELL"], "metric": "strongest RHCP local max 4-20 deg from the peak, dB re peak; main patterns only (no pickup)",
                                     "ideal_aperture_reference_lobe_db": ref, "uniform_amplitude_same_phases_butler": {k: v["tile_lobe_4_20deg"] for k, v in uni.items()},
                                     "candidates": res}, indent=1))


if __name__ == "__main__":
    if "--scan" in sys.argv:
        main(sys.argv[1], "--quick" in sys.argv)
    else:
        curated(sys.argv[1])
