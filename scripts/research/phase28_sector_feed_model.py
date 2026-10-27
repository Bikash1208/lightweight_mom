"""Phase 28 (model side): scalar 1-bit cells fed by a SECTOR source (a point launcher with a 78 deg half-power beamwidth, the Yagi-type launcher of [78]) at the edge of the 150 mm hexagon,
pitch 6.4 mm and 9.0 mm. IDEALISED_CELLS, MODEL_ONLY.

    python scripts/research/phase28_sector_feed_model.py <out.json> [--workers N]

Source (ASSUMPTION, not from the paper): at the middle of the -x flat of the hexagon, field pattern f(psi) = cos(psi)**1.373 (power cos**2.746 = 78 deg half-power width) for |psi| < 90 deg, cylindrical spreading,
phase exp(-j beta rho), amplitude f(psi) exp(-alpha rho) / sqrt(rho) (rho >= 6.5 mm), current along the local propagation direction (radial from the source, the TM0 field), beta = 1.053 k0,
ONE uniform leakage alpha chosen so that the converted fraction of the launched sector power is 0.504 (same conversion as the monopole and plane-feed studies). Cell state: m = cos(phi_ref - k_t . r) quantised to
+-1 (1 bit). Gain bookkeeping D_ap + 10 log10(0.592 eta) (placeholder launch fraction). Co-/cross-polarisation: Ludwig-3 x-component / y-component at the peak. [78]'s scan limits (-30/+60 in xoz, +-45 in yoz) are included.
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from scipy.optimize import brentq  # noqa: E402

from src.hologram import coded_feasibility as cf  # noqa: E402

cf.set_frequency(10.0e9)
R_MM = 150.0
BETA_OVER_K0 = 1.053
HPBW_DEG = 78.0
P_FIELD = math.log(0.5) / (2.0 * math.log(math.cos(math.radians(HPBW_DEG / 2))))  # field exponent: cos**p, power cos**2p
ETA_TARGET = 0.504
RHO_MIN = 6.5e-3
PITCHES_MM = [6.4, 9.0]
# (theta0, phi0): 0/30/45 in the two planes, plus the scan limits of [78]: -30/+60 (xoz), +-45 (yoz)
BEAMS = [(0.0, 0.0), (30.0, 0.0), (45.0, 0.0), (60.0, 0.0), (30.0, 180.0), (30.0, 90.0), (45.0, 90.0), (45.0, 270.0)]


def source_xy():
    poly = cf.hexagon_m(R_MM * 1e-3)
    return poly, (-R_MM * 1e-3 * math.cos(math.radians(30.0)), 0.0)


def chord_length(poly, src, psi):
    """Distance from `src` along direction `psi` to the polygon boundary (array)."""
    dx, dy = np.cos(psi), np.sin(psi)
    best = np.full(np.shape(psi), np.inf)
    for a, b in zip(poly, np.roll(poly, -1, axis=0), strict=True):
        ex, ey = b[0] - a[0], b[1] - a[1]
        den = dx * ey - dy * ex
        with np.errstate(divide="ignore", invalid="ignore"):
            t = ((a[0] - src[0]) * ey - (a[1] - src[1]) * ex) / den
            s = ((a[0] - src[0]) * dy - (a[1] - src[1]) * dx) / den
        ok = (t > 1e-9) & (s >= -1e-9) & (s <= 1 + 1e-9)
        best = np.where(ok, np.minimum(best, t), best)
    return best


def alpha_for_eta(eta):
    poly, src = source_xy()
    psi = np.linspace(-math.pi / 2 + 1e-4, math.pi / 2 - 1e-4, 2001)
    w = np.cos(psi) ** (2 * P_FIELD)
    L = chord_length(poly, src, psi)
    f = lambda a: float(np.sum(w * (1 - np.exp(-2 * a * L))) / np.sum(w)) - eta  # noqa: E731
    return brentq(f, 1e-3, 500.0), float(np.sum(w * (L < 1.0)) / np.sum(w))


def case(job):
    pitch_mm, beam = job
    pitch = pitch_mm * 1e-3
    poly, src = source_xy()
    green = cf.green_function()
    alpha, _ = alpha_for_eta(ETA_TARGET)
    x, y = cf.lattice(poly, pitch, keepout_m=0.0)
    rho = np.hypot(x - src[0], y - src[1])
    keep = rho >= RHO_MIN
    x, y, rho = x[keep], y[keep], rho[keep]
    psi = np.arctan2(y - src[1], x - src[0])
    f = np.maximum(np.cos(psi), 0.0) ** P_FIELD
    amp = f * np.exp(-alpha * rho) / np.sqrt(rho / RHO_MIN)
    beta = BETA_OVER_K0 * cf.K0
    phi_ref = beta * rho
    th0, ph0 = math.radians(beam[0]), math.radians(beam[1])
    kt = cf.K0 * math.sin(th0)
    m = np.cos(phi_ref - kt * (x * math.cos(ph0) + y * math.sin(ph0)))
    _, mq = cf.lp_states(m, 1)
    j = amp * np.exp(-1j * phi_ref) * mq
    jx, jy = j * np.cos(psi), j * np.sin(psi)
    r = cf.analyse(jx, jy, x, y, green, pitch, "lp", beam)
    eta = ETA_TARGET
    th, ph = np.array([r["peak_theta_deg"]]), np.array([r["peak_phi_deg"]])
    et, ep = cf.far_field(jx, jy, x, y, green, th, ph, pitch)
    p = math.radians(ph[0])
    co = et[0] * math.cos(p) - ep[0] * math.sin(p)
    cr = et[0] * math.sin(p) + ep[0] * math.cos(p)
    # strongest sidelobe farther than 15 deg from the peak, and its position (coarse 2 x 4 deg grid)
    tg, pg = np.meshgrid(np.arange(0.0, 90.0 + 1e-9, 2.0), np.arange(0.0, 360.0, 4.0), indexing="ij")
    e1, e2 = cf.far_field(jx, jy, x, y, green, tg, pg, pitch)
    pw = np.abs(e1) ** 2 + np.abs(e2) ** 2
    far = np.where(cf._angle_deg(tg, pg, r["peak_theta_deg"], r["peak_phi_deg"]) > 15.0, pw, -1.0)
    i = np.unravel_index(np.argmax(far), far.shape)
    return {
        "pitch_mm": pitch_mm,
        "beam": list(beam),
        "n_cells": int(x.size),
        "n_diodes": int(x.size),
        "alpha_np_per_m": alpha,
        "eta_conv": eta,
        "d_ap_dbi": r["d_ap_dbi"],
        "gain_dbi": cf.launch_gain_dbi(r["d_ap_dbi"], eta),
        "peak_theta_deg": r["peak_theta_deg"],
        "peak_phi_deg": r["peak_phi_deg"],
        "pointing_error_deg": r["pointing_error_deg"],
        "sll_db": r["sll_db"],
        "sidelobe_theta_phi_deg": [float(tg[i]), float(pg[i])],
        "xpd_at_peak_db": float(20 * math.log10(abs(co) / max(abs(cr), 1e-30))),
    }


def main(out, workers):
    jobs = [(p, b) for p in PITCHES_MM for b in BEAMS]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        rows = list(ex.map(case, jobs))
    alpha, _ = alpha_for_eta(ETA_TARGET)
    doc = {
        "labels": ["MODEL_ONLY", "IDEALISED_CELLS", "NOT_HFSS_VALIDATED"],
        "frequency_ghz": 10.0,
        "source": {"hpbw_deg": HPBW_DEG, "field_exponent": P_FIELD, "alpha_np_per_m_for_eta": alpha, "eta_target": ETA_TARGET, "note": "ASSUMED sector source, not from [78]"},
        "rows": rows,
    }
    Path(out).write_text(json.dumps(doc, indent=1), encoding="utf-8")
    print("wrote", out, len(rows), "alpha", alpha)


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], int(a[a.index("--workers") + 1]) if "--workers" in a else 16)
