"""Phase 28 (loop cell R4): sector-fed 1-bit hologram with the TRUE per-cell tensor response, unrotated vs rotated cells, 9.0 mm pitch, RT/duroid 5880 3.175 mm, 10.0 GHz. MODEL_ONLY, IDEALISED, NOT_HFSS_VALIDATED.

    python scripts/research/phase28_loop_sector_model.py <out.json> [--workers N]

Per-cell sheet tensors (ON, OFF; normal-incidence MoM, mesh 90, `results/phase28/f10/ring_cell/confirm.jsonl`): mean Ybar = (Y_on + Y_off)/2, modulation dY = (Y_on - Y_off)/2, state s = +-1 (1 bit).
Surface wave at the cell: propagation direction e = (cos psi, sin psi) from the sector source (same source and pattern as `phase28_sector_feed_model.py`).
* UNROTATED (diode axis along x everywhere): phase constant beta(psi) from the mean sheet along e, e^T Ybar e; induced modulation current J = s A exp(-j beta rho) dY e / dY_xx (a vector: amplitude and cross-polar part from the tensor);
  local leakage alpha_loc(psi) = alpha0 |e^T dY e / dY_xx|^2 (first-order scaling of the exact alpha0 of the aligned cell: ASSUMPTION).
* ROTATED (each cell turned so that its diode axis lies along e): the tensor is aligned with the wave, so beta = beta(0), alpha_loc = alpha0, J along e. The 40 deg rotation needed at the sector edge would make the 6.6 mm loops
  overlap their 9 mm neighbours (corner reach 4.65 mm > 4.5 mm); this is an upper bound, NOT a realisable layout (see the module docstring of `phase28_loop_rotation.py`).
alpha0 = exact sheet-model leakage of the aligned cell (`ring_cell/report.json`, mesh 90, 10 GHz). The converted fraction eta follows from alpha_loc along each ray; gain = D_ap + 10 log10(0.592 eta) (placeholder launch fraction).
Also reported: the fraction of cells whose 1-bit state is aliased (modulation phase advance per cell step above pi) = the sampling limit at 9 mm.
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

import phase28_sector_feed_model as sf  # noqa: E402

from src.dispersion.sheet_model import find_beta_from_sheet_admittance  # noqa: E402
from src.hologram import coded_feasibility as cf  # noqa: E402
from src.materials import Substrate  # noqa: E402

SUB = Substrate("rt5880_3p175", 2.2, 0.0009, 3.175e-3)
PITCH = 9.0e-3
BEAMS = [(0.0, 0.0), (30.0, 0.0), (30.0, 180.0), (45.0, 90.0), (45.0, 180.0)]
CASES = ["MEAS_10mA", "MEAS_5mA"]
CONFIRM = ROOT / "results/phase28/f10/ring_cell/confirm.jsonl"
REPORT = ROOT / "results/phase28/f10/ring_cell/report.json"


def _rows():
    import phase28_loop_measured as lm

    rows = []
    for name in ("measured_rows.jsonl", "measured_extra_rows.jsonl"):
        rows += lm.rows_of(ROOT / "results/phase28/f10/loop_R4" / name)
    return lm, rows


def tensors(case):
    """`case` = MEAS_10mA or MEAS_5mA: the measured MA4AGP907 (full Z(f)), mesh 90, 10 GHz, loop R4."""
    lm, rows = _rows()
    on_label = "ON_" + case.split("_")[1]
    idx = {lm.key(r): r for r in rows}
    out = {}
    for lab, flag in ((on_label, True), ("OFF_0V", False)):
        r = idx[lm.key({"state": lab, "mesh": 90, "f_hz": 10e9, "scale": None})]
        out[flag] = np.array([[complex(*c) for c in row] for row in r["y_s"]])
    return out[True], out[False]


def alpha0(case):
    lm, rows = _rows()
    return lm.report(rows)["cases"][case.split("_")[1]]["mesh90_10ghz"]["alpha_exact"]


def beta_of_b(b_ms):
    return find_beta_from_sheet_admittance(1j * b_ms * 1e-3, SUB, cf.F).beta_rad_per_m


def case_run(job):
    case, mode, beam = job
    cf.SUB = SUB
    cf.set_frequency(10.0e9)
    yon, yoff = tensors(case)
    ybar, dy = 0.5 * (yon + yoff), 0.5 * (yon - yoff)
    a0 = alpha0(case)
    poly, src = sf.source_xy()
    x, y = cf.lattice(poly, PITCH, keepout_m=0.0)
    rho = np.hypot(x - src[0], y - src[1])
    keep = rho >= sf.RHO_MIN
    x, y, rho = x[keep], y[keep], rho[keep]
    psi = np.arctan2(y - src[1], x - src[0])
    e = np.stack([np.cos(psi), np.sin(psi)], axis=1)

    def per_angle(p):
        ee = np.array([math.cos(p), math.sin(p)])
        if mode == "rotated":
            return beta_of_b(ybar[0, 0].imag * 1e3), 1.0, np.array([math.cos(p), math.sin(p)], dtype=complex)
        ratio = (ee @ dy @ ee) / dy[0, 0]
        return beta_of_b((ee @ ybar @ ee).imag * 1e3), abs(ratio) ** 2, (dy @ ee) / dy[0, 0]

    # per-cell response (the angle grid keeps the root solves cheap)
    grid = np.radians(np.linspace(-90, 90, 91))
    tab = [per_angle(p) for p in grid]
    beta_g = np.array([t[0] for t in tab])
    scale_g = np.array([t[1] for t in tab])
    vx = np.array([t[2][0] for t in tab])
    vy = np.array([t[2][1] for t in tab])
    beta = np.interp(psi, grid, beta_g)
    scale = np.interp(psi, grid, scale_g)
    v = np.stack([np.interp(psi, grid, vx.real) + 1j * np.interp(psi, grid, vx.imag), np.interp(psi, grid, vy.real) + 1j * np.interp(psi, grid, vy.imag)], axis=1)
    alpha_loc = a0 * scale
    f = np.maximum(np.cos(psi), 0.0) ** sf.P_FIELD
    amp = f * np.exp(-alpha_loc * rho) / np.sqrt(rho / sf.RHO_MIN)
    phi_ref = beta * rho
    th0, ph0 = math.radians(beam[0]), math.radians(beam[1])
    kt = cf.K0 * math.sin(th0)
    ktx, kty = kt * math.cos(ph0), kt * math.sin(ph0)
    m = np.cos(phi_ref - (ktx * x + kty * y))
    _, mq = cf.lp_states(m, 1)
    j = (amp * np.exp(-1j * phi_ref) * mq)[:, None] * v
    jx, jy = j[:, 0], j[:, 1]
    green = cf.green_function()
    r = cf.analyse(jx, jy, x, y, green, PITCH, "lp", beam)
    # converted fraction from alpha_loc along each ray
    psis = np.linspace(-math.pi / 2 + 1e-4, math.pi / 2 - 1e-4, 721)
    w = np.cos(psis) ** (2 * sf.P_FIELD)
    L = sf.chord_length(poly, src, psis)
    al = np.interp(psis, grid, a0 * scale_g)
    eta = float(np.sum(w * (1 - np.exp(-2 * al * L))) / np.sum(w))
    th, ph = np.array([r["peak_theta_deg"]]), np.array([r["peak_phi_deg"]])
    et, ep = cf.far_field(jx, jy, x, y, green, th, ph, PITCH)
    p = math.radians(ph[0])
    co, cr = et[0] * math.cos(p) - ep[0] * math.sin(p), et[0] * math.sin(p) + ep[0] * math.cos(p)
    xpd = cf.bounded_db(20 * math.log10(abs(co) / max(abs(cr), 1e-30)))
    # aliasing of the 1-bit sampling: modulation phase advance per lattice step along x and y
    kx_mod = beta * np.cos(psi) - ktx
    ky_mod = beta * np.sin(psi) - kty
    aliased = float(np.mean((np.abs(kx_mod) * PITCH > math.pi) | (np.abs(ky_mod) * PITCH > math.pi)))
    return {
        "case": case, "mode": mode, "beam": list(beam), "n_cells": int(x.size), "alpha0": a0, "eta_conv": eta,
        "d_ap_dbi": r["d_ap_dbi"], "gain_dbi": cf.launch_gain_dbi(r["d_ap_dbi"], eta), "peak_theta_deg": r["peak_theta_deg"], "peak_phi_deg": r["peak_phi_deg"],
        "pointing_error_deg": r["pointing_error_deg"], "sll_db": r["sll_db"], "xpd_at_peak": xpd, "aliased_cell_fraction": aliased,
        "beta_over_k0_range": [float(beta.min() / cf.K0), float(beta.max() / cf.K0)],
    }


def main(out, workers):
    jobs = [(c, m, b) for c in CASES for m in ("unrotated", "rotated") for b in BEAMS]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        rows = list(ex.map(case_run, jobs))
    Path(out).write_text(json.dumps({"labels": ["MODEL_ONLY", "IDEALISED_CELLS", "NOT_HFSS_VALIDATED"], "rows": rows}, indent=1), encoding="utf-8")
    print("wrote", out, len(rows))


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], int(a[a.index("--workers") + 1]) if "--workers" in a else 8)
