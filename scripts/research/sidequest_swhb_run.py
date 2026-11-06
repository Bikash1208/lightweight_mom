"""Run candidate-B board designs.  MODEL_ONLY."""
import json
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_sll_def as SLD
import sidequest_swhb as S
import sidequest_swhb_design as DS

RES = ROOT / "results/sidequest_hybrid"


def patch_rows(path, mesh, psi, gap, lam="RO4003C", p=6.4):
    rows = json.loads(Path(path).read_text())
    return [r for r in rows if r["mesh"] == mesh and r["psi"] == psi and abs(r["gap"] - gap) < 1e-6 and r["lam"] == lam and abs(r["p"] - p) < 1e-6]


def metrics(des, outs, eta_l=1.0):
    th, ph = outs[0]["th"], outs[0]["ph"]
    thd, phd = np.degrees(th), np.degrees(ph)
    TH = np.meshgrid(th, ph, indexing="ij")[0]
    res = []
    for s, r in enumerate(outs):
        u = r["u"]  # W/sr per 1 W launched (r^2 |E|^2 / 2 eta)
        a = SLD.analyze(thd, phd, u, peak=None)
        if des["variant"] == "add" and s == 1:  # two wanted beams: SLL outside BOTH main lobes, relative to the stronger beam
            ins = np.zeros(u.shape, bool)
            u0 = 0.0
            for bm in des["beams"]:
                ii = np.unravel_index(np.argmin(np.abs(thd[:, None] - bm[0]) + np.abs(phd[None, :] - bm[1]) * 0.0 + 0.0), u.shape)
                mk, uu = SLD.analyze(thd, phd, u, peak=(float(bm[0]), float(bm[1])), return_mask=True)
                ins |= mk
                u0 = max(u0, uu)
            TT = np.meshgrid(thd, phd, indexing="ij")[0]
            out_ = np.where(ins | (TT > 90.0), -np.inf, u)
            a = dict(a)
            a["SLL_db"] = float(10 * math.log10(out_.max() / u0))
        want = des["beams"][s] if (des["variant"] == "swap" or s == 0) else des["beams"][0]
        i = np.unravel_index(np.argmax(u), u.shape)
        g = lambda t, p_: 10 * math.log10(4 * math.pi * u[np.argmin(abs(thd - t)), np.argmin(abs(phd - p_))] * eta_l + 1e-30)
        entry = {"state": s + 1, "beta_k0": r["beta_k0"], "alpha_loss_Np_m": float(r["alpha_loss"]), "peak_dir": a["peak_dir"], "G_peak_dBi": 10 * math.log10(4 * math.pi * u.max() * eta_l), "D_peak_dBi": 10 * math.log10(4 * math.pi * u.max() / r["P_rad_unit"]),
                 "P_rad": r["P_rad_unit"], "P_edge": r["P_edge"], "P_abs_other": 1 - r["P_rad_unit"] - r["P_edge"], "SLL_db": a["SLL_db"], "HPBW_deg": a["HPBW_mean_deg"], "kappa": r["kappa"],
                 "G_at_beam1_dBi": g(*des["beams"][0]), "G_at_beam2_dBi": g(*des["beams"][1]), "frac_60_90": float(np.sum(u * np.sin(TH) * (np.degrees(TH) >= 60)) / np.sum(u * np.sin(TH)))}
        # XPD at the peak: |E_phi|^2 vs |E_theta|^2 for y-pol beams
        eth, eph = r["eth"], r["eph"]
        ci = np.unravel_index(np.argmax(u), u.shape)
        t_, p_ = th[ci[0]], ph[ci[1]]
        ey = eth[ci] * math.cos(t_) * math.sin(p_) + eph[ci] * math.cos(p_)
        ex = eth[ci] * math.cos(t_) * math.cos(p_) - eph[ci] * math.sin(p_)
        entry["XPD_db"] = 20 * math.log10(abs(ey) / max(abs(ex), 1e-30))
        res.append(entry)
    return res


if __name__ == "__main__":
    rows = patch_rows(RES / "cellB_ypatch_final_m32_m64.json", 64, 0.0, 0.8)
    lib = DS.CellLib(rows, 3.38, 0.0027, 1.524e-3, 6.4e-3, L_ser=1.6e-9)
    for name, y0, y1, nd in lib.opts[:3]:
        print(name, y0 * 1e3, y1 * 1e3)
    t0 = time.time()
    des = DS.design(lib, 0.075, "swap", full=False, verbose=True)
    print("design time %.0f s n=%d nd=%d" % (time.time() - t0, des["n"], des["nd"]), flush=True)
    outs = DS.evaluate(des)
    for m in metrics(des, outs):
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in m.items()})
