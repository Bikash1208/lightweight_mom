"""Phase 29: ghost-beam diagnostic (MODEL_ONLY, no synthesis): one feed at a time, strongest sidelobe (outside 1.5 HPBW) and whether it coincides with another feed's main beam; cell-family anisotropy.
Designs: 4-feed 300 mm (GA genes of tensor_ga_n4_300.json) and two-feed 150 mm (psi2 = 180 deg, the scan optimum). Writes results/phase29/ghost_diag.json."""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts/research"))
RES = ROOT / "results/phase29"


def analyse(d, green, grid=(1.0, 1.0)):
    from phase29_pointB_sweep import _hpbw

    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa

    th, ph = np.arange(0.0, 90.0 + 1e-9, grid[0]), np.arange(0.0, 360.0, grid[1])
    F = sa.port_fields(d, green, th, ph)
    G = [sa.gain_db(*f) for f in F]
    tg, _ = np.meshgrid(th, ph, indexing="ij")
    pk = []
    for g in G:
        lb = c.lobes(g, th, ph, within_db=3.0)
        pk.append((lb[0][1], lb[0][2]))
    rows = []
    for j, g in enumerate(G):
        hp = _hpbw(g, th, ph)
        t2, p2 = np.meshgrid(th, ph, indexing="ij")
        ang = np.degrees(np.arccos(np.clip(np.sum(c._unit(t2, p2) * c._unit(*pk[j]), axis=-1), -1, 1)))
        m = np.where(ang > 1.5 * hp, g, -1e9)
        it, ip = np.unravel_index(np.argmax(m), m.shape)
        sl = (float(th[it]), float(ph[ip]))
        rel = float(g[it, ip] - g.max())
        others = {
            str(k + 1): {
                "dist_deg": round(c.angle_between(sl, pk[k]), 1),
                "dist_hpbw": round(c.angle_between(sl, pk[k]) / hp, 1),
                "level_at_other_beam_db": round(float(g[int(round(pk[k][0] / grid[0])), int(round(pk[k][1] / grid[1])) % len(ph)] - g.max()), 1),
            }
            for k in range(len(G))
            if k != j
        }
        ghost = [k for k, v in others.items() if v["dist_deg"] <= 2 * hp]
        rows.append(
            {
                "port": j + 1,
                "peak_theta_phi": [round(pk[j][0], 1), round(pk[j][1], 1)],
                "hpbw": round(hp, 1),
                "sidelobe_theta_phi": [sl[0], sl[1]],
                "sidelobe_db": round(rel, 1),
                "others": others,
                "ghost_with_ports": ghost,
            }
        )
    return rows


def anisotropy(lut):
    out = {}
    for s in range(4):
        ey = []
        for ia in (0, 6):
            for ip in (0, 6):
                y = lut.Y[s, ia, ip]
                ys = 0.5 * (y + y.T)
                ev = np.linalg.eigvalsh(ys.imag)
                ey.append(round(float(abs(ev).max() / max(abs(ev).min(), 1e-12)), 2))
        y0 = lut.Y[s, 0, 0]
        out[f"slot{s}"] = {
            "Imag_Y_diag_a0_psi0_mS": [round(float(y0[0, 0].imag * 1e3), 3), round(float(y0[1, 1].imag * 1e3), 3)],
            "X_diag_ratio_a0_psi0": round(float(abs(y0[1, 1].imag / y0[0, 0].imag)), 2),
            "principal_Im_Y_ratio_(a0psi0,a0psi90,a90psi0,a90psi90)": ey,
            "dY_over_Ybar_M": round(float(lut.M_slot[s]), 3),
        }
    return out


def main():
    import phase29_tensor_ladder as L
    import phase29_weakmod as W

    from src.hologram import shared_aperture_tensor as st

    L.init()
    g4 = np.array(json.loads((RES / "tensor_ga_n4_300.json").read_text())["best_genes"])
    S = {"genes_n4": g4.tolist(), "anisotropy": anisotropy(L._S["lut"])}
    cfg = {"mode": "table", "M": None}
    for kind, g, key in (("n4", g4, "n4_300mm"), ("n2", [180.0], "n2_150mm")):
        d = W.design(kind, g, cfg, 240)
        S[key] = analyse(d, L._S["green"])
        print(key, json.dumps(S[key]))
    (RES / "ghost_diag.json").write_text(json.dumps(S, indent=1, default=float))
    print(json.dumps(S["anisotropy"]))


if __name__ == "__main__":
    main()
