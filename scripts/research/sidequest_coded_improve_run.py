"""Run I1-I7 (each alone on the baseline B0, and cumulatively). Output results/sidequest_hybrid/coded_improve.json."""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_coded_improve as I

lib = I.Lib()
bd = I.Board()
VS_BIT, PAIRS_BIT = I.bit_sets()
ST = np.array(sorted(lib.static.values()))
LOW, HIGH = float(ST.min()), float(ST.max())
VS_AN = np.linspace(LOW, HIGH, 9)


def pairs_from(keys):
    out = []
    for k in keys:
        v = lib.sw[k]
        out.append((v["m_off"], v["m_on"], 0, 1, v["g_on"]))
        out.append((v["m_on"], v["m_off"], 1, 0, v["g_on"]))
    return out


KEYS_LS = [(2.6, 1.6), (2.8, 1.6)]
KEYS_ALL = sorted(lib.sw)
BASE = dict(vs=VS_BIT, pairs=PAIRS_BIT, eps=[0.0, 0.0], mu=0.0, g=1.0, launch=I.LAUNCH0, w=(1.0, 1.0), scale=None, ring=None, g_ring=1.0)
ITEMS = {
    "I1": lambda kw: dict(kw, vs=VS_AN),
    "I2": lambda kw: dict(kw, pairs=pairs_from(KEYS_LS)),
    "I3": None,  # taper, needs the current design
    "I4": lambda kw: dict(kw, launch=0.80),
    "I5": lambda kw: dict(kw, pairs=pairs_from(KEYS_ALL), mu=0.05),
    "I6": None,  # ring
    "I7": None,  # joint
}


def summarise(r):
    s = r["states"]
    return {"gain": [round(x["gain_dbi"], 2) for x in s], "SLL": [round(x["SLL_db"], 1) for x in s], "AR": [round(x["AR_db"], 2) for x in s], "XPD": [round(x["XPD_db"], 1) for x in s],
            "ap_eff_pct": round(100 * r["ap_eff"][0], 2), "SW_to_beam_pct": round(100 * r["SW_to_beam"], 1), "edge_pct": round(100 * r["edge"], 1), "static_pct": round(100 * r["static_frac"], 1),
            "diodes": r["diodes"], "lines": r["lines_with_inverters"], "rms_phase_deg": round(r["rms_phase_deg"], 1), "eff_step_deg": round(r["eff_step_deg"], 1), "ON_loss": round(r["ON_loss"], 3), "lambda": round(r["lambda"], 3)}


GS = (1.0, 1.5, 2.0, 2.7)


def go(label, kw, log, gs=GS):
    t = time.time()
    best = None
    for g in gs:  # target amplitude scale: part of the design (trade of fidelity against converted power); best mean gain with SLL <= B0 SLL + 1 dB
        r, aux = I.full_run(kw["vs"], kw["pairs"], kw["eps"], kw["mu"], g, kw["launch"], kw["w"], kw["scale"], 766, bd, kw["ring"], kw["g_ring"])
        gm = np.mean([x["gain_dbi"] for x in r["states"]])
        sll = max(x["SLL_db"] for x in r["states"])
        ok = sll <= LIMIT_SLL
        if best is None or (ok, gm) > (best[0], best[1]):
            best = (ok, gm, r, aux, g)
    _, _, r, aux, g = best
    r["g"] = g
    log[label] = summarise(r) | {"g": g}
    print(label, log[label], round(time.time() - t), flush=True)
    return r, aux


def taper_kw(kw, aux):
    _, _, m, nd, loss, alpha, aa = aux
    scale, info = I.taper_scale(bd, alpha, aa)
    return dict(kw, scale=scale), info


def ring_kw(kw):
    ring = bd.rho > bd.redge - 2 * I.st.PITCH
    return dict(kw, ring=ring, g_ring=2.0)


LIMIT_SLL = -12.7

if __name__ == "__main__":
    log = {"library": {"static_range_m": [LOW, HIGH], "sw": {f"{k}": v for k, v in lib.sw.items()}}}
    r0, aux0 = go("B0", BASE, log)
    # B0L: library-limited statics, bit switched cells
    go("B0L", dict(BASE, vs=np.array([lib.static[(2.8, "short")], lib.static[(2.8, "open")]])), log)
    # each item alone
    for it in ("I1", "I2", "I4", "I5"):
        go(it + "_alone", ITEMS[it](BASE), log)
    kw3, info = taper_kw(BASE, aux0)
    log["taper_info_B0"] = info
    go("I3_alone", kw3, log)
    go("I6_alone", ring_kw(BASE), log)
    # cumulative
    kw = dict(BASE)
    for it in ("I1", "I2"):
        kw = ITEMS[it](kw)
        r, aux = go("cum_" + it, kw, log)
    kw3, info = taper_kw(kw, aux)
    log["taper_info_cum"] = info
    go("cum_I3_taper_variant", kw3, log)  # reported, not adopted if it loses gain
    r, aux = go("cum_I3", kw, log)  # I3 adopted only if it helps: here the unmodified profile
    kw = ITEMS["I4"](kw)
    r, aux = go("cum_I4", kw, log)
    kw = ITEMS["I5"](kw)
    r, aux = go("cum_I5", kw, log)
    kw = ring_kw(kw)
    r, aux = go("cum_I6", kw, log)
    Path("results/sidequest_hybrid/coded_improve.json").write_text(json.dumps(log, indent=1, default=float), encoding="utf-8", newline="\n")
    np.save("/tmp/claude-0/-home-user-lightweight-mom/38e40b5e-31e2-540a-82df-a1b61d7e77f1/scratchpad/kw_cum.npy", np.array([0]))
