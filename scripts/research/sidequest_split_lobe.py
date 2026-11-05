"""Where could the HFSS state-2 lobe at (20, 154) [-0.1 dB] come from?  Decompose the MODEL far field of the baseline board_K2_default (RT5880, lam 0.28, mu 30) into the static-cell part (same m in both states) and the switched part,
and test a state-dependent phase constant (beta) error.  MODEL_ONLY.  -> results/sidequest_hybrid/k2_split_lobe.json"""
import json
import math
import os
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
os.environ.setdefault("K4_A0", "0.844")
os.environ.setdefault("K4_BETA_K0", "1.0731")
import numpy as np
import sidequest_coded_k4 as K
import sidequest_k2_default_design as D

st = D.st


def pattern(bd, amp, m, state, mask=None):
    j, _ = bd.unit(amp)
    mm = m[state] if mask is None else m[state] * mask[None, :]
    cur = (mm[0] * j["x"][0] + mm[1] * j["y"][0], mm[0] * j["x"][1] + mm[1] * j["y"][1])
    et, ep = st.far_field(bd.x, bd.y, *cur)
    return np.abs(et) ** 2 + np.abs(ep) ** 2


def info(inten, ref):
    i = np.unravel_index(np.argmax(inten), inten.shape)
    return {"peak_dir_deg": [float(np.degrees(st.TG[i])), float(np.degrees(st.PG[i]))], "peak_rel_to_main_dB": float(10 * np.log10(inten[i] / ref + 1e-30))}


def level_at(inten, ref, th, ph):
    ii = np.argmin(np.abs(np.degrees(st.TG[:, 0]) - th))
    jj = np.argmin(np.abs(np.degrees(st.PG[0]) - ph))
    return float(10 * np.log10(inten[ii, jj] / ref + 1e-30))


if __name__ == "__main__" and not (len(sys.argv) > 1 and sys.argv[1] == "sens"):
    stat, _p, _n = D.library()
    levels = sorted({round(v[0], 4) for v in stat.values()})
    D.setup(10.5e9, 0.15, D.MARGIN)
    bd = K.KBoard([0, 180], ["LP-V", "LP-V"])
    r = D.run_design(bd, D.options(levels, [(1.0, -1.0, 1.0)], [(0, 1)]), [0.0, np.radians(180)], 0.28, 2.0, mu=30.0)
    m, amp, nd = r["m"], r["amp"], r["nd"]
    sw = nd > 0
    out = {"diodes": int(nd.sum()), "gain_model": [s["gain_dbi"] for s in r["states"]], "states": []}
    for k in range(2):
        full = pattern(bd, amp, m, k)
        ref = full.max()
        stt = pattern(bd, amp, m, k, mask=(~sw).astype(float))
        swo = pattern(bd, amp, m, k, mask=sw.astype(float))
        # strongest sidelobe of the full pattern outside 1.2 HPBW of the peak (as the metric), and its direction
        ipk = np.unravel_index(np.argmax(full), full.shape)
        cosd = np.cos(st.TG) * np.cos(st.TG[ipk]) + np.sin(st.TG) * np.sin(st.TG[ipk]) * np.cos(st.PG - st.PG[ipk])
        dist = np.degrees(np.arccos(np.clip(cosd, -1, 1)))
        d_dbi = 10 * math.log10(4 * math.pi * full[ipk] / float(np.sum(full * np.sin(st.TG) * math.radians(1) * math.radians(2))))
        hp = math.degrees(math.sqrt(16 / 10 ** (d_dbi / 10)))
        sl = np.where(dist > 1.2 * hp, full, 0)
        isl = np.unravel_index(np.argmax(sl), sl.shape)
        rec = {"state": k + 1, "main": info(full, ref), "strongest_sidelobe": {"dir_deg": [float(np.degrees(st.TG[isl])), float(np.degrees(st.PG[isl]))], "rel_dB": float(10 * np.log10(sl[isl] / ref))}, "HPBW_est_deg": hp,
               "static_cells_only": info(stt, ref), "switched_cells_only": info(swo, ref),
               "static_only_level_at_(25,0)": level_at(stt, ref, 25, 0), "static_only_level_at_(25,180)": level_at(stt, ref, 25, 180), "full_level_at_(20,154)": level_at(full, ref, 20, 154), "static_only_level_at_(20,154)": level_at(stt, ref, 20, 154)}
        out["states"].append(rec)
        print(rec, flush=True)
    # state-dependent beta: state 2 (diodes ON) has lower mean loading; beta/k0 shifts by dbeta in state 2 only
    base_beta = K.BETA
    res = []
    for dk in (-0.04, -0.02, 0.0, 0.02, 0.04, 0.08):
        K.BETA = base_beta + dk * st.K0
        full = pattern(bd, amp, m, 1)
        ref0 = None
        ipk = np.unravel_index(np.argmax(full), full.shape)
        res.append({"dbeta_over_k0": dk, "peak_dir_deg": [float(np.degrees(st.TG[ipk])), float(np.degrees(st.PG[ipk]))], "peak_intensity_rel": float(full.max()), "level_at_(20,154)_rel_to_own_peak_dB": level_at(full, full.max(), 20, 154)})
        print(res[-1], flush=True)
    K.BETA = base_beta
    out["state2_beta_shift"] = res
    Path("results/sidequest_hybrid/k2_split_lobe.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline="\n")


def sensitivity():
    """Isolation (other-beam level relative to the own peak) when the static levels are scaled by (1+d) or the switched cells' OFF level by (1+d): the wrong-beam suppression is an interference of the static and the switched fields."""
    stat, _p, _n = D.library()
    levels = sorted({round(v[0], 4) for v in stat.values()})
    D.setup(10.5e9, 0.15, D.MARGIN)
    bd = K.KBoard([0, 180], ["LP-V", "LP-V"])
    r = D.run_design(bd, D.options(levels, [(1.0, -1.0, 1.0)], [(0, 1)]), [0.0, np.radians(180)], 0.28, 2.0, mu=30.0)
    m, amp, nd = r["m"], r["amp"], r["nd"]
    sw = nd > 0
    rows = []
    for what in ("static", "switched_OFF", "switched_ON"):
        for d in (-0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15):
            mp = m.copy()
            if what == "static":
                mp[:, :, ~sw] *= 1 + d
            elif what == "switched_OFF":
                mp[0][:, sw] = m[0][:, sw] * (1 + d)
            else:
                mp[1][:, sw] = m[1][:, sw] * (1 + d)
            iso = []
            for k in range(2):
                full = pattern(bd, amp, mp, k)
                ref = full.max()
                iso.append(level_at(full, ref, 25, 180 if k == 0 else 0))
            rows.append({"perturbed": what, "d": d, "other_beam_rel_dB_state1_state2": iso})
            print(what, d, [round(x, 1) for x in iso], flush=True)
    return rows


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "sens":
    rows = sensitivity()
    Path("results/sidequest_hybrid/k2_split_lobe_sens.json").write_text(json.dumps(rows, indent=1), encoding="utf-8", newline=chr(10))
