"""SLL under the new definition: the highest peak OUTSIDE the main lobe, where the main lobe ends at its first null (first local minimum along each radial cut from the peak), or at 2 x HPBW if no null is found within that radius.
Also the 3 dB beamwidth (elevation cut and azimuth cut, angular).  Applied to (a) the MODEL of board v2 (r 2.30, state-dependent beta) and (b) the HFSS far field of board_K2_default (RT5880) for reference.
    -> results/sidequest_hybrid/sll_def_v2.json"""
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def analyze(th_deg, ph_deg, u, peak=None, r_max=70.0, dr=0.5):
    """th_deg (nt,), ph_deg (np,) with phi covering 360 (periodic), u (nt, np) intensity (linear). Returns peak, HPBW cuts, first-null radii, SLL (dB below the peak) and its direction."""
    ph = np.asarray(ph_deg, float)
    phe = np.concatenate([ph, [ph[0] + 360.0]]) if abs((ph[-1] - ph[0]) - 360.0) > 1e-6 else ph.copy()
    ue = np.concatenate([u, u[:, :1]], axis=1) if phe.size == ph.size + 1 else u
    f = RegularGridInterpolator((np.asarray(th_deg, float), phe), ue, bounds_error=False, fill_value=None)

    def uval(t, p):
        t = np.asarray(t, float)
        p = (np.asarray(p, float) - ph[0]) % 360.0 + ph[0]
        out = f(np.stack([np.clip(t, th_deg[0], th_deg[-1]), p], axis=-1))
        return np.where((t >= th_deg[0]) & (t <= th_deg[-1]), out, np.nan)

    i = np.unravel_index(np.nanargmax(u), u.shape)
    tp, pp = float(th_deg[i[0]]), float(ph[i[1]])
    if peak is not None:
        tp, pp = peak
    tr, pr = math.radians(tp), math.radians(pp)
    p0 = np.array([math.sin(tr) * math.cos(pr), math.sin(tr) * math.sin(pr), math.cos(tr)])
    e_th = np.array([math.cos(tr) * math.cos(pr), math.cos(tr) * math.sin(pr), -math.sin(tr)])
    e_ph = np.array([-math.sin(pr), math.cos(pr), 0.0])
    u0 = float(np.asarray(uval(tp, pp)).reshape(-1)[0])
    alphas = np.arange(0, 360, 10.0)
    deltas = np.arange(0.0, r_max + dr, dr)
    r_null, r_half = {}, {}
    cuts = {}
    for a in alphas:
        ar = math.radians(a)
        d = np.radians(deltas)[:, None]
        pts = np.cos(d) * p0 + np.sin(d) * (math.cos(ar) * e_th + math.sin(ar) * e_ph)
        t = np.degrees(np.arccos(np.clip(pts[:, 2], -1, 1)))
        p = np.degrees(np.arctan2(pts[:, 1], pts[:, 0]))
        v = uval(t, p) / u0
        cuts[a] = v
        half = np.where(v < 0.5)[0]
        r_half[a] = deltas[half[0]] if half.size else np.nan
        # first null: first local minimum after the main lobe has fallen below 0.5 and then rises for 2 consecutive samples
        null = np.nan
        start = half[0] if half.size else 0
        for k in range(start + 1, len(v) - 2):
            if np.isnan(v[k + 2]):
                break
            if v[k + 1] > v[k] * 1.0000001 and v[k + 2] > v[k + 1] * 1.0000001 and v[k] <= v[k - 1]:
                null = deltas[k]
                break
        r_null[a] = null
    el = r_half[0.0] + r_half[180.0] if not (np.isnan(r_half[0.0]) or np.isnan(r_half[180.0])) else float("nan")
    az = r_half[90.0] + r_half[270.0] if not (np.isnan(r_half[90.0]) or np.isnan(r_half[270.0])) else float("nan")
    hp_mean = np.nanmean([2 * r_half[a] for a in alphas])
    # radius of the main lobe per azimuth angle alpha: first null, else 2 x HPBW
    cap = 2.0 * hp_mean if not np.isnan(hp_mean) else 40.0
    rr = np.array([r_null[a] if not np.isnan(r_null[a]) else cap for a in alphas])
    TT, PP = np.meshgrid(np.asarray(th_deg, float), ph, indexing="ij")
    tt, pp_ = np.radians(TT), np.radians(PP)
    pts = np.stack([np.sin(tt) * np.cos(pp_), np.sin(tt) * np.sin(pp_), np.cos(tt)], axis=-1)
    cosd = np.clip(pts @ p0, -1, 1)
    dist = np.degrees(np.arccos(cosd))
    ang = (np.degrees(np.arctan2(pts @ e_ph, pts @ e_th)) % 360.0)
    r_dir = np.interp(ang, np.append(alphas, 360.0), np.append(rr, rr[0]))
    inside = dist <= r_dir
    outside = np.where(inside | (TT > 90.0), -np.inf, u)
    j = np.unravel_index(np.argmax(outside), outside.shape)
    sll = 10 * math.log10(outside[j] / u0)
    return {"peak_dir": [tp, pp], "HPBW_elevation_deg": el, "HPBW_azimuth_angular_deg": az, "HPBW_mean_deg": float(hp_mean), "first_null_radius_deg": {int(a): (None if np.isnan(r_null[a]) else float(r_null[a])) for a in alphas[::3]},
            "n_dirs_without_null": int(sum(np.isnan(r_null[a]) for a in alphas)), "SLL_db": float(sll), "SLL_dir": [float(TT[j]), float(PP[j])], "mainlobe_radius_mean_deg": float(rr.mean())}


def model_v2():
    import sidequest_ro3006_design as RD
    import sidequest_ro3006_v2 as V
    import sidequest_a0_direct as A

    cells = json.loads(Path("results/sidequest_hybrid/ro3006_board_design_r2.3_v2.json").read_text())["cells"]
    lib = RD.load_lib(96, True)
    betas = (1.33, 1.26)
    a0s = [A.a0_direct(A.ybar_for_beta(b, RD.ER, RD.H * 1e-3, RD.F), lib["dY_cal_mS"] * 1e-3, RD.ER, RD.H * 1e-3, RD.F, tand=RD.TD)[0] for b in betas]
    lib["A0"] = float(sum(a0s) / 2)
    np_, K, D, bd = RD.engine(lib, 0.71)
    st = D.st
    lv, pairs = lib["levels"], lib["pairs"]
    n = bd.n
    m = np.zeros((2, 2, n))
    nd = np.zeros(n)
    for i, (ix, iy, kd, tk, fr) in enumerate(cells):
        if kd == "SW":
            m[0, :, i], m[1, :, i] = pairs[fr][0], pairs[fr][1]
            nd[i] = 2
        else:
            v = 1.0 if (tk == "open" and fr == 1.0) else (lv["solid"] if tk == "none" else (lv["strap_%s" % ("0.79" if fr < 1 else "1.0")] if tk == "strap" else lv["open_%s" % fr]))
            m[:, :, i] = v
    sw = nd > 0
    loss = np.zeros((2, n))
    loss[1, sw] = 1.0
    amp = np.sqrt(2 * K.A0 / bd.rho) * np.exp(-0.5 * (K.A0 + K.AA0) * (bd.rho - st.RHO_KEEP))
    for _ in range(4):
        amp, conv, edge = K.amp_profile(bd, m, loss)
    out = {"edge": float(edge), "conv": float(conv), "states": []}
    th = np.degrees(st.TG[:, 0])
    ph = np.degrees(st.PG[0])
    for k in range(2):
        K.BETA = betas[k] * st.K0
        j, _ = bd.unit(amp)
        cur = (m[k, 0] * j["x"][0] + m[k, 1] * j["y"][0], m[k, 0] * j["x"][1] + m[k, 1] * j["y"][1])
        et, ep = st.far_field(bd.x, bd.y, *cur)
        u = np.abs(et) ** 2 + np.abs(ep) ** 2
        a = analyze(th, ph, u)
        out["states"].append(a)
        print("model v2 state", k + 1, {x: (round(y, 1) if isinstance(y, float) else y) for x, y in a.items() if x in ("peak_dir", "HPBW_elevation_deg", "HPBW_azimuth_angular_deg", "HPBW_mean_deg", "SLL_db", "SLL_dir", "n_dirs_without_null", "mainlobe_radius_mean_deg")}, flush=True)
    return out


def hfss_ref():
    import sidequest_ffd_analysis as FA

    out = {}
    th = np.arange(0, 181, 2.0)
    ph = np.arange(-180, 180, 2.0)
    for k, p in FA.FILES.items():
        et, ep = FA.load(p)
        u = (np.abs(et) ** 2 + np.abs(ep) ** 2)[:, :-1]
        sel = th <= 90.0
        a = analyze(th[sel], ph, u[sel])
        out[k] = a
        print("HFSS board_K2_default state", k, {x: (round(y, 1) if isinstance(y, float) else y) for x, y in a.items() if x in ("peak_dir", "HPBW_elevation_deg", "HPBW_azimuth_angular_deg", "HPBW_mean_deg", "SLL_db", "SLL_dir", "n_dirs_without_null", "mainlobe_radius_mean_deg")}, flush=True)
    return out


if __name__ == "__main__":
    res = {"model_v2": model_v2(), "hfss_board_K2_default": hfss_ref()}
    Path("results/sidequest_hybrid/sll_def_v2.json").write_text(json.dumps(res, indent=1), encoding="utf-8", newline=chr(10))
