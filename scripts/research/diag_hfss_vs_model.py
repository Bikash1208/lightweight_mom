"""DIAGNOSTIC (model only, no design change): HFSS full-board far field (BOARD_N240_RF_full_win.step) vs the board model.  MODEL_ONLY / NOT_HFSS_VALIDATED.
    SUB6_MESH=96 python diag_hfss_vs_model.py <ffd_off> <ffd_on> <out_dir>
HFSS .ffd: header '0 180 91' / '-180 180 181' / 'Frequencies 1' / 'Frequency f'; then theta-outer (0..180, 91) x phi-inner (-180..180, 181) rows of Re/Im Etheta, Re/Im Ephi (V, rE), e^{+jwt}.
RHCP = (Etheta + j Ephi)/sqrt2 (IEEE, e^{+jwt}) = the project's copol_component(.., 'rhcp') convention; LHCP = (Etheta - j Ephi)/sqrt2.
Model field: the same engine as the design (currents of sub6_board_model4.currents for the OFF / ON state, layered-medium spectral dyad), evaluated on the HFSS (theta, phi) grid for theta <= 90.
Variants: dihedral transform of the board (identity, rot 90/180/270, mirror about the lines at 0/45/90/135 deg), x conj (time-convention flip), x switched-cell OFF/ON response swapped.
"""
import json
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("SUB6_MESH", "96")
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np

import sub6_final_design as FD
import sub6_board_model4 as R4

SQ2 = math.sqrt(2.0)
TH_DEG = np.arange(0, 181, 2)   # the HFSS grid: theta 0 ... 180 step 2 deg (91 points), phi -180 ... 180 step 2 deg (181 points)
PH_DEG = np.arange(-180, 181, 2)


def load_ffd(path):
    a = np.loadtxt(path, skiprows=4)
    assert a.shape == (91 * 181, 4), a.shape
    et = (a[:, 0] + 1j * a[:, 1]).reshape(91, 181)
    ep = (a[:, 2] + 1j * a[:, 3]).reshape(91, 181)
    return et, ep


def model_field(S, jx, jy):
    th = np.radians(np.arange(0, 91, 2))
    phi = np.radians(PH_DEG)
    TH, PH = np.meshgrid(th, phi, indexing="ij")
    kx, ky = (R4.K0 * np.sin(TH) * np.cos(PH)).ravel(), (R4.K0 * np.sin(TH) * np.sin(PH)).ravel()
    g = S.green.spectral_dyad(kx, ky) * np.cos(TH).ravel()[:, None, None]
    Jx = np.zeros(kx.size, complex)
    Jy = np.zeros(kx.size, complex)
    for i0 in range(0, kx.size, 2000):
        sl = slice(i0, i0 + 2000)
        d = np.exp(1j * (np.outer(kx[sl], S.x) + np.outer(ky[sl], S.y))) * R4.M.PITCH**2
        Jx[sl], Jy[sl] = d @ jx, d @ jy
    ex = g[:, 0, 0] * Jx + g[:, 0, 1] * Jy
    ey = g[:, 1, 0] * Jx + g[:, 1, 1] * Jy
    c, s = np.cos(PH).ravel(), np.sin(PH).ravel()
    ct = np.cos(TH).ravel()
    et = (ex * ct * c + ey * ct * s).reshape(46, 181)
    ep = (-ex * s + ey * c).reshape(46, 181)
    return et, ep


def pad_to_180(e):
    out = np.zeros((91, 181), complex)
    out[:46] = e
    return out


def phi_index(phi_deg):
    """index of phi (deg, any real) on the grid -180 ... 180 step 2 (phi = 180 maps to the +180 end for exactly 180, else wrapped into [-180, 180))"""
    p = np.mod(np.asarray(phi_deg) + 180.0, 360.0) - 180.0
    idx = np.rint((p + 180.0) / 2.0).astype(int)
    return np.where((np.asarray(phi_deg) == 180.0), 180, idx)


def transform(et, ep, kind, arg):
    """fields of the board after a geometric transform of the board: kind 'rot' (CCW by arg deg) or 'mir' (mirror about the line through the origin at angle arg deg)."""
    ph = PH_DEG.astype(float)
    if kind == "rot":
        src = phi_index(ph - arg)
        return et[:, src], ep[:, src]
    src = phi_index(2 * arg - ph)
    return et[:, src], -ep[:, src]


def circ(et, ep):
    return (et + 1j * ep) / SQ2, (et - 1j * ep) / SQ2    # RHCP, LHCP (IEEE, e^{+jwt})


def prad(et, ep, tmax):
    """radiated 'power' integral of |Et|^2 + |Ep|^2 over theta in 0 ... tmax (deg, 2 deg grid), full phi circle (the duplicate +180 column is dropped)"""
    nt = tmax // 2 + 1
    u = (np.abs(et[:nt, :180]) ** 2 + np.abs(ep[:nt, :180]) ** 2)
    w = np.sin(np.radians(TH_DEG[:nt]))
    wt = np.full(nt, 1.0)
    wt[0] = wt[-1] = 0.5
    return float(np.sum(u * (w * wt)[:, None]) * math.radians(2.0) * math.radians(2.0))


def dir_db(e, ptot):
    return 10 * np.log10(np.maximum(4 * math.pi * np.abs(e) ** 2 / ptot, 1e-12))


def peaks(rh, tmin=6, tmax=60):
    """peak RHCP power on the phi ~ 0 side (|phi| < 90) and the phi ~ 180 side (|phi| >= 90) within theta tmin ... tmax (deg)"""
    i0, i1 = tmin // 2, tmax // 2
    p = np.abs(rh[i0:i1 + 1]) ** 2
    ph = PH_DEG
    out = {}
    for nm, sel in (("A_side", np.abs(ph) < 90), ("B_side", np.abs(ph) >= 90)):
        q = np.where(sel[None, :], p, 0.0)
        i = np.unravel_index(np.argmax(q), q.shape)
        out[nm] = (int(TH_DEG[i[0] + i0]), int(ph[i[1]]))
    return out


def at(db, th, phd):
    """value (dB array on the grid) at theta th (deg, linear interpolation of the power between the 2 deg rows) and phi phd (deg)"""
    i = int(th // 2)
    f = (th - 2 * i) / 2.0
    j = phi_index(phd)
    p = (1 - f) * 10 ** (db[i, j] / 10) + f * 10 ** (db[i + 1, j] / 10)
    return float(10 * np.log10(p))


def corr(a, b, tmax=60):
    n = tmax // 2 + 1
    x, y = a[:n, :180].ravel(), b[:n, :180].ravel()
    x, y = x - x.mean(), y - y.mean()
    return float(np.dot(x, y) / math.sqrt(np.dot(x, x) * np.dot(y, y) + 1e-300))


def ccorr(a, b, tmax=60):
    n = tmax // 2 + 1
    x, y = a[:n, :180].ravel(), b[:n, :180].ravel()
    return float(abs(np.vdot(x, y)) / math.sqrt(np.vdot(x, x).real * np.vdot(y, y).real + 1e-300))


def main():
    ffd_off, ffd_on, out_dir = sys.argv[1], sys.argv[2], Path(sys.argv[3])
    out_dir.mkdir(parents=True, exist_ok=True)
    H = {}
    for st, p in (("off", ffd_off), ("on", ffd_on)):
        et, ep = load_ffd(p)
        H[st] = (et, ep)
    fd = FD.final_design("cal")
    D, sa, ph = fd["D"], fd["sa"], fd["ph"]
    S = D["S"]
    fl = fd["flips"].astype(float)
    mask = D["mask"] & D["live"]
    cur = {}
    for st in ("off", "on"):
        jx, jy, a, r = R4.currents(D, sa, ph, fd["direction"], fl, st)
        cur[st] = (jx, jy)
    # model fields: [conjJ][mode][state]; mode 'none': own state; 'sw': switched cells take the OTHER state's currents (static cells keep their own)
    FM = {}
    JT = {"J-id": lambda jx, jy: (jx, jy), "J mirror(jy->-jy)": lambda jx, jy: (jx, -jy), "J mirror(jx->-jx)": lambda jx, jy: (-jx, jy), "J rot90": lambda jx, jy: (-jy, jx), "J swap(jx<->jy)": lambda jx, jy: (jy, jx)}
    for jt, fn in JT.items():
        for conjJ in (0, 1):
            for mode in ("none", "sw"):
                for st in ("off", "on"):
                    oth = {"off": "on", "on": "off"}[st]
                    jx = np.where(mask, cur[oth][0], cur[st][0]) if mode == "sw" else cur[st][0]
                    jy = np.where(mask, cur[oth][1], cur[st][1]) if mode == "sw" else cur[st][1]
                    if conjJ:
                        jx, jy = np.conj(jx), np.conj(jy)
                    jx, jy = fn(jx, jy)
                    et, ep = model_field(S, jx, jy)
                    FM[(jt, conjJ, mode, st)] = (pad_to_180(et), pad_to_180(ep))
    M = {st: FM[("J-id", 0, "none", st)] for st in ("off", "on")}
    # ---- directivity / beams, HFSS and model (identity)
    def summary(e2, tmax):
        et, ep = e2
        rh, lh = circ(et, ep)
        pt = prad(et, ep, tmax)
        dR, dL = dir_db(rh, pt), dir_db(lh, pt)
        pk = peaks(rh)
        res = dict(peaks=pk, D_RHCP_peak_dBi=float(dR[:31].max()), D_total_peak_dBi=float(dir_db(np.sqrt(np.abs(et) ** 2 + np.abs(ep) ** 2), pt)[:46].max()))
        for nm, (t, p) in pk.items():
            res[nm + "_D_RHCP_dBi"] = float(dR[t // 2, phi_index(p)])
            res[nm + "_D_LHCP_dBi"] = float(dL[t // 2, phi_index(p)])
        res["at_25_0"] = dict(RHCP=at(dR, 25, 0), LHCP=at(dL, 25, 0))
        res["at_25_180"] = dict(RHCP=at(dR, 25, 180), LHCP=at(dL, 25, 180))
        # fraction of the radiated power within 10 deg of the RHCP peak, and the half-power width of the strongest beam along phi and theta
        t0, p0 = max(pk.values(), key=lambda tp: np.abs(rh[tp[0] // 2, phi_index(tp[1])]))
        nt_ = tmax // 2 + 1
        TT, PP = np.meshgrid(np.radians(TH_DEG[:nt_]), np.radians(PH_DEG[:180]), indexing="ij")
        cd = np.cos(TT) * math.cos(math.radians(t0)) + np.sin(TT) * math.sin(math.radians(t0)) * np.cos(PP - math.radians(p0))
        near = np.degrees(np.arccos(np.clip(cd, -1, 1))) <= 10.0
        u = (np.abs(et[:nt_, :180]) ** 2 + np.abs(ep[:nt_, :180]) ** 2) * np.sin(np.radians(TH_DEG[:nt_]))[:, None]
        res["frac_power_within_10deg_of_main_beam"] = float(u[near].sum() / u.sum())
        res["main_beam"] = (int(t0), int(p0))
        cutp = np.abs(rh[:, phi_index(p0)]) ** 2
        res["hpbw_theta_deg"] = hpbw(cutp[:nt_], t0 // 2) * 2.0
        # SLL: strongest RHCP+LHCP power outside 10 deg of the main beam relative to the peak
        tot = np.abs(et[:nt_, :180]) ** 2 + np.abs(ep[:nt_, :180]) ** 2
        res["SLL_dB"] = float(10 * np.log10(tot[~near].max() / tot.max()))
        return res

    def hpbw(p, t0):
        pk = p[t0]
        lo = t0
        while lo > 0 and p[lo] > 0.5 * pk:
            lo -= 1
        hi = t0
        while hi < len(p) - 1 and p[hi] > 0.5 * pk:
            hi += 1
        return float(hi - lo)

    res = dict(label="MODEL_ONLY diagnostic", hfss={st: summary(H[st], 180) for st in H}, model_identity={st: summary(M[st], 90) for st in M})
    # ---- variants: geometry (8 dihedral) x conj of the radiated field (time-convention flip of the field) x conj of the cell excitation phase (sign of the cell phase / propagation convention)
    #               x switched-cell / state swap ('none', 'sw' = only the switched cells exchange OFF and ON, 'full' = both states exchanged)
    rows = []
    geoms = [("id", None, None), ("rot90", "rot", 90), ("rot180", "rot", 180), ("rot270", "rot", 270), ("mirror x->-x", "mir", 90), ("mirror y->-y", "mir", 0), ("swap x<->y", "mir", 45), ("swap x<->-y", "mir", 135)]
    hc = {st: circ(*H[st]) for st in ("off", "on")}
    hsum, hdel = hc["on"][0] + hc["off"][0], hc["on"][0] - hc["off"][0]
    for jt in JT:
      for (gname, kind, arg) in geoms:
        for conjE in (0, 1):
            for conjJ in (0, 1):
                for mode in ("none", "sw", "full"):
                    row = dict(variant=("" if jt == "J-id" else jt + " | ") + gname + (" + conj(E)" if conjE else "") + (" + conj(cell phase)" if conjJ else "") + {"none": "", "sw": " + switched-cell OFF/ON swapped", "full": " + both states swapped"}[mode],
                               geom=gname, jt=jt, conjE=conjE, conjJ=conjJ, swap=mode)
                    rhs = {}
                    for st in ("off", "on"):
                        src = {"off": "on", "on": "off"}[st] if mode == "full" else st
                        et, ep = FM[(jt, conjJ, "sw" if mode == "sw" else "none", src)]
                        if kind is not None:
                            et, ep = transform(et, ep, kind, arg)
                        if conjE:
                            et, ep = np.conj(et), np.conj(ep)
                        rh, lh = circ(et, ep)
                        rhs[st] = rh
                        hrh, hlh = hc[st]
                        row["corr_" + st] = corr(np.abs(rh) ** 2, np.abs(hrh) ** 2)
                        row["ccorr_" + st] = ccorr(rh, hrh)
                        row["corrL_" + st] = corr(np.abs(lh) ** 2, np.abs(hlh) ** 2)
                        row["pk_" + st] = peaks(rh)
                        pt = prad(et, ep, 90)
                        row["A25_RHCP_dBi_" + st] = at(dir_db(rh, pt), 25, 0)
                        row["B25_RHCP_dBi_" + st] = at(dir_db(rh, pt), 25, 180)
                    row["ccorr_sum"] = ccorr(rhs["on"] + rhs["off"], hsum)
                    row["ccorr_delta"] = ccorr(rhs["on"] - rhs["off"], hdel)
                    row["score"] = row["corr_off"] + row["corr_on"] + row["ccorr_sum"] + row["ccorr_delta"]
                    rows.append(row)
    rows.sort(key=lambda r: -r["score"])
    res["variants_sorted"] = rows
    (out_dir / "diag_results.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8", newline="\n")
    np.savez_compressed(out_dir / "fields.npz", hoff_t=H["off"][0], hoff_p=H["off"][1], hon_t=H["on"][0], hon_p=H["on"][1], moff_t=M["off"][0], moff_p=M["off"][1], mon_t=M["on"][0], mon_p=M["on"][1])
    print(json.dumps(dict(hfss=res["hfss"], model=res["model_identity"]), indent=1, default=str)[:3500])
    def pr(r):
        print("%-78s corr OFF %.3f ON %.3f | ccorr OFF %.3f ON %.3f sum %.3f delta %.3f | A@25/B@25 OFF %.1f/%.1f ON %.1f/%.1f | pk OFF %s ON %s" % (
            r["variant"], r["corr_off"], r["corr_on"], r["ccorr_off"], r["ccorr_on"], r["ccorr_sum"], r["ccorr_delta"], r["A25_RHCP_dBi_off"], r["B25_RHCP_dBi_off"], r["A25_RHCP_dBi_on"], r["B25_RHCP_dBi_on"], r["pk_off"], r["pk_on"]))
    print("--- top 12 by score")
    for r in rows[:12]:
        pr(r)
    print("--- requested single variants")
    want = ("id", "mirror x->-x", "mirror y->-y", "rot180", "id + conj(E)", "swap x<->y", "id + switched-cell OFF/ON swapped", "id + both states swapped", "id + conj(cell phase)")
    for r in rows:
        if r["variant"] in want:
            pr(r)

if __name__ == "__main__":
    main()
