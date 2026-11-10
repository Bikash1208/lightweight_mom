"""PARALLEL TEST items 1 + 2 (model only, no design change): full Bloch sheet tensors (diag_full_tensor_run.py -> tensors.json) -> report (item 1) and the board model for the EXACT STEP cell set with the
full tensors instead of the reflection law (item 2), vs the HFSS OFF / ON far field.   MODEL_ONLY / NOT_HFSS_VALIDATED.
    SUB6_MESH=96 python diag_full_tensor_board.py <tensors.json> <fields.npz (HFSS H fields)> <out.json>
Cell-frame tensor Y_c(theta') (u = slot axis, v = across; kt = beta k0 (cos theta', sin theta')) is interpolated over theta' by the exact symmetric series (even: Y_uu, Y_vv; odd: Y_uv; period 180).
Site current:  J = E0 exp(-jP) R(axis) Y_c(theta') (cos theta', sin theta')^T,  theta' = phi_R - axis,  E0 = a / (|Im Delta|/2),  a = sqrt(2 alpha W / rho) (the board model's own amplitude),
Delta = Y_uu(0) - Y_vv(90) of that cell and state.  The same alpha, W, P (phase), axes and flips as the reflection-law model; only the direction/sign law is replaced by the measured tensor
(the island ON/OFF ratio then comes out of the data).
"""
import json
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("SUB6_MESH", "96")
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np

import diag_hfss_vs_model as G
import sub6_final_design as FD
import sub6_board_model4 as R4

ANG = np.array([0.0, 22.5, 45.0, 67.5, 90.0])


def load_tensors(path):
    d = json.load(open(path))
    out = {}
    for k, v in d.items():
        row, st, th = k.split("|")
        out.setdefault((row, st), {})[float(th)] = np.array(v["Y_re"]) + 1j * np.array(v["Y_im"])
    return out


class Tensor:
    """Y_c(theta') of one cell/state from the 5 nodes 0 ... 90 deg"""

    def __init__(self, nodes, mode="full"):
        self.nodes = nodes
        self.mode = mode
        a = np.array([nodes[t] for t in ANG])                       # (5, 2, 2)
        m = np.arange(5)
        cosb = np.cos(2 * np.radians(ANG)[:, None] * m[None, :])    # (5, 5)
        self.cu = np.linalg.solve(cosb, a[:, 0, 0]); self.cv = np.linalg.solve(cosb, a[:, 1, 1])
        sinb = np.sin(2 * np.radians(ANG[1:4])[:, None] * np.arange(1, 4)[None, :])
        self.su = np.linalg.solve(sinb, 0.5 * (a[1:4, 0, 1] + a[1:4, 1, 0]))
        self.Yu, self.Yv = a[0, 0, 0], a[4, 1, 1]                   # kt along the slot (TM along u) / across (TM along v)
        self.Ybar = 0.5 * (self.Yu + self.Yv)
        self.Delta = self.Yu - self.Yv

    def with_mode(self, mode):
        t = Tensor.__new__(Tensor)
        t.__dict__.update(self.__dict__)
        t.mode = mode
        return t

    def Y(self, th):
        th = np.asarray(th, float)
        t2 = 2 * th
        out = np.zeros(th.shape + (2, 2), complex)
        if self.mode == "ideal":                                     # tensor sheet: diag(Yu, Yv) independent of the kt direction
            out[..., 0, 0], out[..., 1, 1] = self.Yu, self.Yv
            return out
        out[..., 0, 0] = sum(self.cu[m] * np.cos(m * t2) for m in range(5))
        out[..., 1, 1] = sum(self.cv[m] * np.cos(m * t2) for m in range(5))
        if self.mode != "nooff":
            off = sum(self.su[m - 1] * np.sin(m * t2) for m in range(1, 4))
            out[..., 0, 1] = out[..., 1, 0] = off
        return out


def _reflect_check(t, th):
    r = math.radians(th)
    e = np.array([math.cos(r), math.sin(r)]); m = np.array([math.cos(r), -math.sin(r)])
    j = t.Y(np.array(r)) @ e
    jm = j - t.Ybar * e
    law = 0.5 * t.Delta * m
    c = np.vdot(law, jm) / np.vdot(law, law)
    return dict(abs_ratio=float(abs(c)), phase_deg=float(np.degrees(np.angle(c))), residual_rel=float(np.linalg.norm(jm - c * law) / np.linalg.norm(law)))


def report_item1(T):
    rep = {}
    for key in sorted(T):
        nodes = T[key]
        if any(t not in nodes for t in ANG):
            continue
        t = Tensor(nodes)
        n = nodes
        ya, yc = t.Yu, t.Yv
        rep["|".join(key)] = dict(
            Y_along_mS=[ya.real * 1e3, ya.imag * 1e3], Y_across_mS=[yc.real * 1e3, yc.imag * 1e3],
            ratio_across_over_minus_along_abs=float(abs(yc) / abs(ya)), phase_diff_deg=float(np.degrees(np.angle(yc / (-ya)))),
            mean_mS=[t.Ybar.real * 1e3, t.Ybar.imag * 1e3], Delta_mS=[t.Delta.real * 1e3, t.Delta.imag * 1e3],
            Y_vv_with_kt_along_u_mS=[n[0.0][1, 1].real * 1e3, n[0.0][1, 1].imag * 1e3], Y_uu_with_kt_along_v_mS=[n[90.0][0, 0].real * 1e3, n[90.0][0, 0].imag * 1e3],
            tensor_sheet_violation_Yuu=float(abs(n[90.0][0, 0] - ya) / abs(ya)), tensor_sheet_violation_Yvv=float(abs(n[0.0][1, 1] - yc) / abs(yc)),
            Y_off_mS_by_theta={str(th): [n[th][0, 1].real * 1e3, n[th][0, 1].imag * 1e3] for th in ANG},
            Y_off_over_Delta_by_theta={str(th): float(abs(n[th][0, 1]) / abs(t.Delta)) for th in ANG},
            Y_uu_by_theta_mS={str(th): [n[th][0, 0].real * 1e3, n[th][0, 0].imag * 1e3] for th in ANG},
            Y_vv_by_theta_mS={str(th): [n[th][1, 1].real * 1e3, n[th][1, 1].imag * 1e3] for th in ANG},
            reflect_check={str(th): _reflect_check(t, th) for th in (22.5, 45.0, 67.5)})
    return rep


def site_current(D, sa, ph, fd, st, cl, mode_static, mode_island, subtract_mean, S, mask, live):
    """full-tensor site currents of state st ('off'/'on')"""
    phiR = np.arctan2(S.y, S.x)
    axis = fd["axis"]
    th = phiR - axis
    _, _, a, _ = R4.currents(D, sa, ph, fd["direction"], fd["flips"].astype(float), st)
    P = ph["Poff"] if st == "off" else ph["Pon"]
    jx = np.zeros(S.n, complex); jy = np.zeros(S.n, complex)
    ca, sn = np.cos(axis), np.sin(axis)
    for i in range(S.n):
        if not live[i]:
            continue
        cls = cl[(i, st)]
        t = cls.with_mode(mode_island if mask[i] else mode_static)
        e0 = a[i] / (0.5 * abs(cls.Delta.imag))
        e = np.array([math.cos(th[i]), math.sin(th[i])])
        jc = t.Y(np.array(th[i])) @ e
        if subtract_mean:
            jc = jc - cls.Ybar * e
        jc = e0 * np.exp(-1j * P[i]) * jc
        jx[i] = ca[i] * jc[0] - sn[i] * jc[1]
        jy[i] = sn[i] * jc[0] + ca[i] * jc[1]
    return jx, jy


def main():
    tpath, npz, outp = sys.argv[1:4]
    T = load_tensors(tpath)
    Hz = np.load(npz)
    H = {"off": (Hz["hoff_t"], Hz["hoff_p"]), "on": (Hz["hon_t"], Hz["hon_p"])}
    fd = FD.final_design("cal")
    D, sa, ph = fd["D"], fd["sa"], fd["ph"]
    S = D["S"]
    mask, live = D["mask"], D["live"]
    ist = D["ist"]
    ids = D["T"]["st"]["id"]
    isl_id = D["T"]["sw"]["id"][0]
    base = {}
    cl = {}
    for i in range(S.n):
        if not live[i]:
            continue
        for st in ("off", "on"):
            key = (isl_id, "OFF_BAR50" if st == "off" else "ON_BAR50") if mask[i] else (ids[ist[i]], "NOCUT")
            if key not in base:
                base[key] = Tensor(T[key])
            cl[(i, st)] = base[key]
    res = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED full-tensor test")
    res["item1"] = report_item1(T)
    kk = (isl_id, "OFF_BAR50")
    if -45.0 in T[kk]:
        y45, ym45 = T[kk][45.0], T[kk][-45.0]
        res["island_symmetry_check_OFF"] = dict(Y_uu_pm45_mS=[y45[0, 0] * 1e3, ym45[0, 0] * 1e3], Y_vv_pm45_mS=[y45[1, 1] * 1e3, ym45[1, 1] * 1e3], Y_uv_pm45_mS=[y45[0, 1] * 1e3, ym45[0, 1] * 1e3])

    def hfss_ref(st):
        hrh = G.circ(*H[st])[0]
        hpt = G.prad(H[st][0], H[st][1], 180)
        hd = G.dir_db(hrh, hpt)
        return dict(A25_dBi=G.at(hd, 25, 0), B25_dBi=G.at(hd, 25, 180), Dmax_RHCP=float(hd[:31].max()), peaks=G.peaks(hrh))
    res["hfss"] = {st: hfss_ref(st) for st in ("off", "on")}
    he = {st: G.circ(*H[st])[0][12, 180] for st in ("off", "on")}
    res["hfss_Eon_over_Eoff_24_180"] = dict(abs=float(abs(he["on"] / he["off"])), dB=float(20 * math.log10(abs(he["on"] / he["off"]))), phase_deg=float(np.degrees(np.angle(he["on"] / he["off"]))))
    fl = fd["flips"].astype(float)
    variants = {"V0 reflection law (existing model)": {st: R4.currents(D, sa, ph, fd["direction"], fl, st)[:2] for st in ("off", "on")}}
    spec = [("V1 ideal tensor sheet diag(Yu,Yv), ALL current", "ideal", "ideal", False),
            ("V2 ideal tensor sheet, modulated part only", "ideal", "ideal", True),
            ("V3 FULL measured tensor, ALL current", "full", "full", False),
            ("V4 FULL measured tensor, modulated part only", "full", "full", True),
            ("V5 static FULL / island ideal, mod only", "full", "ideal", True),
            ("V6 static ideal / island FULL, mod only", "ideal", "full", True),
            ("V7 FULL without Y_uv, mod only", "nooff", "nooff", True),
            ("V8 static FULL / island ideal, ALL current", "full", "ideal", False),
            ("V9 static ideal / island FULL, ALL current", "ideal", "full", False)]
    for name, ms, mi, sub in spec:
        variants[name] = {st: site_current(D, sa, ph, fd, st, cl, ms, mi, sub, S, mask, live) for st in ("off", "on")}
        print("currents", name, flush=True)
    res["variants"] = {}
    fields = {}
    for name, cur in variants.items():
        row = {}
        for st in ("off", "on"):
            et, ep = G.model_field(S, *cur[st])
            e2 = (G.pad_to_180(et), G.pad_to_180(ep))
            rh = G.circ(*e2)[0]
            pt = G.prad(e2[0], e2[1], 90)
            dR = G.dir_db(rh, pt)
            hrh = G.circ(*H[st])[0]
            row[st] = dict(A25_dBi=G.at(dR, 25, 0), B25_dBi=G.at(dR, 25, 180), Dmax_RHCP=float(dR[:31].max()), peaks=G.peaks(rh), corr_pow=G.corr(np.abs(rh) ** 2, np.abs(hrh) ** 2),
                           ccorr=G.ccorr(rh, hrh), E24_180=rh[12, 180], mean_current_over_modulated=None)
            fields[(name, st)] = rh
        eo, ef = row["on"]["E24_180"], row["off"]["E24_180"]
        row["Eon_over_Eoff_24_180"] = dict(abs=float(abs(eo / ef)), dB=float(20 * math.log10(max(abs(eo / ef), 1e-9))), phase_deg=float(np.degrees(np.angle(eo / ef))))
        res["variants"][name] = row
        print(name, {st: (round(row[st]["A25_dBi"], 1), round(row[st]["B25_dBi"], 1), round(row[st]["ccorr"], 3), round(row[st]["corr_pow"], 3)) for st in ("off", "on")}, "Eon/Eoff dB", round(row["Eon_over_Eoff_24_180"]["dB"], 1), flush=True)
    Path(outp).write_text(json.dumps(res, indent=1, default=lambda o: [o.real, o.imag] if isinstance(o, complex) else str(o)), encoding="utf-8", newline="\n")
    print("HFSS", json.dumps({st: {k: v for k, v in res["hfss"][st].items() if k in ("A25_dBi", "B25_dBi")} for st in res["hfss"]}), res["hfss_Eon_over_Eoff_24_180"])


if __name__ == "__main__":
    main()
