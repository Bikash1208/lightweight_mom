"""CP option (i): fixed RHCP beam 1 at (25, 0) from Phase 26 rotating-slot tensor cells + LP beam 2 at (25, 180) from the 180 GB crossed-slit cells (2 diodes each, switched on one line).  MODEL_ONLY.
Engine: Phase 26 ray-march / cell_design / vector far field (src/hologram/tensor_*), gain = D_ap x launch x eta_conv x anchor (anchor = HFSS 23.03 - model 23.51 dBi of the single hexagon).
GB cells: sheet-current element J = +-amp exp(-j phi_ref) rho_hat with amp = sqrt(q_GB), q_GB the Born radiated power per area of the cell's susceptance step against the local mean sheet (same units as the tensor 'amp' = sqrt(q)).
    python sidequest_cp_hybrid.py"""
import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_a0_direct as A0
import sidequest_sll_def as SLD
import sidequest_swhb as S
from src.green_function import LayeredGreenFunction
from src.hologram import tensor_cp_design as T
from src.hologram.leaky_forward import pattern_budget
from src.hologram.tensor_leaky_cp import SupercellTable, cell_design, march_rays
from src.pattern.vector_aperture_pattern import circular_components, copol_component

PKG = ROOT / "hfss_exports/sidequest_board_k2_lp"
RES = ROOT / "results/sidequest_hybrid"
P = 6.4e-3
F = T.F
ANCHOR_DB = 23.03 - 23.511963450519115  # HFSS - model of the single-hexagon RHCP board (Phase 26 C)
ENV = S.Env(2.2, 0.0009, 1.575e-3, F)
BEAM1 = (25.0, 0.0)
BEAM2 = (25.0, 180.0)
SLL_OFFSET_DB = 2.5  # the A/B calibration: model SLL is 2.5 dB optimistic


def load_cells():
    rows = list(csv.DictReader(open(PKG / "cells.csv")))
    x = np.array([float(r["x_mm"]) * 1e-3 for r in rows])
    y = np.array([float(r["y_mm"]) * 1e-3 for r in rows])
    kind = np.array([r["kind"] for r in rows])
    grp = np.array([r["group"] for r in rows])
    return x, y, kind, grp


def lib_y():
    lib = json.loads((RES / "k2_default_lib.json").read_text())["rows"]
    g = lambda m: complex(*next(r for r in lib if r["family"] == "A" and abs(r["ls"] - 2.8) < 1e-9 and r["mode"] == m)["Y_mS"][0]) * 1e-3
    return {"open": g("open"), "ON": g("ON")}


def tensor_design(x, y, theta0=25.0, phi0=0.0, eta_design=0.9, taper=0.5, hand=+1, n_psi=720, table=None):
    table = table or SupercellTable.from_files(T.TABLE_GLOB)
    th0, ph0 = np.radians(theta0), np.radians(phi0)
    m = march_rays(table, T.hexagon(), (0.0, 0.0), T.K0, th0, ph0, eta_design, taper_c=taper, n_psi=n_psi)
    d = cell_design(m, table, x, y, (0.0, 0.0), T.K0, th0, ph0, hand=hand)
    return table, m, d


def gb_currents(d, gb_mask, state_on, rel_sign=+1.0):
    """J (jx, jy) of the GB cells (in the amp = sqrt(q) units of the tensor design) for GB OFF / ON."""
    Yd = lib_y()
    y_state = Yd["ON"] if state_on else Yd["open"]
    rho, phi = d["rho"][gb_mask], d["phi"][gb_mask]
    beta = d["beta_over_k0"][gb_mask]
    w = d["w"][gb_mask]
    phi_ref = d["phi_ref"][gb_mask]
    ybar = np.array([A0.ybar_for_beta(max(b, 1.001), ENV.er, ENV.h, F) for b in beta])  # local mean-sheet susceptance [S]
    dy = (y_state.imag - ybar) * 1.0
    s = np.array([ENV.tm_power(max(b, 1.001)) for b in beta])
    q = np.array([S.element_q(ENV, 1j * yb, 1j * yb) for yb in ybar])  # W per (A m)^2 for the local mean sheet
    qa = q * P ** 2 * dy ** 2 * w / (np.maximum(rho, 1e-9) * s)  # radiated power per area (same units as the tensor q = 2 alpha W / rho)
    amp = np.sqrt(np.maximum(qa, 0.0))
    j = rel_sign * np.sign(dy) * amp * np.exp(-1j * phi_ref)
    return j * np.cos(phi), j * np.sin(phi), qa


def build(x, y, kind, grp, d, state_on, gb_mask=None, rel_sign=+1.0):
    """Hybrid aperture current: tensor cells everywhere except the GB sites."""
    gbm = (grp == "GB") if gb_mask is None else gb_mask
    jx, jy = d["jx"].copy(), d["jy"].copy()
    jx[gbm] = 0.0
    jy[gbm] = 0.0
    gx, gy, qa = gb_currents(d, gbm, state_on, rel_sign)
    jx[gbm], jy[gbm] = gx, gy
    # radiated power bookkeeping (per area units, total launched 2 pi): tensor q = amp^2 (removed at GB sites), GB q
    q_t = d["amp"] ** 2
    area = P * P
    eta_hyb_add = (np.sum(qa) * area - np.sum(q_t[gbm]) * area) / (2 * math.pi)
    return jx, jy, eta_hyb_add


def pattern_maps(x, y, jx, jy, green, theta, phi):
    dd = {"jx": jx, "jy": jy}
    return T.pattern(dd, x, y, green, theta, phi)


def analyse(et, ep, theta, phi, eta_conv, label, beams):
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    u_tot = np.abs(et) ** 2 + np.abs(ep) ** 2
    rh = np.abs(copol_component(et, ep, "rhcp")) ** 2
    lh = np.abs(copol_component(et, ep, "lhcp")) ** 2
    dth, dph = np.radians(theta[1] - theta[0]), np.radians(phi[1] - phi[0])
    prad = float(np.sum(u_tot * np.sin(np.radians(tg))) * dth * dph)
    g0 = 4 * math.pi / prad * T.LAUNCH_BARE * eta_conv * 10 ** (ANCHOR_DB / 10)  # gain per unit pattern intensity
    gt, gr, gl = g0 * u_tot, g0 * rh, g0 * lh
    out = {"label": label, "eta_conv": float(eta_conv)}
    i = np.unravel_index(np.argmax(np.where(np.hypot(tg - beams[0][0], ((pg - beams[0][1] + 180) % 360) - 180) <= 6.0, gr, -1)), gr.shape)
    ar = 20 * np.log10((np.sqrt(gr) + np.sqrt(gl)) / np.maximum(np.abs(np.sqrt(gr) - np.sqrt(gl)), 1e-12))
    out["beam1"] = {"peak_dir": [float(tg[i]), float(pg[i])], "rhcp_gain_dbi": float(10 * np.log10(gr[i])), "lhcp_at_peak_dbi": float(10 * np.log10(gl[i] + 1e-30)), "xpd_db": float(10 * np.log10(gr[i] / (gl[i] + 1e-30))), "ar_peak_db": float(ar[i])}
    cd = np.cos(np.radians(tg)) * np.cos(np.radians(tg[i])) + np.sin(np.radians(tg)) * np.sin(np.radians(tg[i])) * np.cos(np.radians(pg - pg[i]))
    dpk = np.degrees(np.arccos(np.clip(cd, -1, 1)))
    # 3 dB-AR beamwidth: extent of the region around the peak with AR <= 3 dB and RHCP within 3 dB of the peak, measured along theta and phi cuts through the peak
    ok = (ar <= 3.0) & (gr >= 0.5 * gr[i])
    def extent(line_ok, step):
        idx = np.nonzero(line_ok)[0]
        return None if idx.size == 0 else float((idx.max() - idx.min() + 1) * step)
    out["beam1"]["ar_le_3dB_width_theta_deg"] = extent(ok[:, i[1]], theta[1] - theta[0])
    out["beam1"]["ar_le_3dB_width_phi_deg"] = extent(ok[i[0], :], phi[1] - phi[0])
    out["beam1"]["ar_worst_in_3dB_beam_db"] = float(ar[(gr >= 0.5 * gr[i]) & (dpk <= 10.0)].max())
    # beam 2 (LP): total gain at the design direction, polarisation ellipse
    j2 = (np.argmin(abs(theta - beams[1][0])), np.argmin(abs(phi - beams[1][1])))
    out["beam2"] = {"total_gain_dbi_at_design_dir": float(10 * np.log10(gt[j2])), "rhcp_dbi": float(10 * np.log10(gr[j2] + 1e-30)), "lhcp_dbi": float(10 * np.log10(gl[j2] + 1e-30))}
    ex, ey = et[j2], ep[j2]
    # tilt of the polarisation ellipse in the (theta_hat, phi_hat) basis and the axial ratio / XPD relative to the phi-hat (LP-V) orientation
    tilt = 0.5 * math.degrees(math.atan2(2 * np.real(ex * np.conj(ey)), abs(ex) ** 2 - abs(ey) ** 2))
    out["beam2"]["pol_tilt_from_theta_hat_deg"] = float(tilt)
    out["beam2"]["xpd_phi_hat_vs_theta_hat_db"] = float(20 * math.log10(max(abs(ey), 1e-30) / max(abs(ex), 1e-30)))
    out["beam2"]["ar_db"] = float(ar[j2])
    # SLL of the total pattern outside the main lobes (first-null definition), per state
    mk = []
    u0 = 0.0
    ins = np.zeros(gt.shape, bool)
    for bm in beams:
        try:
            m_, uu = SLD.analyze(theta, phi, gt, peak=(float(bm[0]), float(bm[1])), return_mask=True)
            ins |= m_
            u0 = max(u0, uu)
        except Exception:
            pass
    TT = tg
    outs = np.where(ins | (TT > 90.0), -np.inf, gt)
    out["sll_total_db_model"] = float(10 * np.log10(outs.max() / max(gt[j2], gr[i], u0)))
    out["sll_total_db_with_offset"] = out["sll_total_db_model"] + SLL_OFFSET_DB
    # RHCP-only SLL of beam 1 (outside its main lobe)
    m1, u1 = SLD.analyze(theta, phi, gr, peak=(float(tg[i]), float(pg[i])), return_mask=True)
    out["sll_rhcp_beam1_db_model"] = float(10 * np.log10(np.where(m1 | (TT > 90.0), -np.inf, gr).max() / u1))
    out["sll_rhcp_beam1_db_with_offset"] = out["sll_rhcp_beam1_db_model"] + SLL_OFFSET_DB
    a1 = SLD.analyze(theta, phi, gr, peak=(float(tg[i]), float(pg[i])))
    out["beam1"]["hpbw_deg"] = float(a1["HPBW_mean_deg"])
    a2 = SLD.analyze(theta, phi, gt, peak=(float(beams[1][0]), float(beams[1][1])))
    out["beam2"]["hpbw_deg"] = float(a2["HPBW_mean_deg"])
    out["frac"] = {"within_30deg": float(np.sum(u_tot * np.sin(np.radians(tg)) * (tg <= 30)) / np.sum(u_tot * np.sin(np.radians(tg))))}
    out["_maps"] = (gt, gr, gl)
    return out


def run(rel_sign=+1.0, gb_mask_override=None, verbose=True, flip=None, cache=None):
    x, y, kind, grp = load_cells()
    table, m, d = cache if cache is not None else tensor_design(x, y)
    if flip is not None:
        d = dict(d)
        d['jx'] = d['jx'] * flip
        d['jy'] = d['jy'] * flip
    green = LayeredGreenFunction(T.SUB, F, backing="pec")
    theta = np.arange(0.0, 90.0 + 1e-9, 1.0)
    phi = np.arange(0.0, 360.0, 2.0)
    res = {"table_beta_over_k0_used": [float(np.nanmin(d["beta_over_k0"])), float(np.nanmean(d["beta_over_k0"])), float(np.nanmax(d["beta_over_k0"]))], "eta_conv_march": float(m["eta_conv"])}
    gbm = (grp == "GB") if gb_mask_override is None else gb_mask_override
    for name, on in (("GB_OFF", False), ("GB_ON", True)):
        jx, jy, add = build(x, y, kind, grp, d, on, gbm, rel_sign)
        et, ep = pattern_maps(x, y, jx, jy, green, theta, phi)
        eta = m["eta_conv"] + add
        res[name] = analyse(et, ep, theta, phi, eta, name, (BEAM1, BEAM2))
        res[name]["eta_hyb_add"] = float(add)
        if verbose:
            r = res[name]
            print(name, "eta_conv %.3f" % eta, json.dumps({k: v for k, v in r.items() if k not in ("_maps", "label")}, default=float), flush=True)
    return res, (x, y, kind, grp, d, m)


if __name__ == "__main__" and (len(sys.argv) < 2 or sys.argv[1] == "base"):
    for rs in (+1.0, -1.0):
        print("relative sign GB vs tensor:", rs)
        res, _ = run(rs)
        res = {k: ({a: b for a, b in v.items() if a != "_maps"} if isinstance(v, dict) else v) for k, v in res.items()}
        (RES / f"cp_hybrid_rel{int(rs):+d}.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")


# ------------------------------------------------------------------------------------------------ GB placement study (step 5)
def site_fields(x, y, d, green, beam, rel_sign, sites):
    """Far-field (E_theta, E_phi) at `beam` of the GB ON-OFF difference current of EACH candidate site alone (complex 2-vectors) and the copol (y-polarised) scalar."""
    from src.pattern.vector_aperture_pattern import _vector_far_field

    allm = np.zeros(len(x), bool)
    allm[:] = True
    gxo, gyo, _ = gb_currents(d, allm, False, rel_sign)
    gxn, gyn, _ = gb_currents(d, allm, True, rel_sign)
    out = np.zeros(len(sites), complex)
    th, ph = np.radians(beam[0]), np.radians(beam[1])
    u = np.array([np.cos(th) * np.sin(ph), np.cos(ph)])
    for n, i in enumerate(sites):
        et, ep = _vector_far_field(np.array([beam[0]]), np.array([beam[1]]), x[i:i + 1], y[i:i + 1], np.array([gxn[i] - gxo[i]]), np.array([gyn[i] - gyo[i]]), F, green, P, P)
        out[n] = et[0] * u[0] + ep[0] * u[1]
    return out


def placement(rel_sign=-1.0, counts=(180, 150, 120, 90)):
    x, y, kind, grp, d, m, = None, None, None, None, None, None
    x, y, kind, grp = load_cells()
    table, m, d = tensor_design(x, y)
    green = LayeredGreenFunction(T.SUB, F, backing="pec")
    n = len(x)
    allsites = np.arange(n)
    c = site_fields(x, y, d, green, BEAM2, rel_sign, allsites)
    cur = np.nonzero(grp == "GB")[0]
    S = c[cur].sum()
    res = {"rel_sign": rel_sign, "S_current_abs": float(abs(S)), "incoherent_sum_current": float(np.sqrt(np.sum(np.abs(c[cur]) ** 2))), "coherence_current": float(abs(S) / np.sqrt(np.sum(np.abs(c[cur]) ** 2)) / np.sqrt(len(cur)))}
    proj_cur = np.real(c * np.conj(S)) / abs(S)
    order_cur = cur[np.argsort(-proj_cur[cur])]
    # optimal sites among all: iterate reference phase
    ref = S
    for _ in range(6):
        pr = np.real(c * np.conj(ref))
        best = np.argsort(-pr)[:180]
        ref = c[best].sum()
    res["S_best180_abs"] = float(abs(ref))
    cases = {}
    for ncell in counts:
        for tag, sel in (("current_top", order_cur[:ncell]), ("optimal_sites", np.argsort(-np.real(c * np.conj(ref)))[:ncell])):
            mask = np.zeros(n, bool)
            mask[sel] = True
            rr, _ = run(rel_sign, gb_mask_override=mask, verbose=False)
            cases[f"{tag}_{ncell}cells_{2 * ncell}diodes"] = {k: rr[k] for k in ("GB_OFF", "GB_ON")}
            g = lambda s: (rr[s]["beam1"]["rhcp_gain_dbi"], rr[s]["beam1"]["ar_peak_db"], rr[s]["beam2"]["total_gain_dbi_at_design_dir"], rr[s]["sll_total_db_with_offset"])
            go, gn = g("GB_OFF"), g("GB_ON")
            print("%-34s rel %+d | OFF: b1 %.1f dBi AR %.2f b2 %.1f SLL %.1f | ON: b1 %.1f AR %.2f b2 %.1f SLL %.1f | contrast %.1f dB balance %.1f dB | |S| %.3g" % (
                f"{tag} {ncell} cells/{2 * ncell} diodes", rel_sign, go[0], go[1], go[2], go[3], gn[0], gn[1], gn[2], gn[3], gn[2] - go[2], gn[0] - gn[2], abs(c[sel].sum())), flush=True)
    res["cases"] = cases
    return res, (x, y, kind, grp, d, m, c, ref)


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "placement":
    for rs in (-1.0, +1.0):
        r, _ = placement(rs)
        for k, v in r["cases"].items():
            for s in v:
                v[s] = {a: b for a, b in v[s].items() if a != "_maps"}
        (RES / f"cp_hybrid_placement_rel{int(rs):+d}.json").write_text(json.dumps(r, indent=1, default=float), encoding="utf-8", newline="\n")


# ------------------------------------------------------------------------------------------------ joint optimisation (GB sites + tensor slot 90-degree rotations = sign flips)
def joint_opt(rel_sign=-1.0, ncell=180, mu=4.0, loss_budget_db=0.4, max_iter=400, log=print):
    from src.pattern.vector_aperture_pattern import _vector_far_field

    x, y, kind, grp = load_cells()
    cache = tensor_design(x, y)
    table, m, d = cache
    green = LayeredGreenFunction(T.SUB, F, backing="pec")
    n = len(x)
    th2, ph2 = np.radians(BEAM2[0]), np.radians(BEAM2[1])
    u2 = np.array([np.cos(th2) * np.sin(ph2), np.cos(ph2)])
    allm = np.ones(n, bool)
    gxo, gyo, _ = gb_currents(d, allm, False, rel_sign)
    gxn, gyn, _ = gb_currents(d, allm, True, rel_sign)
    jt = (d["jx"], d["jy"])

    def fld(jx, jy, beam):
        out = np.zeros((n, 2), complex)
        for i in range(n):
            et, ep = _vector_far_field(np.array([beam[0]]), np.array([beam[1]]), x[i:i + 1], y[i:i + 1], np.array([jx[i]]), np.array([jy[i]]), F, green, P, P)
            out[i] = (et[0], ep[0])
        return out

    v2 = {nm: fld(*j, BEAM2) for nm, j in (("t", jt), ("off", (gxo, gyo)), ("on", (gxn, gyn)))}
    v1 = {nm: fld(*j, BEAM1) for nm, j in (("t", jt), ("off", (gxo, gyo)), ("on", (gxn, gyn)))}
    c2 = {k: v[:, 0] * u2[0] + v[:, 1] * u2[1] for k, v in v2.items()}
    b1 = {k: np.array([circular_components(v[i, 0], v[i, 1])[0] for i in range(n)]) for k, v in v1.items()}  # RHCP complex
    gb = (grp == "GB").copy()
    if ncell != 180:
        order = np.argsort(-np.real(c2["on"] * np.conj(c2["on"][gb].sum())))
        gb = np.zeros(n, bool)
        gb[[i for i in order if grp[i] == "GB"][:ncell]] = True
    s = np.ones(n)
    B0 = abs(b1["t"][~(grp == "GB")].sum() + b1["off"][grp == "GB"].sum())

    def state(gb, s):
        FT, BT = np.sum((s * c2["t"])[~gb]), np.sum((s * b1["t"])[~gb])
        Fo, Fn = FT + c2["off"][gb].sum(), FT + c2["on"][gb].sum()
        Bo, Bn = BT + b1["off"][gb].sum(), BT + b1["on"][gb].sum()
        return Fo, Fn, Bo, Bn

    def cost(Fo, Fn, Bo, Bn):
        floor = (B0 ** 2) * 10 ** (-loss_budget_db / 10)
        pen = lambda B: np.maximum(floor - np.abs(B) ** 2, 0.0) * 50.0
        return -np.abs(Fn) ** 2 + mu * np.abs(Fo) ** 2 + (pen(Bo) + pen(Bn)) * 1.0

    scale = abs(c2["on"][gb].sum()) ** 2
    cur = cost(*state(gb, s)) / scale
    log("start: cost %.3f  F_OFF %.3g F_ON %.3g" % (cur, abs(state(gb, s)[0]), abs(state(gb, s)[1])))
    for it in range(max_iter):
        Fo, Fn, Bo, Bn = state(gb, s)
        T_idx = np.nonzero(~gb)[0]
        # all single flips of tensor cells
        dF = -2 * s[T_idx] * c2["t"][T_idx]
        dB = -2 * s[T_idx] * b1["t"][T_idx]
        cf = cost(Fo + dF, Fn + dF, Bo + dB, Bn + dB) / scale
        # swaps: GB j -> tensor (sign +1), tensor k -> GB
        G_idx = np.nonzero(gb)[0]
        # restrict swap candidates to keep the problem small: all pairs (j in G, k in T)
        dFo = (c2["t"][G_idx][:, None] - c2["off"][G_idx][:, None]) + (-s[T_idx] * c2["t"][T_idx] + c2["off"][T_idx])[None, :]
        dFn = (c2["t"][G_idx][:, None] - c2["on"][G_idx][:, None]) + (-s[T_idx] * c2["t"][T_idx] + c2["on"][T_idx])[None, :]
        dBo = (b1["t"][G_idx][:, None] - b1["off"][G_idx][:, None]) + (-s[T_idx] * b1["t"][T_idx] + b1["off"][T_idx])[None, :]
        dBn = (b1["t"][G_idx][:, None] - b1["on"][G_idx][:, None]) + (-s[T_idx] * b1["t"][T_idx] + b1["on"][T_idx])[None, :]
        cs = cost(Fo + dFo, Fn + dFn, Bo + dBo, Bn + dBn) / scale
        bf, bs = cf.min(), cs.min()
        if min(bf, bs) >= cur - 1e-9:
            break
        if bf <= bs:
            k = T_idx[int(np.argmin(cf))]
            s[k] *= -1
            cur = bf
        else:
            a, b = np.unravel_index(int(np.argmin(cs)), cs.shape)
            j, k = G_idx[a], T_idx[b]
            gb[j] = False
            s[j] = 1.0
            gb[k] = True
            s[k] = 1.0
            cur = bs
    Fo, Fn, Bo, Bn = state(gb, s)
    log("joint opt done after %d iterations: |F_OFF|^2 %.3g |F_ON|^2 %.3g (x%.2f of the start ON) | beam-1 amplitude^2 OFF/ON vs start: %.3f / %.3f | flips %d | GB sites changed %d" % (
        it, abs(Fo) ** 2, abs(Fn) ** 2, abs(Fn) ** 2 / scale, abs(Bo) ** 2 / B0 ** 2, abs(Bn) ** 2 / B0 ** 2, int(np.sum(s < 0)), int(np.sum(gb & (grp != "GB")))))
    return {"gb_mask": gb, "flip": s, "cache": cache, "x": x, "y": y, "grp": grp}


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "joint":
    out = {}
    for rs in (-1.0, +1.0):
        jo = joint_opt(rs)
        for eval_sign in (-1.0, +1.0):
            print("design for rel", rs, "evaluated under rel", eval_sign, flush=True)
            rr, _ = run(eval_sign, gb_mask_override=jo["gb_mask"], flip=jo["flip"], cache=jo["cache"])
            out[f"design_rel{int(rs):+d}_eval_rel{int(eval_sign):+d}"] = {k: {a: b for a, b in v.items() if a != "_maps"} for k, v in rr.items() if k in ("GB_OFF", "GB_ON")}
        np.savez(RES / f"cp_hybrid_joint_design_rel{int(rs):+d}.npz", gb=jo["gb_mask"], flip=jo["flip"], x=jo["x"], y=jo["y"])
    (RES / "cp_hybrid_joint.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")


def summarize_run(rr):
    o, n = rr["GB_OFF"], rr["GB_ON"]
    return {"b1_off_dbi": o["beam1"]["rhcp_gain_dbi"], "b1_on_dbi": n["beam1"]["rhcp_gain_dbi"], "ar_off": o["beam1"]["ar_peak_db"], "ar_on": n["beam1"]["ar_peak_db"], "ar_worst_3dB_beam_off": o["beam1"]["ar_worst_in_3dB_beam_db"],
            "b2_off_dbi": o["beam2"]["total_gain_dbi_at_design_dir"], "b2_on_dbi": n["beam2"]["total_gain_dbi_at_design_dir"], "contrast_db": n["beam2"]["total_gain_dbi_at_design_dir"] - o["beam2"]["total_gain_dbi_at_design_dir"],
            "balance_b1_minus_b2_on_db": n["beam1"]["rhcp_gain_dbi"] - n["beam2"]["total_gain_dbi_at_design_dir"], "sll_off_with_offset": o["sll_total_db_with_offset"], "sll_on_with_offset": n["sll_total_db_with_offset"],
            "hpbw_b1": o["beam1"]["hpbw_deg"], "hpbw_b2": n["beam2"]["hpbw_deg"], "b2_xpd_phi_vs_theta_db": n["beam2"]["xpd_phi_hat_vs_theta_hat_db"], "b2_tilt_deg": n["beam2"]["pol_tilt_from_theta_hat_deg"]}


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "sweep":
    out = {}
    for ncell, lb, mu in ((180, 0.4, 4.0), (150, 0.4, 4.0), (120, 0.4, 4.0), (90, 0.4, 4.0), (180, 1.0, 4.0), (180, 2.0, 4.0)):
        jo = joint_opt(-1.0, ncell=ncell, mu=mu, loss_budget_db=lb, log=lambda *a: None)
        rr, _ = run(-1.0, gb_mask_override=jo["gb_mask"], flip=jo["flip"], cache=jo["cache"], verbose=False)
        sm = summarize_run(rr)
        sm.update({"ncell": int(jo["gb_mask"].sum()), "diodes": int(2 * jo["gb_mask"].sum()), "loss_budget_db": lb, "mu": mu, "flips": int(np.sum(jo["flip"] < 0)), "gb_sites_changed": int(np.sum(jo["gb_mask"] & (jo["grp"] != "GB")))})
        out[f"n{ncell}_lb{lb}_mu{mu}"] = sm
        print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in sm.items()}, flush=True)
        np.savez(RES / f"cp_hybrid_design_n{ncell}_lb{lb}.npz", gb=jo["gb_mask"], flip=jo["flip"], x=jo["x"], y=jo["y"])
    (RES / "cp_hybrid_sweep.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
