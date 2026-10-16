"""Phase 27 item 5: HFSS-export comparison for the tensor-CP boards (single hexagon, four hexagons) -- moved from
`scripts/research/phase26_compare_single_hex.py` / `phase26_compare_fourhex.py` (now thin wrappers). No modelling: every
number is read from the user's exports or from `predictions.json`; named states, never a bare pass/fail.

    python scripts/research/phase26_compare_single_hex.py --gain gain.csv [--s11 file.s1p|.s4p [--port N]]
        [--params antenna_params.txt] [--ar ar_cut.csv] [--ez ez_line.csv] [--pin scripted|connector|<label>]
        [--pred predictions.json] [--out report.json]

No modelling: every number is read from the user's exports or from `predictions.json`; a missing export
gives `NOT_PROVIDED`, a missing prediction `NO_PREDICTION`. The feed may be the scripted pin of the AEDT
script or the user's own connector model from `hex150_gap0`: the comparison is the same; the label is only
recorded (the loaded-slab launch fraction, hence the predicted gain bracket, depends on the feed).

Inputs and tolerated formats
* gain CSV: HFSS "Realized Gain" export with columns containing Freq, Phi, Theta, RealizedGainLHCP, RealizedGainRHCP
  (dB), any theta/phi step (1 deg expected), phi in -180..180 or 0..360, optional extra frequencies (the one
  nearest 10.54 GHz is used).
* antenna parameters: the HFSS "Antenna Parameters" text table (rows `Peak Realized Gain  f1 v1 f2 v2`); the last
  value column is used.
* AR cut: CSV with a theta (or angle) column and a column containing "Axial" or "AR" (dB).
* Ez line: CSV with a distance column (mm or m) and `Mag_E` or complex Ez (re, im) columns -> alpha fit.
"""

import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]  # repo root (src/hologram/ -> parents[2])
PRED = ROOT / "hfss_exports/hex150_tensor_cp_single/predictions.json"
BEAM = (30.0, 225.0)
F_DESIGN = 10.54
GAIN_FALSIFIER_DBI = 19.0
GAIN_TOL_DB = 1.5
POINTING_TOL_DEG = 2.0


# ------------------------------------------------------------------------------------------ parsers
def read_gain_csv(path, f_ghz=F_DESIGN):
    rows = list(csv.reader(Path(path).read_text().splitlines()))
    head = [h.strip().lower() for h in rows[0]]

    def col(*keys):
        for i, h in enumerate(head):
            if all(k in h for k in keys):
                return i
        raise ValueError(f"{path}: no column containing {keys} in {rows[0]}")

    i_f = col("freq") if any("freq" in h for h in head) else None
    i_p, i_t = col("phi"), col("theta")
    i_l, i_r = col("lhcp"), col("rhcp")
    a = np.array(
        [[float(r[i]) for i in ([i_f] if i_f is not None else []) + [i_p, i_t, i_l, i_r]] for r in rows[1:] if r]
    )
    if i_f is not None:
        fs = np.unique(a[:, 0])
        f_use = float(fs[np.argmin(np.abs(fs - f_ghz))])
        a = a[np.abs(a[:, 0] - f_use) < 1e-9][:, 1:]
    else:
        f_use = None
    phi = np.mod(a[:, 0], 360.0)
    th, gl, gr = a[:, 1], a[:, 2], a[:, 3]
    u_th, u_ph = np.unique(np.round(th, 6)), np.unique(np.round(phi, 6))
    pl = np.full((u_ph.size, u_th.size), np.nan)
    pr = np.full_like(pl, np.nan)
    ip, it = np.searchsorted(u_ph, np.round(phi, 6)), np.searchsorted(u_th, np.round(th, 6))
    pl[ip, it] = 10 ** (gl / 10.0)
    pr[ip, it] = 10 ** (gr / 10.0)
    return {"theta": u_th, "phi": u_ph, "pl": pl, "pr": pr, "freq_ghz": f_use}


def read_params(path):
    """`{label: last numeric value}` of an HFSS antenna-parameter table (`Peak Realized Gain   4.09   97.16`)."""
    out = {}
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        m = re.match(r"^\s*([A-Za-z][A-Za-z /()0-9_-]*?)\s{2,}(.*)$", line)
        if not m:
            continue
        nums = re.findall(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?", m.group(2))
        if nums:
            out[m.group(1).strip()] = float(nums[-1])
    return out


def read_s_parameter(path, port=1, f_ghz=F_DESIGN):
    """Touchstone `.sNp` (S[port, port]) or an HFSS sweep CSV (`Freq [GHz]`, `dB(S(1,1))`)."""
    if str(path).lower().endswith(".csv"):
        rows = list(csv.reader(Path(path).read_text().splitlines()))[1:]
        a = np.array([[float(v) for v in r[:2]] for r in rows if r])
        f, db = a[:, 0], a[:, 1]
    else:
        from ..phase9.hfss_reference import read_touchstone

        net = read_touchstone(path)
        s = net.component(port, port)
        f = net.frequency_hz / 1e9
        db = 20 * np.log10(np.maximum(np.abs(s), 1e-12))
    i = int(np.argmin(np.abs(f - f_ghz)))
    below = f[db <= -10.0]
    out = {
        "s_db_at_design": float(db[i]),
        "f_nearest_ghz": float(f[i]),
        "min_s_db": float(db.min()),
        "f_at_min_ghz": float(f[np.argmin(db)]),
    }
    if below.size:
        out["band_below_minus10_db_ghz"] = [float(below.min()), float(below.max())]
        out["design_inside_minus10_db"] = bool(db[i] <= -10.0)
    else:
        out["band_below_minus10_db_ghz"] = None
        out["design_inside_minus10_db"] = False
    return out


def read_ar_sweep(path, theta0=BEAM[0], f_ghz=F_DESIGN):
    """HFSS AR export with `Freq`, `Phi`, `Theta`, `dB(AxialRatioValue)` columns (a frequency sweep at one phi cut):
    AR vs frequency at the beam theta, the 3 dB AR band there, and the cut at the frequency nearest `f_ghz`.
    Returns None if the file has no frequency column."""
    rows = list(csv.reader(Path(path).read_text().splitlines()))
    head = [h.strip().lower() for h in rows[0]]
    if not any("freq" in h for h in head):
        return None
    ifq = next(i for i, h in enumerate(head) if "freq" in h)
    ipp = next(i for i, h in enumerate(head) if "phi" in h)
    ith = next(i for i, h in enumerate(head) if "theta" in h)
    iar = next(i for i, h in enumerate(head) if "axial" in h)
    a = np.array([[float(r[i]) for i in (ifq, ipp, ith, iar)] for r in rows[1:] if r])
    fs = np.unique(a[:, 0])
    ar_at_theta = []
    for f in fs:
        m = a[a[:, 0] == f]
        ar_at_theta.append(float(np.interp(theta0, m[:, 2], m[:, 3])))
    ar_at_theta = np.array(ar_at_theta)
    f_use = float(fs[np.argmin(np.abs(fs - f_ghz))])
    cut = a[a[:, 0] == f_use]
    ok = fs[ar_at_theta <= 3.0]
    return {
        "phi_cut_deg": float(a[0, 1]),
        "frequencies_ghz": fs.tolist(),
        "ar_at_beam_theta_db": ar_at_theta.tolist(),
        "frequency_used_ghz": f_use,
        "ar_at_beam_theta_at_used_freq_db": float(np.interp(theta0, cut[:, 2], cut[:, 3])),
        "min_ar_in_cut_db": float(cut[:, 3].min()),
        "min_ar_theta_deg": float(cut[np.argmin(cut[:, 3]), 2]),
        "freqs_with_ar_le_3db_at_beam_theta_ghz": ok.tolist(),
        "contiguous_3db_band_around_design_ghz": _contiguous_band(fs, ar_at_theta <= 3.0, f_use),
    }


def _contiguous_band(fs, ok, f0):
    i0 = int(np.argmin(np.abs(fs - f0)))
    if not ok[i0]:
        return None
    lo = hi = i0
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    while hi < fs.size - 1 and ok[hi + 1]:
        hi += 1
    return [float(fs[lo]), float(fs[hi])]


def read_re_field(path):
    """`{label: (V, theta_deg, phi_deg)}` from the `rE Field` block of an HFSS antenna-parameter table (last column)."""
    out = {}
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        m = re.match(
            r"^\s*([A-Za-z0-9/ ]+?)\s+(?:[\d.]+GHz\s+)?([-+]?\d+\.?\d*)\s*V\s+\((\d+)deg,(\d+)deg\)\s*$", line.rstrip()
        )
        if m:
            out[m.group(1).strip()] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
    return out


def read_ar_cut(path):
    rows = list(csv.reader(Path(path).read_text().splitlines()))
    head = [h.strip().lower() for h in rows[0]]
    it = next(i for i, h in enumerate(head) if "theta" in h or "angle" in h or "deg" in h)
    ia = next(i for i, h in enumerate(head) if "axial" in h or re.search(r"\bar\b", h) or h.startswith("db(ar"))
    a = np.array([[float(r[it]), float(r[ia])] for r in rows[1:] if r])
    return a[:, 0], a[:, 1]


def fit_decay(path):
    rows = list(csv.reader(Path(path).read_text().splitlines()))
    head = [h.strip().lower() for h in rows[0]]
    data = np.array([[float(v) for v in r] for r in rows[1:] if r])
    d = data[:, 0]
    rho = d * (1e-3 if d.max() > 5.0 else 1.0)  # mm or m
    im = next((i for i, h in enumerate(head) if "mag" in h), None)
    if im is not None:
        e = np.abs(data[:, im])
    else:
        ir = next(i for i, h in enumerate(head) if "re" in h and "ez" in h)
        ii = next(i for i, h in enumerate(head) if "im" in h and "ez" in h)
        e = np.hypot(data[:, ir], data[:, ii])
    m = (rho > 0.02) & (rho < 0.13) & (e > 0)
    # |E| ~ rho^-1/2 exp(-alpha rho): ln(|E| sqrt(rho)) linear in rho
    y = np.log(e[m] * np.sqrt(rho[m]))
    p = np.polyfit(rho[m], y, 1)
    return {
        "alpha_np_per_m": float(-p[0]),
        "points": int(m.sum()),
        "rho_range_mm": [float(rho[m].min() * 1e3), float(rho[m].max() * 1e3)],
    }


# ------------------------------------------------------------------------------------------ analysis
def _angle_to(theta_deg, phi_deg, t0, p0):
    t, p = np.radians(theta_deg), np.radians(phi_deg)
    t0, p0 = np.radians(t0), np.radians(p0)
    c = np.cos(t) * np.cos(t0) + np.sin(t) * np.sin(t0) * np.cos(p - p0)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


def analyse(g, beam=None):
    beam = BEAM if beam is None else beam
    th, ph, pl, pr = g["theta"], g["phi"], g["pl"], g["pr"]
    tg, pg = np.meshgrid(th, ph)
    valid = np.isfinite(pl) & np.isfinite(pr)
    dth = np.gradient(th) if th.size > 1 else np.array([1.0])
    dph = float(np.median(np.diff(ph))) if ph.size > 1 else 360.0
    w = np.sin(np.radians(tg)) * np.radians(dth)[None, :] * np.radians(dph)
    w = np.where(valid, w, 0.0)
    pl0, pr0 = np.where(valid, pl, 0.0), np.where(valid, pr, 0.0)
    tot = pl0 + pr0
    p_int = float(np.sum(tot * w) / (4 * np.pi))  # = total efficiency (realized gain integrates to e_total)
    i = np.unravel_index(np.nanargmax(np.where(valid, pr, -1)), pr.shape)
    peak_r = float(pr[i])
    peak_t, peak_p = float(tg[i]), float(pg[i])
    d_beam = _angle_to(tg, pg, *beam)
    d_peak = _angle_to(tg, pg, peak_t, peak_p)
    out = {
        "peak_rhcp_gain_dbi": 10 * np.log10(peak_r),
        "peak_theta_phi_deg": [peak_t, peak_p],
        "pointing_error_deg": float(_angle_to(peak_t, peak_p, *beam)),
        "xpd_at_peak_db": float(10 * np.log10(peak_r / pl[i])),
        "integrated_realized_gain_over_4pi_total_efficiency": p_int,
        "grid_complete": bool(valid.all()),
    }
    er, el = np.sqrt(pr0), np.sqrt(pl0)
    ar = 20 * np.log10((er + el) / np.maximum(np.abs(er - el), 1e-12))
    beam3 = valid & (pr >= 0.5 * peak_r) & (d_peak <= 10.0)
    out["ar_at_peak_db"] = float(ar[i])
    out["ar_max_over_minus3db_beam_db"] = float(ar[beam3].max())
    out["minus3db_beam_points"] = int(beam3.sum())
    side = valid & (d_peak > 15.0)
    out["rhcp_sll_db"] = float(10 * np.log10(pr[side].max() / peak_r))

    # power budget: fractions of the integrated (radiated) power, regions disjoint
    def frac(mask, arr=tot):
        return float(np.sum(arr * w * mask) / (4 * np.pi) / p_int)

    upper = tg <= 90.0
    main = d_beam <= 8.0
    out["budget_of_radiated"] = {
        "main_beam_le_8deg_total": frac(main),
        "main_beam_le_8deg_rhcp": frac(main, pr0),
        "near_sidelobes_8_to_20deg": frac((d_beam > 8.0) & (d_beam <= 20.0)),
        "far_sidelobes_theta_le_70": frac((d_beam > 20.0) & (tg <= 70.0)),
        "grazing_70_to_90": frac((tg > 70.0) & upper),
        "back_hemisphere_gt_90": frac(tg > 90.0),
        "lhcp_total_all_angles": frac(np.ones_like(tg, dtype=bool), pl0),
        "rhcp_total_all_angles": frac(np.ones_like(tg, dtype=bool), pr0),
        "rhcp_within_beam_over_total_within_beam": float(np.sum(pr0 * w * main) / max(np.sum(tot * w * main), 1e-30)),
    }
    return out, (tg, pg, ar)


def compare_params_only(args, pred):
    """Checks that need only the antenna-parameter table (no gain CSV): peak realized gain, directivity, pointing,
    a LOWER BOUND of the XPD (RHCP peak vs the largest LHCP anywhere) and the efficiencies."""
    p = read_params(args.params)
    re_f = read_re_field(args.params)
    ap = pred["aperture_pattern"]
    gp = pred["realized_gain_rhcp_dbi"]
    lo, hi = gp["launch_bare_0.592"], gp["launch_loaded_0.66_to_0.88"][1]
    rep = {
        "labels": ["HFSS_USER_EXPORT", "COMPARISON_ONLY_NO_MODELLING"],
        "pin_model": args.pin,
        "gain_csv": "NOT_PROVIDED",
        "note": "antenna-parameter table only: no gain CSV, so no pattern budget, SLL or beam-purity integrals",
        "measured": {},
        "checks": {},
    }
    pk = 10 * np.log10(p["Peak Realized Gain"])
    state = (
        "MODEL_FALSIFIED_PEAK_GAIN_BELOW_1BIT_BOARD"
        if pk < GAIN_FALSIFIER_DBI
        else "GAIN_BELOW_PREDICTION"
        if pk < lo - GAIN_TOL_DB
        else "GAIN_ABOVE_PREDICTION"
        if pk > hi + GAIN_TOL_DB
        else "GAIN_WITHIN_PREDICTED_RANGE"
    )
    rep["checks"]["peak_gain"] = {
        "hfss_dbi": pk,
        "predicted_dbi_bare_launch": lo,
        "predicted_dbi_loaded_launch_max": hi,
        "tolerance_db": GAIN_TOL_DB,
        "falsifier_dbi": GAIN_FALSIFIER_DBI,
        "state": state,
        "delta_vs_bare_launch_db": pk - lo,
    }
    d_dbi = 10 * np.log10(p["Peak Directivity"])
    rep["checks"]["directivity_dbi"] = {
        "hfss": d_dbi,
        "model_aperture_only": ap["directivity_ap_rhcp_dbi"],
        "ratio_to_aperture_directivity": 10 ** ((d_dbi - ap["directivity_ap_rhcp_dbi"]) / 10),
        "model_main_beam_share_of_accepted": pred["accepted_power_budget"]["launch_bare_0.592"]["main_beam_8deg"],
    }
    if "RHCP" in re_f:
        _, t, ph = re_f["RHCP"]
        err = float(_angle_to(t, ph, *BEAM))
        rep["checks"]["pointing"] = {
            "hfss_peak_rhcp_theta_phi_deg": [t, ph],
            "hfss_deg": err,
            "tolerance_deg": POINTING_TOL_DEG,
            "state": "OK" if err <= POINTING_TOL_DEG else "POINTING_FALSIFIED",
        }
    if "RHCP" in re_f and "LHCP" in re_f:
        rep["checks"]["xpd_lower_bound_db"] = {
            "hfss_rhcp_peak_over_largest_lhcp_anywhere_db": 20 * np.log10(re_f["RHCP"][0] / re_f["LHCP"][0]),
            "lhcp_max_at_theta_phi_deg": [re_f["LHCP"][1], re_f["LHCP"][2]],
            "predicted_xpd_at_peak_db": ap["xpd_at_peak_db"],
            "brief_target_db": 25.0,
            "note": "XPD at the beam peak is >= this bound (the largest LHCP lies elsewhere)",
        }
    for k in ("Total Efficiency", "Radiation Efficiency"):
        if k in p:
            rep["checks"][k.lower().replace(" ", "_")] = p[k]
    rep["antenna_params_frequency_note"] = "table computed at 10.5 GHz (design 10.54)"
    rep["antenna_params"] = p
    return rep


def compare(args):
    pred = json.loads(Path(args.pred).read_text())
    if not args.gain:
        rep = compare_params_only(args, pred)
        _attach_sweeps(rep, args)
        return rep
    g = read_gain_csv(args.gain)
    a, (tg, pg, ar) = analyse(g)
    rep = {
        "labels": ["HFSS_USER_EXPORT", "COMPARISON_ONLY_NO_MODELLING"],
        "pin_model": args.pin,
        "gain_csv": str(args.gain),
        "frequency_used_ghz": g["freq_ghz"],
        "measured": a,
        "checks": {},
    }
    ap = pred["aperture_pattern"]
    gp = pred["realized_gain_rhcp_dbi"]
    lo, hi = gp["launch_bare_0.592"], gp["launch_loaded_0.66_to_0.88"][1]
    chk = rep["checks"]
    pk = a["peak_rhcp_gain_dbi"]
    if pk < GAIN_FALSIFIER_DBI:
        state = "MODEL_FALSIFIED_PEAK_GAIN_BELOW_1BIT_BOARD"
    elif pk < lo - GAIN_TOL_DB:
        state = "GAIN_BELOW_PREDICTION"
    elif pk > hi + GAIN_TOL_DB:
        state = "GAIN_ABOVE_PREDICTION"
    else:
        state = "GAIN_WITHIN_PREDICTED_RANGE"
    chk["peak_gain"] = {
        "hfss_dbi": pk,
        "predicted_dbi_bare_launch": lo,
        "predicted_dbi_loaded_launch_max": hi,
        "tolerance_db": GAIN_TOL_DB,
        "falsifier_dbi": GAIN_FALSIFIER_DBI,
        "state": state,
    }
    chk["pointing"] = {
        "hfss_deg": a["pointing_error_deg"],
        "tolerance_deg": POINTING_TOL_DEG,
        "state": "OK" if a["pointing_error_deg"] <= POINTING_TOL_DEG else "POINTING_FALSIFIED",
    }
    chk["xpd_at_peak_db"] = {
        "hfss": a["xpd_at_peak_db"],
        "predicted": ap["xpd_at_peak_db"],
        "brief_target": 25.0,
        "delta_db": a["xpd_at_peak_db"] - ap["xpd_at_peak_db"],
    }
    chk["ar_at_peak_db"] = {
        "hfss": a["ar_at_peak_db"],
        "predicted": ap["ar_at_peak_db"],
        "brief_target_le": 3.0,
        "worst_in_minus3db_beam_hfss": a["ar_max_over_minus3db_beam_db"],
        "worst_predicted": ap["ar_max_over_minus3db_beam_db"],
    }
    bb = a["budget_of_radiated"]
    chk["rhcp_fraction_in_beam"] = {
        "hfss": bb["rhcp_within_beam_over_total_within_beam"],
        "predicted": ap["within_8deg_rhcp_fraction"] / ap["within_8deg_fraction"],
        "falsifier_lt": 0.8,
        "state": "OK" if bb["rhcp_within_beam_over_total_within_beam"] >= 0.8 else "BEAM_PURITY_FALSIFIED",
    }
    chk["rhcp_sll_db"] = {"hfss": a["rhcp_sll_db"], "predicted": ap["rhcp_sll_db"], "brief_target_le": -15.0}
    # power budget in the model's normalisation: share of the ACCEPTED power. radiation efficiency ~ 1: of-radiated = of-accepted
    pb = pred["accepted_power_budget"]["launch_bare_0.592"]
    chk["power_budget_fractions_of_radiated"] = {
        "hfss": bb,
        "model_main_beam_8deg": pb["main_beam_8deg"],
        "model_lhcp_aperture_only": pb["lhcp_total_in_aperture_pattern"],
        "model_grazing_plus_back_upper_bound": pb["grazing_plus_back_upper_bound_unmodelled"],
        "hfss_grazing_plus_back": bb["grazing_70_to_90"] + bb["back_hemisphere_gt_90"],
        "note": "The model gives no pattern to the pin's direct space wave (0.41) or the unconverted remainder (0.27): they are lumped "
        "into the grazing+back upper bound. HFSS LHCP includes the pin's direct radiation, the model's LHCP does not.",
    }
    chk["integral_consistency"] = {
        "integrated_gain_over_4pi": a["integrated_realized_gain_over_4pi_total_efficiency"],
        "grid_complete": a["grid_complete"],
    }
    if args.params:
        p = read_params(args.params)
        rep["antenna_params"] = p
        if "Peak Realized Gain" in p:
            chk["peak_realized_gain_table_dbi"] = 10 * np.log10(
                p["Peak Realized Gain"]
            )  # HFSS tables list linear values
        for k in ("Total Efficiency", "Radiation Efficiency", "Accepted Power", "Radiated Power", "Incident Power"):
            if k in p:
                rep.setdefault("antenna_params_selected", {})[k] = p[k]
    else:
        rep["antenna_params"] = "NOT_PROVIDED"
    _attach_sweeps(rep, args)
    return rep


def _attach_sweeps(rep, args):
    if args.s11:
        rep["s_parameter"] = read_s_parameter(args.s11, args.port)
        rep["s_parameter"]["prediction"] = "NO_PREDICTION (match not modelled; the pin is the scripted/connector feed)"
        rep["s_parameter"]["brief_target"] = "S11 <= -10 dB at 10.54 GHz (launcher redesign constraint)"
    else:
        rep["s_parameter"] = "NOT_PROVIDED"
    if args.ar:
        sw = read_ar_sweep(args.ar)
        if sw is not None:
            rep["ar_sweep"] = sw
            rep["ar_cut"] = "SEE ar_sweep"
        else:
            t, arv = read_ar_cut(args.ar)
            rep["ar_cut"] = {"min_db": float(arv.min()), "at_theta_deg": float(t[np.argmin(arv)])}
    else:
        rep["ar_cut"] = "NOT_PROVIDED"
    _ez(rep, args)


def _ez(rep, args):
    if args.ez:
        d = fit_decay(args.ez)
        d["predicted_alpha_range_np_per_m"] = [2.0, 4.0]
        d["state"] = "WITHIN_PREDICTED_RANGE" if 2.0 <= d["alpha_np_per_m"] <= 4.0 else "OUTSIDE_PREDICTED_RANGE"
        rep["surface_wave_decay"] = d
    else:
        rep["surface_wave_decay"] = "NOT_PROVIDED"


def markdown_params_only(rep):
    c = rep["checks"]
    g, d = c["peak_gain"], c["directivity_dbi"]
    L = [
        f"# hex150_tensor_cp_single: HFSS vs prediction (pin model: {rep['pin_model']}) -- antenna-parameter table only",
        "",
        "| Quantity | HFSS | Predicted | Target / falsifier | State |",
        "|---|---|---|---|---|",
        f"| Peak RHCP realized gain [dBi] | {g['hfss_dbi']:.2f} | {g['predicted_dbi_bare_launch']:.1f} (bare launch) .. {g['predicted_dbi_loaded_launch_max']:.1f} (loaded) | >= 24; falsified < {g['falsifier_dbi']:g} | {g['state']} ({g['delta_vs_bare_launch_db']:+.2f} dB vs bare) |",
        f"| Peak directivity [dBi] | {d['hfss']:.2f} | aperture-only {d['model_aperture_only']:.1f} | | HFSS/aperture = {d['ratio_to_aperture_directivity']:.2f} (model main-beam share {d['model_main_beam_share_of_accepted']:.2f}) |",
    ]
    if "pointing" in c:
        L.append(
            f"| Peak RHCP direction [deg] | {c['pointing']['hfss_peak_rhcp_theta_phi_deg']} | 30, 224 | within {POINTING_TOL_DEG:g} of 30, 225 | {c['pointing']['state']} ({c['pointing']['hfss_deg']:.1f}) |"
        )
    if "xpd_lower_bound_db" in c:
        x = c["xpd_lower_bound_db"]
        L.append(
            f"| XPD at peak [dB] | >= {x['hfss_rhcp_peak_over_largest_lhcp_anywhere_db']:.1f} (bound) | {x['predicted_xpd_at_peak_db']:.1f} | >= 25 | |"
        )
    if "total_efficiency" in c:
        L.append(
            f"| Total / radiation efficiency | {c['total_efficiency']:.3f} / {c.get('radiation_efficiency', float('nan')):.3f} | | | |"
        )
    L.append("")
    for k in ("s_parameter", "ar_sweep", "ar_cut", "surface_wave_decay"):
        if k in rep:
            L.append(f"- {k}: `{json.dumps(rep[k], default=float)}`")
    L.append("")
    L.append(
        "Power budget (beam / sidelobes / grazing / back / LHCP): NOT_AVAILABLE_FROM_CURRENT_DATA -- needs the realized-gain RHCP/LHCP CSV."
    )
    return "\n".join(L) + "\n"


def markdown(rep):
    if rep["gain_csv"] == "NOT_PROVIDED":
        return markdown_params_only(rep)
    c, m = rep["checks"], rep["measured"]
    b = m["budget_of_radiated"]
    L = [f"# hex150_tensor_cp_single: HFSS vs prediction (pin model: {rep['pin_model']})", ""]
    L.append("| Quantity | HFSS | Predicted | Target / falsifier | State |")
    L.append("|---|---|---|---|---|")
    g = c["peak_gain"]
    L.append(
        f"| Peak RHCP realized gain [dBi] | {g['hfss_dbi']:.2f} | {g['predicted_dbi_bare_launch']:.1f} (bare) .. {g['predicted_dbi_loaded_launch_max']:.1f} (loaded) | >= 24; falsified < {g['falsifier_dbi']:g} | {g['state']} |"
    )
    L.append(
        f"| Peak direction (theta, phi) [deg] | {m['peak_theta_phi_deg'][0]:.0f}, {m['peak_theta_phi_deg'][1]:.0f} | 30, 224 | within {POINTING_TOL_DEG:g} deg of 30, 225 | {c['pointing']['state']} ({m['pointing_error_deg']:.1f} deg) |"
    )
    L.append(
        f"| XPD at peak [dB] | {m['xpd_at_peak_db']:.1f} | {c['xpd_at_peak_db']['predicted']:.1f} | >= 25 | delta {c['xpd_at_peak_db']['delta_db']:+.1f} dB |"
    )
    L.append(
        f"| AR at peak / worst in -3 dB beam [dB] | {m['ar_at_peak_db']:.2f} / {m['ar_max_over_minus3db_beam_db']:.2f} | {c['ar_at_peak_db']['predicted']:.2f} / {c['ar_at_peak_db']['worst_predicted']:.2f} | <= 3 | |"
    )
    L.append(
        f"| RHCP / total power within 8 deg | {c['rhcp_fraction_in_beam']['hfss']:.3f} | {c['rhcp_fraction_in_beam']['predicted']:.3f} | falsified < 0.8 | {c['rhcp_fraction_in_beam']['state']} |"
    )
    L.append(f"| RHCP SLL [dB] | {m['rhcp_sll_db']:.1f} | {c['rhcp_sll_db']['predicted']:.1f} | <= -15 | |")
    L.append("")
    L.append("## Power budget (fractions of the radiated power; HFSS radiation efficiency ~ 1 so ~ share of accepted)")
    L.append("| Region | HFSS | Model |")
    L.append("|---|---|---|")
    pb = c["power_budget_fractions_of_radiated"]
    L.append(
        f"| Main beam <= 8 deg (RHCP) | {b['main_beam_le_8deg_rhcp']:.3f} ({b['main_beam_le_8deg_total']:.3f} total) | {pb['model_main_beam_8deg']:.3f} |"
    )
    L.append(f"| Near sidelobes 8-20 deg | {b['near_sidelobes_8_to_20deg']:.3f} | in aperture pattern only |")
    L.append(f"| Far sidelobes (theta <= 70) | {b['far_sidelobes_theta_le_70']:.3f} | in aperture pattern only |")
    L.append(
        f"| Grazing 70-90 deg | {b['grazing_70_to_90']:.3f} | {pb['model_grazing_plus_back_upper_bound']:.3f} (upper bound incl. back, not modelled) |"
    )
    L.append(f"| Back hemisphere | {b['back_hemisphere_gt_90']:.3f} | (in the above) |")
    L.append(
        f"| LHCP, all angles | {b['lhcp_total_all_angles']:.3f} | {pb['model_lhcp_aperture_only']:.3f} (aperture only) |"
    )
    L.append("")
    for k in ("s_parameter", "ar_sweep", "ar_cut", "surface_wave_decay", "antenna_params_selected"):
        if k in rep:
            L.append(f"- {k}: `{json.dumps(rep[k], default=float)}`")
    L.append("")
    L.append(pb["note"])
    return "\n".join(L) + "\n"


def main_single(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gain")
    ap.add_argument("--s11")
    ap.add_argument("--port", type=int, default=1)
    ap.add_argument("--params")
    ap.add_argument("--ar")
    ap.add_argument("--ez")
    ap.add_argument("--pin", default="scripted")
    ap.add_argument("--pred", default=str(PRED))
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    rep = compare(args)
    txt = markdown(rep)
    print(txt)
    if args.out:
        Path(args.out).write_text(json.dumps(rep, indent=1, default=float))
        Path(args.out).with_suffix(".md").write_text(txt)
    return rep


# ------------------------------------------------------------------------------------ four hexagons (board A / B)
PRED4 = ROOT / "hfss_exports/hex4_tensor_cp_static/predictions.json"
PHI4 = {1: 45.0, 2: 135.0, 3: 225.0, 4: 315.0}
GAIN4_TOL_DB = 1.0
PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (3, 4)]


def _kv(items):
    out = {}
    for it in items or []:
        k, v = it.split("=", 1)
        out[int(k)] = v
    return out


def port_report(port, gain_path, params_path, pred):
    beam = (30.0, PHI4[port])
    g = read_gain_csv(gain_path, F_DESIGN)
    a, (tg, pg, _ar) = analyse(g, beam=beam)
    th, ph, pl, pr = g["theta"], g["phi"], g["pl"], g["pr"]
    tgrid, pgrid = np.meshgrid(th, ph)
    peak = float(np.nanmax(pr))
    d_beam = _angle_to(tgrid, pgrid, *beam)
    ok = np.isfinite(pl) & (d_beam > 9.0) & (tgrid <= 90.0)
    k = np.unravel_index(np.argmax(np.where(ok, pl, -1.0)), pl.shape)
    lobe = {
        "level_rel_rhcp_peak_db": float(10 * np.log10(pl[k] / peak)),
        "theta_phi_deg": [float(tgrid[k]), float(pgrid[k])],
        "deg_from_beam": float(d_beam[k]),
    }
    rep = {"measured": a, "strongest_lhcp_lobe_outside_9deg": lobe, "freq_ghz": g["freq_ghz"]}
    if params_path:
        rep["antenna_params"] = read_params(params_path)
    p = pred.get("per_port_alone_by_edge_power_scale")
    if not p:
        rep["checks"] = {"state": "NO_PREDICTION"}
        return rep
    p1, p17 = p["s=1.0"][str(port)], p["s=1.7"][str(port)]
    lo, hi = (
        min(p1["realized_gain_rhcp_dbi_anchored"], p17["realized_gain_rhcp_dbi_anchored"]) - GAIN4_TOL_DB,
        max(p1["realized_gain_rhcp_dbi_anchored"], p17["realized_gain_rhcp_dbi_anchored"]) + GAIN4_TOL_DB,
    )
    gdb = a["peak_rhcp_gain_dbi"]
    pl_theta_phi = p1["check_item_neighbour_lhcp_lobe"]["theta_phi_deg"]
    near = _angle_to(tgrid, pgrid, *pl_theta_phi) <= 6.0
    at_pred = float(10 * np.log10(np.nanmax(np.where(near & np.isfinite(pl), pl, 0.0)) / peak))
    pred_lvl = (
        p1["check_item_neighbour_lhcp_lobe"]["rel_to_main_rhcp_peak_db"],
        p17["check_item_neighbour_lhcp_lobe"]["rel_to_main_rhcp_peak_db"],
    )
    rep["checks"] = {
        "peak_gain": {
            "hfss_dbi": gdb,
            "predicted_range_dbi": [lo, hi],
            "state": "GAIN_WITHIN_PREDICTED_RANGE" if lo <= gdb <= hi else "GAIN_OUTSIDE_PREDICTED_RANGE",
        },
        "pointing_error_deg": a["pointing_error_deg"],
        "xpd_db": {"hfss": a["xpd_at_peak_db"], "model_s1": p1["xpd_at_peak_db"]},
        "ar_peak_db": {
            "hfss": a["ar_at_peak_db"],
            "model_s1": p1["ar_at_peak_db"],
            "worst_in_beam_hfss": a["ar_max_over_minus3db_beam_db"],
            "worst_in_beam_model_s1": p1["ar_worst_in_minus3db_beam_db"],
        },
        "neighbour_lhcp_lobe": {
            "hfss_at_predicted_direction_db": at_pred,
            "predicted_direction_theta_phi_deg": pl_theta_phi,
            "predicted_db_s1_s1p7": list(pred_lvl),
            "hfss_strongest_lhcp_outside_9deg": lobe,
            "state": "LHCP_LOBE_AT_OR_BELOW_PREDICTION"
            if at_pred <= max(pred_lvl) + 1.0
            else "LHCP_LOBE_ABOVE_PREDICTION",
        },
    }
    return rep


def coupling_report(s4p, pred):
    from ..phase9.hfss_reference import read_touchstone

    net = read_touchstone(s4p)
    f = net.frequency_hz / 1e9
    i = int(np.argmin(np.abs(f - F_DESIGN)))
    db = lambda a, b: float(20 * np.log10(max(abs(net.component(a, b)[i]), 1e-12)))  # noqa: E731
    out = {"f_nearest_ghz": float(f[i]), "s_ii_db": {str(p): db(p, p) for p in range(1, 5)}, "pairs": {}}
    est = pred.get("coupling", {}).get("adjacent_pairs", {})
    for a, b in PAIRS:
        v = {
            "s_ij_db": db(a, b),
            "s_ji_db": db(b, a),
            "worst_over_sweep_db": float(np.max(20 * np.log10(np.maximum(np.abs(net.component(a, b)), 1e-12)))),
        }
        e = est.get(f"{a}-{b}")
        if e:
            v["raw_estimate_db"], v["calibrated_estimate_db"] = e["raw_estimate_db"], e.get("calibrated_estimate_db")
            v["state"] = (
                "COUPLING_BELOW_RAW_ESTIMATE" if v["s_ij_db"] <= e["raw_estimate_db"] else "COUPLING_ABOVE_RAW_ESTIMATE"
            )
        else:
            v["state"] = "NO_PREDICTION"
        out["pairs"][f"{a}-{b}"] = v
    out["not_adjacent_2-4_db"] = db(2, 4)
    return out


def compare_fourhex(args):
    pred = json.loads(Path(args.pred).read_text())
    rep = {"pred": str(args.pred), "ports": {}, "coupling": "NOT_PROVIDED"}
    gains, params = _kv(args.gain), _kv(args.params)
    for port, path in gains.items():
        rep["ports"][str(port)] = port_report(port, path, params.get(port), pred)
    if args.s4p:
        rep["coupling"] = coupling_report(args.s4p, pred)
    return rep


def markdown_fourhex(rep):
    L = ["# hex4_tensor_cp_static: HFSS vs prediction", ""]
    for p, r in rep["ports"].items():
        c = r["checks"]
        if "peak_gain" not in c:
            L.append(f"* port {p}: {c['state']}")
            continue
        a = r["measured"]
        L += [
            f"## Port {p} alone (beam 30, {PHI4[int(p)]:g})",
            "",
            f"* peak RHCP {c['peak_gain']['hfss_dbi']:.2f} dBi vs predicted {c['peak_gain']['predicted_range_dbi'][0]:.1f}..{c['peak_gain']['predicted_range_dbi'][1]:.1f}: **{c['peak_gain']['state']}**; "
            f"pointing error {a['pointing_error_deg']:.1f} deg",
            f"* XPD {c['xpd_db']['hfss']:.1f} dB (model {c['xpd_db']['model_s1']:.1f}); AR {c['ar_peak_db']['hfss']:.2f} dB at peak / {c['ar_peak_db']['worst_in_beam_hfss']:.2f} worst in beam "
            f"(model {c['ar_peak_db']['model_s1']:.2f} / {c['ar_peak_db']['worst_in_beam_model_s1']:.2f}); RHCP SLL {a['rhcp_sll_db']:.1f} dB",
            f"* neighbour LHCP lobe: HFSS {c['neighbour_lhcp_lobe']['hfss_at_predicted_direction_db']:.1f} dB at the predicted direction "
            f"{c['neighbour_lhcp_lobe']['predicted_direction_theta_phi_deg']} (predicted {c['neighbour_lhcp_lobe']['predicted_db_s1_s1p7'][0]:.1f} / {c['neighbour_lhcp_lobe']['predicted_db_s1_s1p7'][1]:.1f} dB): "
            f"**{c['neighbour_lhcp_lobe']['state']}**; strongest LHCP outside 9 deg: {r['strongest_lhcp_lobe_outside_9deg']['level_rel_rhcp_peak_db']:.1f} dB at "
            f"{r['strongest_lhcp_lobe_outside_9deg']['theta_phi_deg']}",
            f"* budget: main beam <= 8 deg {a['budget_of_radiated']['main_beam_le_8deg_total']:.2f}, grazing {a['budget_of_radiated']['grazing_70_to_90']:.2f}, "
            f"back {a['budget_of_radiated']['back_hemisphere_gt_90']:.2f}, LHCP {a['budget_of_radiated']['lhcp_total_all_angles']:.2f}",
            "",
        ]
    cp_ = rep["coupling"]
    if cp_ == "NOT_PROVIDED":
        L.append("Coupling: NOT_PROVIDED")
    else:
        L += [
            f"## Coupling at {cp_['f_nearest_ghz']:.2f} GHz",
            "",
            "S_ii: " + ", ".join(f"{k}: {v:.1f} dB" for k, v in cp_["s_ii_db"].items()),
            "",
            "| pair | S_ij dB | worst over sweep | raw estimate | calibrated | state |",
            "|---|---|---|---|---|---|",
        ]
        for k, v in cp_["pairs"].items():
            L.append(
                f"| {k} | {v['s_ij_db']:.1f} | {v['worst_over_sweep_db']:.1f} | {v.get('raw_estimate_db', float('nan')):.1f} | {v.get('calibrated_estimate_db') or float('nan'):.1f} | {v['state']} |"
            )
        L.append(f"\n2-4 (not adjacent): {cp_['not_adjacent_2-4_db']:.1f} dB")
    return "\n".join(L) + "\n"


def main_fourhex(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--s4p")
    ap.add_argument("--gain", action="append")
    ap.add_argument("--params", action="append")
    ap.add_argument("--pred", default=str(PRED4))
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    rep = compare_fourhex(args)
    txt = markdown_fourhex(rep)
    print(txt)
    if args.out:
        Path(args.out).write_text(json.dumps(rep, indent=1, default=float))
        Path(args.out).with_suffix(".md").write_text(txt)
    return rep
