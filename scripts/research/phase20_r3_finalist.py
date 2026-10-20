"""Phase 20 R3 item 3: finalist correction and full target table of one cell. Meshes 60/90/120 (pitch 0.25/0.167/0.125 mm for a 15 mm cell; 48/72/96
for 12 mm cells), observed-order extrapolation of the nominal co-polar state responses (R2 practice), tolerance corners at the two coarser
meshes corrected by the nominal extrapolation ratio, theta = 0 / 30 (+ 60 reported), band where the steps stay 90 +- 20 deg (coarse mesh, uncorrected).
Writes JSON. MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.

    python scripts/research/phase20_r3_finalist.py <cp2|lp1|lp2> <madp|smp> <params_mm comma sep> <period_mm> <out.json> [--no-fine]
"""

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from phase20_r3_opt import corner_kits, lattice_fit, solve_s11, step_error
from scipy.optimize import brentq

from src.ris.cells import madp_kit, smp1345_kit
from src.ris.polarization import L, R, circular_reflection


def build(family, x, period_mm):
    if family == "cp2":
        from phase20_r3_opt_cp2 import make_cell

        return make_cell(x, period_mm * 1e-3)
    import phase20_r3_opt_lp as lp

    lp.PERIOD = period_mm * 1e-3
    return lp.make_cell(family, x)


def channels(s, family):
    """(co-polar, cross1, cross2) complex amplitudes: CP -> r_RR, r_LR, r_RL; LP -> r_TM,TM, r_TE,TM, r_TM,TE."""
    s = np.asarray(s)
    if family == "cp2":
        c = circular_reflection(s)
        return c[R, R], c[L, R], c[R, L]
    return s[1, 1], s[0, 1], s[1, 0]


def run_set(cell, family, kits, conds, n, pool):
    """{(corner, theta): [(co, x1, x2) per state]} at mesh n."""
    tasks, keys = [], []
    for kname, th in conds:
        for st in cell.states():
            tasks.append((cell, st, kits[kname], 10e9, th, n))
            keys.append((kname, th))
    res = pool.map(solve_s11, tasks)
    out = {}
    for k, s in zip(keys, res, strict=True):
        out.setdefault(k, []).append(channels(s, family))
    return out


def extrapolate3(ns, z):
    """Observed-order extrapolation of three complex values at mesh counts ns (R2 `extrapolate`): returns (z_ext, p, phase uncertainty deg)."""
    d1, d2 = abs(z[0] - z[1]), abs(z[1] - z[2])
    p = None
    if d1 > 1e-9 and d2 > 1e-9:
        try:
            p = float(
                brentq(lambda q, r=d1 / d2: (ns[0] ** -q - ns[1] ** -q) / (ns[1] ** -q - ns[2] ** -q) - r, 0.1, 8.0)
            )
        except ValueError:
            p = None

    def ext(q):
        return (ns[2] ** q * z[2] - ns[1] ** q * z[1]) / (ns[2] ** q - ns[1] ** q)

    if p is None or not (
        0.5 <= p <= 3.0
    ):  # no usable observed order (non-monotone / not asymptotic): first order, flagged by the uncertainty
        p = None
    e = ext(p) if p else ext(1.0)
    if abs(e) > 1.0:  # passive cell: an extrapolated |r| > 1 is an artefact, keep the finest mesh value
        e = z[2]
    unc = max(abs(np.degrees(np.angle(ext(1.0) / e))), abs(np.degrees(np.angle(z[2] / e))))
    return e, p, float(unc)


def metrics_from(co_by_cond, x_by_cond, nominal_perm=None):
    """Target table from co-polar complex lists per (corner, theta) and the max cross amplitude per (corner, theta)."""
    rows = {}
    ph = np.degrees(np.angle(np.array(co_by_cond[("nominal", 0.0)])))
    dev, perm = lattice_fit(ph)
    for (k, th), co in co_by_cond.items():
        co = np.array(co)
        p = np.degrees(np.angle(co))
        a = 20 * np.log10(np.maximum(abs(co), 1e-9))
        d_th, perm_th = lattice_fit(p)
        rows[f"{k}@{th:g}"] = {
            "phase_deg": p.tolist(),
            "amp_db": a.tolist(),
            "amp_spread_db": float(a.max() - a.min()),
            "amp_min_db": float(a.min()),
            "lattice_dev_deg": d_th,
            "step_error_deg_fixed_assignment": step_error(p, perm),
            "xpol_db": float(20 * np.log10(max(x_by_cond[(k, th)], 1e-9))),
        }
    return rows, perm


def main(family, kit_name, x, period_mm, out, fine=True):
    kit = madp_kit() if kit_name == "madp" else smp1345_kit()
    kits = corner_kits(kit)
    cell = build(family, x, period_mm)
    base = round(period_mm / 0.25)  # 60 for 15 mm, 48 for 12 mm
    ns = (base, int(base * 1.5), base * 2)
    nom_conds = [("nominal", 0.0), ("nominal", 30.0), ("nominal", 60.0)]
    corner_conds = [(k, 0.0) for k in kits if k != "nominal"]
    res = {}
    with Pool(8) as pool:
        res[ns[0]] = run_set(cell, family, kits, nom_conds + corner_conds, ns[0], pool)
    with Pool(2) as pool:
        res[ns[1]] = run_set(cell, family, kits, nom_conds + corner_conds, ns[1], pool)
    if fine:
        with Pool(1) as pool:
            res[ns[2]] = run_set(cell, family, kits, nom_conds, ns[2], pool)
    result = {
        "family": family,
        "kit": kit.name,
        "period_mm": period_mm,
        "params_mm": list(x),
        "meshes": list(ns),
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "PROVISIONAL_DATASHEET", "VIA_OMITTED_APPROXIMATION"],
    }
    # corrected nominal
    corr, orders = {}, {}
    for cond in nom_conds:
        for s in range(4):
            if fine:
                z = [res[n][cond][s][0] for n in ns]
                e, p, unc = extrapolate3(ns, z)
                orders[f"{cond[1]:g}/{s}"] = {"p": p, "phase_unc_deg": unc}
                corr.setdefault(cond, []).append(e)
            else:
                corr.setdefault(cond, []).append(res[ns[1]][cond][s][0])
    # corners: mesh ns[1] scaled by the nominal extrapolation ratio (multiplicative on the complex co-pol, per state)
    co = dict(corr)
    for cond in corner_conds:
        co[cond] = [
            res[ns[1]][cond][s][0] * corr[("nominal", 0.0)][s] / res[ns[1]][("nominal", 0.0)][s][0] for s in range(4)
        ]
    xp = {c: max(max(abs(v[1]), abs(v[2])) for v in res[ns[1]][c]) for c in res[ns[1]]}
    result["corrected"], perm = metrics_from(co, xp)
    result["raw"] = {
        str(n): metrics_from(
            {c: [v[0] for v in vals] for c, vals in res[n].items()},
            {c: max(max(abs(v[1]), abs(v[2])) for v in vals) for c, vals in res[n].items()},
        )[0]
        for n in res
    }
    result["perm"] = list(perm)
    result["extrapolation"] = orders
    result["raw_complex"] = {
        str(n): {f"{c[0]}@{c[1]:g}": [[[v.real, v.imag] for v in tup] for tup in vals] for c, vals in r.items()}
        for n, r in res.items()
    }
    Path(out).write_text(json.dumps(result, indent=1))
    print(json.dumps(result["corrected"], indent=0)[:3000])


if __name__ == "__main__":
    a = sys.argv
    main(a[1], a[2], [float(v) for v in a[3].split(",")], float(a[4]), a[5], "--no-fine" not in a)
