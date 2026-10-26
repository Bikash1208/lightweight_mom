"""Phase 20 R4 full-grid review: old (reduced-grid) vs new (full-grid) library interpolation error, cross-machine node agreement and the mesh-90-only correction residual.

    python scripts/research/phase20_r4_fullgrid_review.py <out.json>

(1) The new library's mesh-60 nodes are direct MoM solves (server): they are held-out truth for the OLD library at every new node that is not an old node; the error of the old
    interpolant of the RAW table (c_raw, i.e. no mesh correction in either side) is reported on r_RR (phase, level), split by theta <= 30 / > 30 and by phi (phi = 90 singled out).
(2) The new library cannot be tested at its own nodes: its interpolation error is bounded by decimation (every second f node / every second theta node / every second phi node
    removed, interpolated back at the removed nodes): an UPPER bound, since the decimated grid is twice as coarse. The brief's held-out gate (>= 50 random off-grid points solved
    directly) is the server job `r4-fullgrid-validate`.
(3) Mac R4 chunks vs server chunks at the common nodes (cross-machine agreement).
(4) Correction map: mesh-90-only Richardson (C90 + 0.8 (C90 - C60)) against the three-mesh one used in the library (C120 + R2 (C120 - C90)), at the 45 correction points.
MODEL_ONLY.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from src.ris.state_library import OutOfDomain, StateLibrary  # noqa: E402

OLD = ROOT / "results/phase20_r4/cp2_library_v1.npz"
NEW = ROOT / "results/phase20_r4/cp2_library_v2_fullgrid.npz"


def err(a, b):
    """r_RR phase [deg] and level [dB] error of a against b."""
    return float(np.degrees(np.angle(a / b))), float(20 * np.log10(abs(a) / abs(b)))


def stats(rows):
    if not rows:
        return None
    ph, db = np.abs([r[0] for r in rows]), np.abs([r[1] for r in rows])
    return {"n": len(rows), "p95_phase_deg": float(np.percentile(ph, 95)), "p95_db": float(np.percentile(db, 95)), "max_phase_deg": float(ph.max()), "max_db": float(db.max()), "median_phase_deg": float(np.median(ph))}


def interp_errors(lib, truth_arrays, skip, key="c_raw"):
    """Errors of `lib` interpolated at the nodes of `truth_arrays` (not in `skip(i, j, k)`), all four states."""
    f, th, ph, c = truth_arrays["f_ghz"], truth_arrays["theta_deg"], truth_arrays["phi_deg"], truth_arrays[key]
    out = {"le30": [], "gt30": []}
    by_phi = {}
    for i, fv in enumerate(f):
        for j, tv in enumerate(th):
            if tv == 0.0:
                continue  # theta = 0 is filled exactly by the frame rotation, not interpolated across phi
            for k, pv in enumerate(ph):
                if skip(i, j, k):
                    continue
                for s in range(4):
                    try:
                        v = lib.structure_circular(fv, tv, pv, s + 1, key=key)[0]
                    except OutOfDomain:
                        continue
                    e = err(v[0, 0], c[i, j, k, s, 0, 0])
                    (out["le30"] if tv <= 30.0 else out["gt30"]).append(e)
                    by_phi.setdefault(float(pv), []).append(e) if tv <= 30.0 else None
    return out, by_phi


def main(out):
    old, new = StateLibrary.load(OLD), StateLibrary.load(NEW)
    na = new.a
    in_old = lambda x, grid: bool(np.any(np.abs(np.asarray(grid) - x) < 1e-9))  # noqa: E731
    skip = lambda i, j, k: in_old(na["f_ghz"][i], old.f) and in_old(na["theta_deg"][j], old.theta) and in_old(na["phi_deg"][k], old.phi)  # noqa: E731
    res = {"labels": ["MODEL_ONLY"]}
    o, by_phi = interp_errors(old, na, skip)
    res["old_library_at_new_nodes"] = {"theta_le_30": stats(o["le30"]), "theta_gt_30": stats(o["gt30"]), "theta_le_30_by_phi": {f"{p:g}": stats(v) for p, v in sorted(by_phi.items())}}

    # decimation of the NEW library
    res["new_library_decimation_upper_bound"] = {}
    for axis, name in ((0, "f_every_2nd"), (1, "theta_every_2nd (0,10,..)"), (2, "phi_every_2nd (30 deg)")):
        keep = [np.arange(0, na["f_ghz"].size, 2), np.arange(0, na["theta_deg"].size, 2), np.arange(0, na["phi_deg"].size, 2)]
        idx = [np.arange(n) for n in (na["f_ghz"].size, na["theta_deg"].size, na["phi_deg"].size)]
        idx[axis] = keep[axis]
        sub = {k: v for k, v in na.items()}
        for key, ax in (("c_raw", 0), ("c_corr", 0)):
            sub[key] = na[key][np.ix_(idx[0], idx[1], idx[2])]
        sub["f_ghz"], sub["theta_deg"], sub["phi_deg"] = na["f_ghz"][idx[0]], na["theta_deg"][idx[1]], na["phi_deg"][idx[2]]
        for key in ("correction", "uncertainty", "s_corr", "s_raw", "higher_order_margin"):
            arr = na[key]
            sub[key] = arr[np.ix_(idx[0], idx[1], idx[2])] if arr.ndim >= 3 else arr
        lib = StateLibrary(sub)
        kept_nodes = [set(np.round(na[n][idx[a]], 9)) for a, n in enumerate(("f_ghz", "theta_deg", "phi_deg"))]
        skip2 = lambda i, j, k: (round(na["f_ghz"][i], 9) in kept_nodes[0] and round(na["theta_deg"][j], 9) in kept_nodes[1] and round(na["phi_deg"][k], 9) in kept_nodes[2])  # noqa: E731
        d, bp = interp_errors(lib, na, skip2)
        res["new_library_decimation_upper_bound"][name] = {"theta_le_30": stats(d["le30"]), "theta_gt_30": stats(d["gt30"]), "theta_le_30_by_phi_90": stats(bp.get(90.0, []))}
    # (3) Mac chunks vs server chunks at common nodes
    chunks = ROOT / "results/phase20_r4/chunks"
    srv = Path("/Volumes/hmsa_sync/data/r4-full-grid/chunks")
    mx = 0.0
    n = 0
    for p in sorted(chunks.glob("grid_f*.npz")):
        q = srv / p.name
        if not q.exists():
            continue
        a, b = np.load(p), np.load(q)
        pa = {tuple(x): i for i, x in enumerate(a["points"])}
        for jb, x in enumerate(b["points"]):
            if tuple(x) in pa:
                ia = pa[tuple(x)]
                mx = max(mx, float(np.abs(a["s"][ia] - b["s"][jb]).max()))
                n += 1
    res["mac_vs_server_common_nodes"] = {"n_points": n, "max_abs_difference_s": mx}
    # (4) correction residual
    import phase20_r4_assemble as asm

    c60, c90, c120 = (asm.load_corr(k) for k in ("corr60", "corr90", "corr120"))
    r2 = (120 / 90) ** 2
    cinf3 = c120 + (c120 - c90) / (r2 - 1.0) if False else None  # the library uses its own R2 (1.2857); see state_library.correction_map
    from src.ris.state_library import correction_map  # noqa: F401

    R2 = 1.2857142857142858
    cinf3 = c120 + R2 * (c120 - c90)
    cinf2 = c90 + (c90 - c60) / ((90 / 60) ** 2 - 1.0)
    d3, d2 = cinf3 - c60, cinf2 - c60
    ph = np.degrees(np.angle((c60 + d2)[..., 0, 0] / (c60 + d3)[..., 0, 0]))
    db = 20 * np.log10(np.abs((c60 + d2)[..., 0, 0]) / np.abs((c60 + d3)[..., 0, 0]))
    res["mesh90_only_correction_residual"] = {
        "definition": "r_RR of (C60 + dC[60,90 Richardson p=2]) vs (C60 + dC[60,90,120]) at the 45 correction points",
        "n": int(ph.size), "max_abs_phase_deg": float(np.abs(ph).max()), "p95_abs_phase_deg": float(np.percentile(np.abs(ph), 95)),
        "max_abs_db": float(np.abs(db).max()), "p95_abs_db": float(np.percentile(np.abs(db), 95)),
        "worst_at_f_ghz_over_points": {f"{f:g}": float(np.abs(ph[i]).max()) for i, f in enumerate(asm.CORR_F)},
    }
    Path(out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1)[:6000])


if __name__ == "__main__":
    main(sys.argv[1])
