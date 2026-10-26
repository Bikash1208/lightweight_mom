"""Phase 20 R4 item 4: held-out validation of the library interpolator.

(a) INTERPOLATION error: N random off-grid (f, theta, phi) points (seeded), all four states solved at the grid mesh (60) and compared with the interpolated RAW table
    (`c_raw`): this isolates the interpolation from the mesh correction.
(b) CORRECTED-LIBRARY error: M of those points re-solved at mesh 120; compared with the interpolated CORRECTED table (`c_corr`) and with the raw one: how much of the mesh
    error the correction map removes and what is left.

    python scripts/research/phase20_r4_validate.py a|b [N] [out.json] [--phi-range LO,HI] [--seed S] [--save-refs refs.npz]
                                                  [--f-range LO,HI] [--use-refs refs.npz]

`PHASE20_R4_LIB` selects the library, `PHASE20_R4_WORKERS` the pool size (default 11 as before). `--phi-range` draws the held-out points only inside that phi window (default 0-360), `--seed` changes the
draw (default the original 20261021), `--save-refs` stores the direct-solve references (points, states, structure-frame circular) so the interpolation can be re-evaluated without new solves.
`--use-refs refs.npz` does exactly that: the points and solved states come from a saved file (no solves; any library via PHASE20_R4_LIB). `--f-range` restricts the drawn frequencies (e.g. 10.4,10.9). The summary also reports the errors inside the resonance band 10.4-10.9 GHz.

Gate (brief): interpolation error <= 3 deg / 0.2 dB at the 95th percentile for theta <= 30 deg; beyond 30 deg reported only. Errors are on r_RR (phase, level) and on the
largest element difference |dC|. MODEL_ONLY.
"""

import json
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from phase20_r4_build_library import solve_point  # noqa: E402

from src.ris.polarization import circular_reflection, to_structure_frame  # noqa: E402
from src.ris.state_library import StateLibrary, higher_order_margin  # noqa: E402

# reference Floquet cutoffs: 2n (the converged choice, same as the library; see FLOQUET_CONVERGENCE.md)
F_REF = {60: 120, 120: 240}
LIB = Path(os.environ.get("PHASE20_R4_LIB", ROOT / "results" / "phase20_r4" / "cp2_library_v1.npz"))


def draw_points(n, seed=20261021, phi_range=(0.0, 360.0), f_range=(8.0, 12.0)):
    rng = np.random.default_rng(seed)
    pts = []
    while len(pts) < n:  # points with a propagating higher order are outside the library domain and are not drawn
        f, th, ph = float(rng.uniform(*f_range)), float(rng.uniform(0, 45)), float(rng.uniform(*phi_range))
        if higher_order_margin(f, th, ph) >= 1.0:
            pts.append((f, th, ph))
    return pts


def to_c(s_states, phi):
    return np.array([to_structure_frame(circular_reflection(s), np.radians(phi)) for s in s_states])


def errors(c_lib, c_ref):
    """Per state: r_RR phase error [deg], level error [dB], max element difference."""
    ph = np.degrees(np.angle(c_lib[:, 0, 0] / c_ref[:, 0, 0]))
    db = 20 * np.log10(np.abs(c_lib[:, 0, 0]) / np.abs(c_ref[:, 0, 0]))
    dmax = np.abs(c_lib - c_ref).reshape(len(c_ref), -1).max(axis=1)
    return ph, db, dmax


BAND = (10.4, 10.9)


def summarise(rows):
    out = {}
    band = [r for r in rows if BAND[0] <= r["f"] <= BAND[1]]
    for tag, sel in (
        ("in_band_10p4_10p9.theta_le_30", [r for r in band if r["theta"] <= 30.0]),
        ("in_band_10p4_10p9.all_theta", band),
    ):
        if sel:
            v = np.abs([x for r in sel for x in r["phase_err_deg"]])
            d = np.abs([x for r in sel for x in r["db_err"]])
            out[f"{tag}.phase_err_deg"] = {
                "n_points": len(sel),
                "median": float(np.median(v)),
                "p95": float(np.percentile(v, 95)),
                "max": float(v.max()),
            }
            out[f"{tag}.db_err"] = {
                "n_points": len(sel),
                "median": float(np.median(d)),
                "p95": float(np.percentile(d, 95)),
                "max": float(d.max()),
            }
    for name, sel in (
        ("theta_le_30", [r for r in rows if r["theta"] <= 30.0]),
        ("theta_gt_30", [r for r in rows if r["theta"] > 30.0]),
    ):
        if not sel:
            continue
        for key in ("phase_err_deg", "db_err", "dmax"):
            v = np.abs([x for r in sel for x in r[key]])
            out[f"{name}.{key}"] = {
                "n": int(v.size),
                "median": float(np.median(v)),
                "p95": float(np.percentile(v, 95)),
                "max": float(v.max()),
            }
    return out


def main(mode, n, out, phi_range=(0.0, 360.0), seed=20261021, save_refs=None, f_range=(8.0, 12.0), use_refs=None):
    lib = StateLibrary.load(LIB)
    if use_refs:
        z = np.load(use_refs)
        pts = [tuple(float(v) for v in p) for p in z["points"]]
        mode, n = str(z["mode"]), len(pts)
    else:
        pts = draw_points(n, seed, phi_range, f_range)
    workers = int(os.environ.get("PHASE20_R4_WORKERS", "11"))
    if use_refs:
        solved = list(z["s"])
        key = "c_raw" if mode == "a" else "c_corr"
    elif mode == "a":
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        with Pool(workers) as pool:
            solved = pool.map(solve_point, [(f, th, ph, 60, F_REF[60]) for f, th, ph in pts], chunksize=1)
        key = "c_raw"
    else:
        with Pool(3) as pool:  # 3 processes x OMP_NUM_THREADS (set by the caller, e.g. 4): peak ~9 GB each
            solved = pool.map(solve_point, [(f, th, ph, 120, F_REF[120]) for f, th, ph in pts], chunksize=1)
        key = "c_corr"
    if save_refs:
        np.savez_compressed(save_refs, points=np.array(pts), s=np.array(solved), mode=mode)
    rows = []
    for (f, th, ph), s in zip(pts, solved, strict=True):
        ref = to_c(s, ph)
        lib_c = np.array([lib.structure_circular(f, th, ph, k + 1, key=key)[0] for k in range(4)])
        p, d, m = errors(lib_c, ref)
        row = {"f": f, "theta": th, "phi": ph, "phase_err_deg": p.tolist(), "db_err": d.tolist(), "dmax": m.tolist()}
        if mode == "b":
            raw = np.array([lib.structure_circular(f, th, ph, k + 1, key="c_raw")[0] for k in range(4)])
            p0, d0, m0 = errors(raw, ref)
            row.update(raw_phase_err_deg=p0.tolist(), raw_db_err=d0.tolist(), raw_dmax=m0.tolist())
        rows.append(row)
    summ = summarise(rows)
    gate = {}
    if "theta_le_30.phase_err_deg" in summ:
        ph95, db95 = summ["theta_le_30.phase_err_deg"]["p95"], summ["theta_le_30.db_err"]["p95"]
        gate = {
            "p95_phase_deg": ph95,
            "p95_db": db95,
            "limit_phase_deg": 3.0,
            "limit_db": 0.2,
            "result": "GATE_PASS" if ph95 <= 3.0 and db95 <= 0.2 else "GATE_FAIL",
        }
    res = {
        "mode": mode,
        "n_points": n,
        "phi_range": list(phi_range),
        "f_range": list(f_range),
        "seed": seed,
        "refs": str(use_refs) if use_refs else None,
        "library": str(LIB),
        "gate_theta_le_30": gate,
        "summary": summ,
    }
    if mode == "b":
        raw_rows = [
            {
                "theta": r["theta"],
                "phase_err_deg": r["raw_phase_err_deg"],
                "db_err": r["raw_db_err"],
                "dmax": r["raw_dmax"],
            }
            for r in rows
        ]
        res["summary_raw_vs_mesh120"] = summarise(raw_rows)
    res["rows"] = rows
    Path(out).write_text(json.dumps(res, indent=1))
    print(json.dumps({"gate_theta_le_30": gate, **summ}, indent=1))


if __name__ == "__main__":
    argv = sys.argv[1:]
    opts = {}
    for key in ("--phi-range", "--seed", "--save-refs", "--f-range", "--use-refs"):
        if key in argv:
            i = argv.index(key)
            opts[key] = argv[i + 1]
            del argv[i : i + 2]
    mode = argv[0]
    n = int(argv[1]) if len(argv) > 1 else (60 if mode == "a" else 20)
    main(
        mode,
        n,
        argv[2] if len(argv) > 2 else ROOT / "results" / "phase20_r4" / f"validation_{mode}.json",
        tuple(float(x) for x in opts["--phi-range"].split(",")) if "--phi-range" in opts else (0.0, 360.0),
        int(opts["--seed"]) if "--seed" in opts else 20261021,
        opts.get("--save-refs"),
        tuple(float(x) for x in opts["--f-range"].split(",")) if "--f-range" in opts else (8.0, 12.0),
        opts.get("--use-refs"),
    )
