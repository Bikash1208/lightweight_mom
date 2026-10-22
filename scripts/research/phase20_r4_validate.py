"""Phase 20 R4 item 4: held-out validation of the library interpolator.

(a) INTERPOLATION error: N random off-grid (f, theta, phi) points (seeded), all four states solved at the grid mesh (60) and compared with the interpolated RAW table
    (`c_raw`): this isolates the interpolation from the mesh correction.
(b) CORRECTED-LIBRARY error: M of those points re-solved at mesh 120; compared with the interpolated CORRECTED table (`c_corr`) and with the raw one: how much of the mesh
    error the correction map removes and what is left.

    python scripts/research/phase20_r4_validate.py a|b [N] [out.json]

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


def draw_points(n, seed=20261021):
    rng = np.random.default_rng(seed)
    pts = []
    while len(pts) < n:  # points with a propagating higher order are outside the library domain and are not drawn
        f, th, ph = float(rng.uniform(8, 12)), float(rng.uniform(0, 45)), float(rng.uniform(0, 360))
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


def summarise(rows):
    out = {}
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


def main(mode, n, out):
    lib = StateLibrary.load(LIB)
    pts = draw_points(n)
    if mode == "a":
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        with Pool(11) as pool:
            solved = pool.map(solve_point, [(f, th, ph, 60, F_REF[60]) for f, th, ph in pts], chunksize=1)
        key = "c_raw"
    else:
        with Pool(3) as pool:  # 3 processes x OMP_NUM_THREADS (set by the caller, e.g. 4): peak ~9 GB each
            solved = pool.map(solve_point, [(f, th, ph, 120, F_REF[120]) for f, th, ph in pts], chunksize=1)
        key = "c_corr"
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
    res = {"mode": mode, "n_points": n, "gate_theta_le_30": gate, "summary": summ}
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
    mode = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else (60 if mode == "a" else 20)
    main(mode, n, sys.argv[3] if len(sys.argv) > 3 else ROOT / "results" / "phase20_r4" / f"validation_{mode}.json")
