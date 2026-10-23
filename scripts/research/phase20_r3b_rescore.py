"""Phase 20 R3b: (1) rebuild the HFSS delta table from the user's run 2 (0.05 GHz steps, residual reference plane corrected) and (2) re-score the finalists of every Nelder-Mead start with it, reporting BOTH rankings
(run-1 delta used during the optimisation vs run-2 delta). The winner on the re-scored objective is written to results/phase20_r3b/cp2_r3b_final.json for the mesh 60/90/120 verification.

    python scripts/research/phase20_r3b_rescore.py [n_top_per_start]

Run-2 plane: the user's note says the port plane moved 1.718 mm while 1.236 mm was needed (residual 0.482 mm). The shift between run 1 and run 2 measured from the data (TE and TM co-polar, |S| > 0.9,
median over the 0.5 GHz common points) is 1.61 mm, so the residual correction used here is d = 1.236 - 1.61 = -0.374 mm (data) and the user's value -0.482 mm is reported as the alternative; a plane error is a COMMON
phase for the four states, so it does not change the ladder steps, amplitudes or the opposite-handed level that the objective scores. HFSS_CALIBRATED (user-approved delta), transferred to the modified cells.
"""

import json
import math
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase20_r3b_opt as opt  # noqa: E402

from src.ris.hfss_calibration import DeltaTable, fit_reference_plane_mm, read_hfss_reim  # noqa: E402
from src.ris.state_library import StateLibrary  # noqa: E402

HF = ROOT / "hfss_results" / "phase20_r3_cp_cell"
C0 = 299792458.0


def run_shift_mm(a, b):
    """Plane shift [mm] of run 2 (`b`) relative to run 1 (`a`): median over TE / TM co-polar entries with |S| > 0.9 at the common frequencies."""
    ds = []
    for k, m1 in a.items():
        if k not in b:
            continue
        for i in (0, 1):
            if abs(m1[i, i]) > 0.9:
                st, f, th, _ = k
                ph = np.angle(b[k][i, i] / m1[i, i])
                ds.append(ph / (2 * (2 * np.pi * f * 1e9 / C0) * math.cos(math.radians(th))) * 1e3)
    return float(np.median(ds)), float(np.percentile(ds, 16)), float(np.percentile(ds, 84))


def build_delta_v2():
    lib = StateLibrary.load(ROOT / "results/phase20_r4/cp2_library_v1.npz")
    d1, _ = fit_reference_plane_mm(
        HF / "hfss_bare_slab_theta0_30.csv", ROOT / "hfss_exports/phase20_r3_cp_cell/bare_grounded_slab_reference.csv"
    )
    r1 = read_hfss_reim(HF / "hfss_cp2_states1-4_theta0_30.csv", 0.0)
    r2raw = read_hfss_reim(HF / "hfss_cp2_states1-4_theta0_30_run2_fine.csv", 0.0)
    shift, lo, hi = run_shift_mm(r1, r2raw)
    d_res = d1 - shift
    r2 = read_hfss_reim(HF / "hfss_cp2_states1-4_theta0_30_run2_fine.csv", d_res)
    src = (
        f"hfss_cp2_states1-4_theta0_30_run2_fine.csv (user HFSS run 2, 0.05 GHz steps); residual plane correction {d_res:.3f} mm = run-1 plane {d1:.3f} - measured run-to-run shift {shift:.3f} mm "
        f"(16-84 % {lo:.3f}..{hi:.3f}; the user's note gives 1.718 mm -> {d1 - 1.718:.3f} mm); cross sign +1; approved by the user 2026-10-23"
    )
    t = DeltaTable.build(lib, r2, src, approved_by_user=True)
    t.save(opt.OUT / "hfss_delta_table_run2.npz")
    return {
        "plane_run1_mm": d1,
        "run_shift_mm": shift,
        "run_shift_16_84": [lo, hi],
        "residual_used_mm": d_res,
        "residual_user_note_mm": d1 - 1.718,
        "source": src,
    }


def main(n_top=5):
    info = build_delta_v2()
    print(json.dumps(info, indent=1))
    cache = opt.cached()
    finalists = {}
    for start in opt.STARTS:
        try:
            traj = opt.trajectory(start)
        except Exception:
            continue
        ranked = sorted((k for k in traj if k in cache and cache[k]["J"] < 1e2), key=lambda k: cache[k]["J"])[:n_top]
        if ranked:
            finalists[start] = ranked
    ev2 = opt.Evaluator(opt.OUT / "hfss_delta_table_run2.npz")
    rows = []
    with Pool(int(__import__("os").environ.get("R3B_POOL", "8"))) as pool:
        for start, xs in finalists.items():
            for rank1, x in enumerate(xs, 1):
                J2, parts2, metrics2 = ev2.evaluate(x, pool)
                rows.append(
                    {
                        "start": start,
                        "rank_run1": rank1,
                        "x_mm": list(x),
                        "J_run1": cache[x]["J"],
                        "J_run2": J2,
                        "parts_run1": cache[x]["parts"],
                        "parts_run2": parts2,
                        "metrics_run2": metrics2,
                    }
                )
    for i, r in enumerate(sorted(rows, key=lambda r: r["J_run1"]), 1):
        r["overall_rank_run1"] = i
    for i, r in enumerate(sorted(rows, key=lambda r: r["J_run2"]), 1):
        r["overall_rank_run2"] = i
    (opt.OUT / "rescore.json").write_text(
        json.dumps(
            {"delta_v2": info, "finalists": rows},
            indent=1,
            default=lambda o: o.item() if hasattr(o, "item") else str(o),
        )
    )
    win = min(rows, key=lambda r: r["J_run2"])
    (opt.OUT / "cp2_r3b_final.json").write_text(
        json.dumps(
            {
                "x_mm": win["x_mm"],
                "names": opt.NAMES,
                "period_mm": 15.0,
                "J_run2": win["J_run2"],
                "J_run1": win["J_run1"],
                "start": win["start"],
                "source": f"Winner of the run-2 re-scored finalists (R3b): start {win['start']}.",
            },
            indent=1,
        )
    )
    print("start rank1 rank2 J_run1 J_run2 x")
    for r in sorted(rows, key=lambda r: r["J_run2"]):
        print(
            f"{r['start']:18s} {r['overall_rank_run1']:3d} {r['overall_rank_run2']:3d} {r['J_run1']:8.3f} {r['J_run2']:8.3f} {np.round(r['x_mm'], 3).tolist()}"
        )


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
