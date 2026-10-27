"""Phase 28 item 3: how much package inductance the diode-switched cell tolerates. The MADP-000907-14020 model has no package inductance in the repo (L_s NOT_AVAILABLE_FROM_CURRENT_DATA,
zero in the nominal model); a diode ON state is a series R_s + j w L_s, which at 10 GHz is j6.3 ohm per 0.1 nH. This evaluates a candidate with ASSUMED package inductances
(a sensitivity, not device data) and, for reference, with the documented SMP1345-079LF model (L_s 0.7 nH).

    python scripts/research/phase28_cell_ls_threshold.py <spec.json> <out.json> [--mesh 40] [--workers 16]
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import json  # noqa: E402
import sys  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase28_cell_characterise as ch  # noqa: E402

LS_NH = (0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7)


def main(spec_path, out_path, mesh, workers):
    spec = json.loads(Path(spec_path).read_text(encoding="utf-8"))
    jobs = [(spec, s, mesh, ch.F0, {"madp_Ls_nH": ls} if ls > 0 else None, 1) for ls in LS_NH for s in (0, 1)]
    jobs += [(dict(spec, kit="smp"), s, mesh, ch.F0, None, 1) for s in (0, 1)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(ch.evaluate, jobs))
    out = []
    for job, r in zip(jobs, rows, strict=True):
        out.append(
            {
                "kit": job[0]["kit"],
                "madp_Ls_nH": 0.0 if job[4] is None and job[0]["kit"] == "madp" else (job[4] or {}).get("madp_Ls_nH"),
                "state": job[1],
                **{
                    k: r[k]
                    for k in (
                        "b_bar_ms",
                        "delta_ms",
                        "dx_aniso_ohm",
                        "dx_aniso_at_design_loading_ohm",
                        "alpha_beam_np_per_m",
                    )
                },
            }
        )
    doc = {
        "label": "PHASE28_ITEM3_PACKAGE_INDUCTANCE_SENSITIVITY",
        "labels": ["ASSUMED_SENSITIVITY_NOT_DEVICE_DATA", "MODEL_ONLY"],
        "spec": spec,
        "mesh": mesh,
        "rows": out,
    }
    Path(out_path).write_text(json.dumps(doc, indent=1), encoding="utf-8")
    for r in out:
        print(
            r["kit"],
            r["madp_Ls_nH"],
            r["state"],
            round(r["delta_ms"], 2),
            round(r["dx_aniso_ohm"] or 0, 1),
            round(r["alpha_beam_np_per_m"] or 0, 3),
        )


if __name__ == "__main__":
    a = sys.argv[1:]
    main(
        a[0],
        a[1],
        int(a[a.index("--mesh") + 1]) if "--mesh" in a else 40,
        int(a[a.index("--workers") + 1]) if "--workers" in a else 16,
    )
