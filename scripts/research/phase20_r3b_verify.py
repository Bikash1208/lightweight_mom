"""Phase 20 R3b: mesh 60 / 90 / 120 verification of the winner cell (F = 2n at every mesh), Richardson correction exactly as in R4 (p = 2 from meshes 90/120, C_inf = C120 + R2 (C120 - C90),
p = 1 variant for the uncertainty), full target table at 10.0-11.0 GHz, theta 0 / 30 deg, phi 0 / 45 / 90 (theta 30), tolerance corners at 10.5 GHz (meshes 60 and 90, nominal correction added),
and the HFSS_CALIBRATED_TRANSFERRED view (run-2 delta added at phi = 0). MODEL_ONLY; the transferred delta is an assumption, not an HFSS number.

    python scripts/research/phase20_r3b_verify.py <final.json> [--workers N] [--smoke]
Run with OPENBLAS_NUM_THREADS=1 per worker; mesh 120 needs ~9 GB per worker.
"""

import json
import math
import os
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

for _k in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase20_r3b_opt as opt  # noqa: E402
from phase20_r4_metrics import point_metrics  # noqa: E402
from provenance import finish, provenance  # noqa: E402

from src.basis_functions import MeshSettings  # noqa: E402
from src.ris.hfss_calibration import DeltaTable  # noqa: E402
from src.ris.polarization import circular_reflection, to_structure_frame  # noqa: E402
from src.ris.state_library import R1_RICHARDSON, R2_RICHARDSON  # noqa: E402
from src.sparameters import compute_s_parameters_oblique_over_loads  # noqa: E402

FREQS = [10.0, 10.25, 10.5, 10.75, 11.0]
POINTS = [(0.0, 0.0), (30.0, 0.0), (30.0, 45.0), (30.0, 90.0)]
CORNER_POINTS = [(0.0, 0.0), (30.0, 45.0)]
OUT = opt.OUT


def solve(task):
    x, f, th, ph, n, kit_names = task
    cell = opt.make_cell(x)
    ks = opt.kits()
    sets = [cell.unit_cell(st, ks[k]).loads for k in kit_names for st in cell.states()]
    base = cell.unit_cell(cell.states()[0], ks["nominal"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = compute_s_parameters_oblique_over_loads(
            base,
            f * 1e9,
            MeshSettings(n, n, 2 * n),
            opt.PL,
            math.radians(th),
            math.radians(ph),
            sets,
            warn_grating=False,
        )
    return np.array([r.s11 for r in res])


def circ(s11, ph, k):
    return np.array([to_structure_frame(circular_reflection(s11[4 * k + i]), math.radians(ph)) for i in range(4)])


def main(final, workers, smoke):
    d = json.loads(Path(final).read_text())
    x = tuple(d["x_mm"])
    meshes = (24, 36, 48) if smoke else (60, 90, 120)
    freqs = FREQS[2:3] if smoke else FREQS
    points = POINTS[:2] if smoke else POINTS
    names = list(opt.kits())
    tasks = [(x, f, th, ph, n, ["nominal"]) for n in sorted(meshes, reverse=True) for f in freqs for th, ph in points]
    corner_tasks = [
        (x, 10.5, th, ph, n, names)
        for n in meshes[:2]
        for th, ph in (CORNER_POINTS if not smoke else CORNER_POINTS[:1])
    ]
    rec = provenance({"final": d, "meshes": meshes, "smoke": smoke}, workers=workers)
    with Pool(workers) as pool:
        res = pool.map(solve, tasks + corner_tasks, chunksize=1)
    by = {(t[1], t[2], t[3], t[4], len(t[5])): r for t, r in zip(tasks + corner_tasks, res, strict=True)}
    delta = DeltaTable.load(OUT / "hfss_delta_table_run2.npz")
    rows = []
    n60, n90, n120 = meshes
    for f in freqs:
        for th, ph in points:
            c60, c90, c120 = (circ(by[(f, th, ph, n, 1)], ph, 0) for n in (n60, n90, n120))
            c_inf2 = c120 + R2_RICHARDSON * (c120 - c90)
            c_inf1 = c120 + R1_RICHARDSON * (c120 - c90)
            row = {
                "f_ghz": f,
                "theta_deg": th,
                "phi_deg": ph,
                "raw60": point_metrics(c60),
                "mesh120": point_metrics(c120),
                "corrected": point_metrics(c_inf2),
                "corrected_p1": point_metrics(c_inf1),
            }
            row["richardson_uncertainty_max_abs"] = float(np.abs(c_inf1 - c_inf2).max())
            row["correction_max_abs"] = float(np.abs(c_inf2 - c60).max())
            if ph == 0.0:
                cal = np.array([delta.apply(c_inf2[s], f, th, s + 1)[0] for s in range(4)])
                row["transferred_hfss_calibrated"] = point_metrics(cal)
                row["transferred_rr_phase_deg"] = np.degrees(np.angle(cal[:, 0, 0])).tolist()
            rows.append(row)
    corners = []
    for th, ph in CORNER_POINTS if not smoke else CORNER_POINTS[:1]:
        c120n = circ(by[(10.5, th, ph, n120, 1)], ph, 0)
        c90n = circ(by[(10.5, th, ph, n90, 1)], ph, 0)
        corr = (c120n + R2_RICHARDSON * (c120n - c90n)) - circ(
            by[(10.5, th, ph, n60, 1)], ph, 0
        )  # nominal correction (state-wise), added to every corner
        for k, nm in enumerate(names):
            c = circ(by[(10.5, th, ph, n60, len(names))], ph, k) + corr
            corners.append({"corner": nm, "theta_deg": th, "phi_deg": ph, **point_metrics(c)})
    out = {
        "provenance": finish(rec),
        "labels": [
            "MODEL_ONLY",
            "HFSS_CALIBRATED_TRANSFERRED(delta of the baseline cell, phi = 0 only)",
            "PROVISIONAL_DATASHEET",
        ],
        "x_mm": list(x),
        "names": opt.NAMES,
        "points": rows,
        "corners": corners,
    }
    tag = "_smoke" if smoke else ""
    (OUT / f"verification{tag}.json").write_text(json.dumps(out, indent=1))
    lines = [
        "| f | theta | phi | dev | step | min amp | spread | opp | (transferred: step / opp) |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        m = r["corrected"]
        t = r.get("transferred_hfss_calibrated")
        lines.append(
            f"| {r['f_ghz']:g} | {r['theta_deg']:g} | {r['phi_deg']:g} | {m['dev']:.1f} | {m['step']:.1f} | {m['amp_min']:.2f} | {m['spread']:.2f} | {m['opp']:.1f} | "
            + (f"{t['step']:.1f} / {t['opp']:.1f}" if t else "-")
            + " |"
        )
    lines += [
        "",
        "| corner (10.5 GHz) | theta | phi | dev | step | min amp | spread | opp |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for c in corners:
        lines.append(
            f"| {c['corner']} | {c['theta_deg']:g} | {c['phi_deg']:g} | {c['dev']:.1f} | {c['step']:.1f} | {c['amp_min']:.2f} | {c['spread']:.2f} | {c['opp']:.1f} |"
        )
    (OUT / f"VERIFICATION{tag}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    wk = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 8
    main(sys.argv[1], wk, "--smoke" in sys.argv)
