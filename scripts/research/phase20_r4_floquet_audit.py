"""Phase 20 R4 Floquet-cutoff audit: differences between the HFSS-package `mom_reference.csv` files at the minimum Floquet cutoff (`*_mincutoff.csv`) and at F = 2n,
and the gradient locality cases at both cutoffs. Writes results/phase20_r4/floquet_audit.json.

    python scripts/research/phase20_r4_floquet_audit.py <gradient_f2n.json>
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from phase20_r3_locality_summary import summarise_run  # noqa: E402


def wrap(x):
    return (x + 180.0) % 360.0 - 180.0


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def package_delta(
    folder,
    new="mom_reference.csv",
    old="mom_reference_mincutoff.csv",
    key_cols=None,
    mag="mom_S11_mag_dB",
    ph="mom_S11_phase_deg",
):
    a, b = read(ROOT / folder / old), read(ROOT / folder / new)
    assert len(a) == len(b), (folder, len(a), len(b))
    dp, dm = [], []
    for ra, rb in zip(a, b, strict=True):
        for k in key_cols or []:
            assert ra[k] == rb[k], (folder, k, ra[k], rb[k])
        dp.append(abs(wrap(float(rb[ph]) - float(ra[ph]))))
        dm.append(abs(float(rb[mag]) - float(ra[mag])))
    dp, dm = np.array(dp), np.array(dm)
    q = lambda v: {"median": float(np.median(v)), "p95": float(np.percentile(v, 95)), "max": float(v.max())}  # noqa: E731
    return {"n_rows": len(a), "phase_change_deg": q(dp), "level_change_dB": q(dm)}


def main(grad_new):
    out = {
        "R1": package_delta("hfss_exports/phase20_r1_loaded_cell", key_cols=["load", "theta_deg", "floquet_mode"]),
        "R2": package_delta(
            "hfss_exports/phase20_r2_diode_cell", key_cols=["state", "model", "f_ghz", "theta_deg", "floquet_mode"]
        ),
        "R3": package_delta(
            "hfss_exports/phase20_r3_cp_cell",
            key_cols=["state", "f_ghz", "theta_deg", "quantity"],
            mag="mom_mag_dB",
            ph="mom_phase_deg",
        ),
        "R3_tolerance_10GHz": package_delta(
            "hfss_exports/phase20_r3_cp_cell",
            new="mom_tolerance_10GHz.csv",
            old="mom_tolerance_10GHz_mincutoff.csv",
            key_cols=["corner", "state", "theta_deg"],
            mag="r_RR_mag_dB",
            ph="r_RR_phase_deg",
        ),
    }
    old = json.loads((ROOT / "results/phase20_r3/locality/gradient.json").read_text())["results"]
    new = json.loads(Path(grad_new).read_text())["results"]
    out["locality_gradient"] = {
        th: [
            {"codes": a["codes"], "min_cutoff": summarise_run(a), "F_2n": summarise_run(b)}
            for a, b in zip(old[th], new[th], strict=True)
        ]
        for th in old
    }
    (ROOT / "results/phase20_r4/floquet_audit.json").write_text(json.dumps(out, indent=1))
    for k, v in out.items():
        if k != "locality_gradient":
            print(k, json.dumps(v))
    for th, runs in out["locality_gradient"].items():
        for r in runs:
            for tag in ("min_cutoff", "F_2n"):
                s = r[tag]
                print(
                    th,
                    r["codes"],
                    tag,
                    f"order {s['strongest_anomalous_order']} pred {s['strongest_pred_db']:.1f} meas {s['strongest_meas_db']:.1f} err {s['strongest_err_db']:+.2f} dB "
                    f"{s['strongest_err_deg']:+.1f} deg gate={s['gate_pass']} opp-handed-in-specular {s['opposite_handed_fraction_in_specular_meas']:.2f}",
                )


if __name__ == "__main__":
    main(sys.argv[1])
