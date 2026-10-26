"""Phase 20 R3 item 4: summary statistics of the locality-test JSONs (`phase20_r3_locality.py`).

    python scripts/research/phase20_r3_locality_summary.py <results/phase20_r3/locality>

Per case: the strongest non-specular order by PREDICTION (the 'anomalous order'): level error [dB] and phase error [deg]; the gate (within 1 dB and 15 deg) result;
the error distribution over all non-specular orders predicted above -25 dB (handedness-preserving r_RR); and where the opposite-handed power goes (fraction of
sum |r_LR|^2 in the specular order (0, 0), measured and predicted).
"""

import json
import sys
from pathlib import Path

import numpy as np

GATE_DB, GATE_DEG = 1.0, 15.0


def wrap(x):
    return (x + 180.0) % 360.0 - 180.0


def summarise_run(run, floor_db=-25.0):
    orders = run["orders"]
    non_spec = [o for o in orders if o["order"] != [0, 0]]
    strongest = max(non_spec, key=lambda o: o["pred_rr_db"])
    sel = [o for o in non_spec if o["pred_rr_db"] > floor_db]
    err_db = [o["err_db"] for o in sel]
    err_deg = [abs(wrap(o["err_deg"])) for o in sel if o["err_deg"] is not None]
    p_meas = np.array([10 ** (o["meas_lr_db"] / 10) for o in orders])
    p_pred = np.array([10 ** (o["pred_lr_db"] / 10) for o in orders])
    i0 = next(i for i, o in enumerate(orders) if o["order"] == [0, 0])
    return {
        "codes": run["codes"],
        "strongest_anomalous_order": strongest["order"],
        "strongest_pred_db": strongest["pred_rr_db"],
        "strongest_meas_db": strongest["meas_rr_db"],
        "strongest_err_db": strongest["err_db"],
        "strongest_err_deg": None if strongest["err_deg"] is None else wrap(strongest["err_deg"]),
        "gate_pass": bool(
            abs(strongest["err_db"]) <= GATE_DB
            and strongest["err_deg"] is not None
            and abs(wrap(strongest["err_deg"])) <= GATE_DEG
        ),
        "n_orders_above_floor": len(sel),
        "abs_err_db_median": float(np.median(np.abs(err_db))) if err_db else None,
        "abs_err_db_max": float(np.max(np.abs(err_db))) if err_db else None,
        "abs_err_deg_median": float(np.median(err_deg)) if err_deg else None,
        "abs_err_deg_max": float(np.max(err_deg)) if err_deg else None,
        "opposite_handed_fraction_in_specular_meas": float(p_meas[i0] / p_meas.sum()),
        "opposite_handed_fraction_in_specular_pred": float(p_pred[i0] / p_pred.sum()),
        "opposite_handed_total_db_meas": float(10 * np.log10(p_meas.sum())),
    }


def main(folder):
    folder = Path(folder)
    out = {}
    for name in ("gradient", "random4x4"):
        f = folder / f"{name}.json"
        if not f.exists():
            continue
        d = json.loads(f.read_text())
        out[name] = {th: [summarise_run(r) for r in runs] for th, runs in d["results"].items()}
    (folder / "summary.json").write_text(json.dumps(out, indent=1))
    for name, by_th in out.items():
        for th, runs in by_th.items():
            for r in runs:
                print(
                    name,
                    th,
                    r["codes"],
                    f"strongest {r['strongest_anomalous_order']} pred {r['strongest_pred_db']:.1f} meas {r['strongest_meas_db']:.1f}"
                    f" err {r['strongest_err_db']:+.2f} dB {r['strongest_err_deg']:+.1f} deg gate={r['gate_pass']}"
                    f" | orders>-25 dB: {r['n_orders_above_floor']} median|err| {r['abs_err_db_median']} dB max {r['abs_err_db_max']} dB"
                    f" | opp-handed in specular: meas {r['opposite_handed_fraction_in_specular_meas']:.2f}"
                    f" pred {r['opposite_handed_fraction_in_specular_pred']:.2f}",
                )


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "results/phase20_r3/locality")
