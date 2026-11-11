"""Write results/sub6_gate/summary_5p8GHz.csv from rows/*.json (no AEDT)."""
import csv
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ROWS = REPO / "results" / "sub6_gate" / "rows"


def at(r, f=5.8):
    return next(d for d in r["per_frequency"] if abs(d["f_ghz"] - f) < 1e-6)


def main():
    rows = {p.stem: json.loads(p.read_text()) for p in sorted(ROWS.glob("*.json"))}
    out = []
    for rid, r in rows.items():
        row = r["row"]
        ref = rows.get(f"REF_H{int(row['lid_mm'])}")
        d = at(r)
        s21, s11 = complex(*d["S21"]), complex(*d["S11"])
        rel = None
        if ref:
            rel = 20 * math.log10(abs(s21) / abs(complex(*at(ref)["S21"])))
        out.append({"run_id": rid, "kind": row["kind"], "n_gate": row["n_gate"], "island_l_mm": row["island_l"] if row["kind"] == "GATE" else "",
                    "state": row["state"] or "", "Rs": row["Rs"] if row["kind"] == "GATE" else "",
                    "Xs": row["Xs"] if row["kind"] == "GATE" else "", "lid_mm": row["lid_mm"],
                    "S21_dB": round(20 * math.log10(abs(s21)), 3),
                    "S21_over_REF_dB": "" if rel is None else round(rel, 3),
                    "abs_S11": round(abs(s11), 4), "absorbed": round(d["absorbed_fraction"], 4),
                    "passes": r.get("passes", ""), })
    p = REPO / "results" / "sub6_gate" / "summary_5p8GHz.csv"
    with p.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(p, len(out))


if __name__ == "__main__":
    main()
