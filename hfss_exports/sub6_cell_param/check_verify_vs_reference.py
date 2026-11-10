"""Compare a lab VERIFY (verify.json) with the reference VERIFY of the main server.

Hard gate: each of the four frequencies within 0.5 % of the reference (user rule 2026-11-10) AND the lab's own VERIFY status PASS
(f <= 1 %, Q <= 30 % against the user's HFSS targets). Q differences against the reference are reported (soft, 10 % flag).

python check_verify_vs_reference.py --lab lab_B_work/verify/verify.json --out lab_B_work/verify/verify_vs_reference.json
Exit code 0 = PASS, 1 = FAIL (do NOT run the queue).
"""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lab", type=Path, required=True)
    ap.add_argument("--ref", type=Path, default=HERE / "verify_reference" / "verify_reference.json")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--f-tol-pct", type=float, default=0.5)
    a = ap.parse_args()
    lab = json.loads(a.lab.read_text()); ref = json.loads(a.ref.read_text())
    rows = []; ok = lab.get("status") == "PASS"
    for r in ref["table"]:
        m = next((x for x in lab["table"] if x["state"] == r["state"] and x["direction"] == r["direction"]), None)
        if m is None:
            ok = False; rows.append({"state": r["state"], "direction": r["direction"], "result": "MISSING"}); continue
        df = 100 * abs(m["f_ghz"] / r["f_ghz"] - 1); dq = 100 * abs(m["q"] / r["q"] - 1)
        row_ok = df <= a.f_tol_pct
        ok = ok and row_ok
        rows.append({"state": r["state"], "direction": r["direction"], "f_ref_ghz": round(r["f_ghz"], 5), "f_lab_ghz": round(m["f_ghz"], 5),
                     "f_diff_pct": round(df, 4), "q_ref": round(r["q"], 1), "q_lab": round(m["q"], 1), "q_diff_pct": round(dq, 2),
                     "q_flag_gt_10pct": dq > 10, "result": "PASS" if row_ok else "FAIL"})
    out = {"status": "PASS" if ok else "FAIL", "lab_verify_status": lab.get("status"), "lab_aedt": lab.get("aedt"), "reference_aedt": ref.get("aedt"),
           "f_tolerance_pct": a.f_tol_pct, "rows": rows}
    if a.out:
        a.out.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=1))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
