"""Compare two independent REF builds and solves; no external reference."""
import json
from pathlib import Path
import shutil

import cell_param_pyaedt as cp
from opt1_cell import configure

configure()
here = Path(__file__).resolve().parent
work = here / "lab_B_opt1_work"
out = here.parents[1] / "results" / "sub6_cell_param" / "lab_B_opt1" / "verify"
out.mkdir(parents=True, exist_ok=True)
name = "OPT1_VERIFY_REF_OFF_along_90"
solves = []
for label in ("A", "B"):
    folder = work / f"verify_{label}"
    row = json.loads((folder / f"{name}.json").read_text())
    summary = json.loads((folder / "build_summary.json").read_text())
    good, reason = cp.readback_ok(row)
    if summary["status"] != "BUILT" or summary["aedt"] != "2025.2" or not good:
        raise RuntimeError(f"VERIFY {label} build/readback failed: {reason}")
    if len(row["modes"]) != 8 or not (folder / f"{name}.prof").exists():
        raise RuntimeError(f"VERIFY {label} missing modes or HFSS profile")
    for file in (folder / f"{name}.json", folder / f"{name}.prof", folder / "build_summary.json"):
        shutil.copy2(file, out / f"{label}_{file.name}")
    solves.append(row)

rows = []
for a, b in zip(solves[0]["modes"], solves[1]["modes"]):
    df = 100 * abs(b["f_ghz"] / a["f_ghz"] - 1)
    dq = 100 * abs(b["q"] / a["q"] - 1)
    rows.append({"mode": a["mode"], "f_A_ghz": a["f_ghz"], "f_B_ghz": b["f_ghz"],
                 "q_A": a["q"], "q_B": b["q"], "f_diff_pct": df, "q_diff_pct": dq,
                 "status": "PASS" if df <= 0.5 and dq <= 5 else "FAIL"})
result = {"status": "PASS" if all(r["status"] == "PASS" for r in rows) else "FAIL",
          "aedt": "2025.2", "geometry": "REF_L4.348", "state": "OFF",
          "direction": "along", "phase_deg": 90, "independent_fresh_builds": 2,
          "frequency_tolerance_pct": 0.5, "q_tolerance_pct": 5,
          "external_reference_used": False, "rows": rows}
(out / "determinism.json").write_text(json.dumps(result, indent=2) + "\n")
print(f"VERIFY {result['status']}: 8 modes, max f difference {max(r['f_diff_pct'] for r in rows):.6g}%, "
      f"max Q difference {max(r['q_diff_pct'] for r in rows):.6g}%")
if result["status"] != "PASS":
    raise SystemExit(1)
