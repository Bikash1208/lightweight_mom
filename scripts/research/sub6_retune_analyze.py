"""Read checkpointed retune tensors; no solves or solver changes."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent))
import sub6_iter1_island_retune as R


def read_stage(stage):
    path = R.OUT/(stage+"_nodes.jsonl")
    groups = {}
    if not path.exists():
        return [],[]
    for line in path.read_text().splitlines():
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue  # writer may currently be appending final line
        ident = (R.key(row["spec"]),row["mesh"],row["frequency_hz"],row["floquet_order"])
        groups.setdefault(ident,[]).append(row)
    rows, rejected = [],[]
    for ident,nodes in groups.items():
        if len({(v["kt"],v["direction"]) for v in nodes}) < 4:
            continue
        try:
            rows.extend(R.summaries(nodes))
        except Exception as err:
            rejected.append(dict(id=ident,error_type=type(err).__name__,error=str(err)))
    return rows,rejected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage")
    parser.add_argument("--select",action="store_true")
    args = parser.parse_args()
    rows,rejected = read_stage(args.stage)
    result = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED",rows=rows,rejected=rejected)
    (R.OUT/(args.stage+"_analysis.json")).write_text(json.dumps(result,indent=2),encoding="utf-8")
    ranked = sorted(rows,key=lambda r:r["states"]["nominal|OFF"]["suppression_dB"],reverse=True)
    eligible = [r for r in ranked if r["states"]["nominal|OFF"]["alpha_np_per_m"] >= 2]
    print(f"complete geometries={len(rows)} rejected={len(rejected)} alpha_OFF>=2={len(eligible)}")
    for r in eligible[:10]:
        a,b = r["states"]["nominal|OFF"],r["states"]["nominal|ON"]
        print(r["spec"],r["mesh"],f"alpha {a['alpha_np_per_m']:.6g}/{b['alpha_np_per_m']:.6g} suppression {a['suppression_dB']:.4f} beta {a['beta_over_k0']:.6g}/{b['beta_over_k0']:.6g}")
    if args.stage == "baseline" and rows:
        original = json.loads((R.ROOT/"results/jobs1/diode_mom_bar50.json").read_text())
        reference = original["rows"]["BAR50-02V C0.11 Rp4k (iteration 1 nominal)"]
        diffs = {}
        for state in ("OFF","ON"):
            for field in ("alpha_np_per_m","beta_over_k0"):
                now = rows[0]["states"]["nominal|"+state][field]
                old = reference[state][field]
                diffs[state+"_"+field] = dict(now=now,reference=old,relative_difference=abs(now-old)/abs(old))
        (R.OUT/"baseline_comparison.json").write_text(json.dumps(diffs,indent=2),encoding="utf-8")
        print(json.dumps(diffs,indent=2))
        if max(v["relative_difference"] for v in diffs.values()) > .05:
            raise RuntimeError("HARD_STOP: reproduced baseline differs >5%")
    if args.select:
        if len(eligible) < 3:
            raise RuntimeError("Fewer than three completed candidates with useful OFF leakage")
        (R.OUT/"selected.json").write_text(json.dumps(eligible[:3],indent=2),encoding="utf-8")
        print("Selected top three among completed candidates; frequency/convergence gates remain pending")


if __name__ == "__main__":
    main()
