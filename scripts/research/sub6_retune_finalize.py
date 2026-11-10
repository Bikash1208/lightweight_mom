"""Package existing retune results/provenance without further MoM runs."""
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/"results/sub6_iter1_island_retune"
report = json.loads((OUT/"mesh_report.json").read_text())
fields = ["id","validation","mesh","slot_l_intent_mm","slot_w_intent_mm",
          "island_l_intent_mm","island_w_intent_mm","physical_gap_mm",
          "mesh_gap_mm","alpha_OFF_np_per_m","alpha_ON_np_per_m",
          "suppression_dB","beta_OFF_over_k0","beta_ON_over_k0",
          "delta_beta_over_k0","delta_beta_rad_per_m",
          "corner_suppression_min_dB","corner_suppression_max_dB",
          "corner_alpha_OFF_min_np_per_m","corner_alpha_OFF_max_np_per_m","mesh_gate"]
with (OUT/"candidate_summary.csv").open("w",encoding="utf-8",newline="") as stream:
    writer = csv.DictWriter(stream,fieldnames=fields,lineterminator="\n")
    writer.writeheader()
    for c in report["candidates"]:
        intent = c["physical_intent_mm"]
        for mesh,data in c["meshes"].items():
            off,on = data["states"]["nominal|OFF"],data["states"]["nominal|ON"]
            writer.writerow(dict(id=c["id"],validation=report["label"],mesh=mesh,
                slot_l_intent_mm=intent["slot_l"],slot_w_intent_mm=intent["slot_w"],
                island_l_intent_mm=intent["island_l"],island_w_intent_mm=intent["island_w"],
                physical_gap_mm=1.,mesh_gap_mm=data["dimensions_mm"]["gap"],
                alpha_OFF_np_per_m=off["alpha_np_per_m"],alpha_ON_np_per_m=on["alpha_np_per_m"],
                suppression_dB=off["suppression_dB"],beta_OFF_over_k0=off["beta_over_k0"],
                beta_ON_over_k0=on["beta_over_k0"],delta_beta_over_k0=off["beta_ON_minus_OFF_over_k0"],
                delta_beta_rad_per_m=on["beta_rad_per_m"]-off["beta_rad_per_m"],
                mesh_gate=c.get("mesh_gate","INCOMPLETE"),
                **{k:data[k] for k in fields if k.startswith("corner_")}))
provenance = json.loads((OUT/"provenance_start.json").read_text())
provenance["end_recorded_utc"] = datetime.now(timezone.utc).isoformat()
provenance["git_status_at_end"] = subprocess.check_output(["git","status","--short"],cwd=ROOT,text=True)
provenance["source_sha256_at_end"] = {str(p.relative_to(ROOT)).replace("\\","/"):
    hashlib.sha256(p.read_bytes()).hexdigest()
    for p in sorted((ROOT/"scripts/research").glob("sub6*retune*.py"))}
one = json.loads((OUT/"search.json").read_text())["rows"]
two = json.loads((OUT/"search2.json").read_text())["rows"]
provenance["geometry_evaluations"] = [len(one),len(two)]
provenance["unique_geometry_count"] = len({r["id"] for r in one+two})
provenance["node_counts"] = {stage:len((OUT/(stage+"_nodes.jsonl")).read_text().splitlines())
    for stage in ("baseline","baseline_wb","baseline_transpose","search","search2","refine")}
provenance["final_status"] = report["status"]
provenance["actual_assembly"] = "existing src.dispersion.leaky_supercell_mom.assemble_z_complex"
provenance["actual_load_solver"] = "exact low-rank update; transpose LU avoids a dense copy"
provenance["worker_cap"] = None
provenance["src_solver_changes"] = False
provenance["frequency_screening"] = "NOT_RUN_AFTER_MANDATED_HARD_STOP"
(OUT/"provenance_end.json").write_text(json.dumps(provenance,indent=2),encoding="utf-8")
# Match committed text bytes before hashing (Windows writers may produce CRLF).
# Trim only log trailing whitespace; preserve the diagnostic text and numbers.
for p in [*OUT.iterdir(),ROOT/"results/retune_initial_fast_suite.log",ROOT/"docs/STATUS.md"]:
    if not p.is_file():
        continue
    original = p.read_bytes()
    normalized = original.replace(b"\r\n",b"\n")
    if p.suffix == ".log":
        normalized = b"\n".join(line.rstrip(b" \t") for line in normalized.split(b"\n"))
    if normalized != original:
        p.write_bytes(normalized)
manifest = dict(label=report["label"],base_commit=provenance["git_commit"],
                producing_source_state="dirty compute branch; source SHA256s in provenance_end.json",
                files=[dict(path=p.name,size_bytes=p.stat().st_size,
                            sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                       for p in sorted(OUT.iterdir()) if p.is_file() and p.name!="manifest.json"])
(OUT/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8",newline="\n")
print("Packaged",len(manifest["files"]),"files; unique geometries",provenance["unique_geometry_count"],report["status"])
