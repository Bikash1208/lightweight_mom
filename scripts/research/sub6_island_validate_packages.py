"""Validate existing STEPs without re-exporting or changing geometry."""
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent))
import sub6_island_cell_packages as P

summary={}
for cell in P.CELLS:
    for via in (False,True):
        name=cell+("_vias_fr4_stubs" if via else "_novia")
        folder=P.OUT/name
        path=folder/(name+".step")
        metadata=json.loads((folder/"geometry_checks.json").read_text(encoding="utf-8"))
        assert hashlib.sha256(path.read_bytes()).hexdigest()==metadata["step_sha256"],"STEP changed since exactification audit"
        parts,bottom=P.build(cell,via)
        pre=P.checks(parts,cell,via,bottom)
        post=P.checks(P.read_names(path),cell,via,bottom)
        post["checks"]["EXACT_DIMS"]=metadata["post"]["checks"]["EXACT_DIMS"]
        post["all_pass"]=post["all_pass"] and post["checks"]["EXACT_DIMS"]["ok"]
        metadata.update(pre=pre,post=post)
        (folder/"geometry_checks.json").write_text(json.dumps(metadata,indent=2),encoding="utf-8",newline="\n")
        summary[name]=pre["all_pass"] and post["all_pass"]
        print(name,summary[name],{k:v["ok"] for k,v in post["checks"].items()},flush=True)
(P.OUT/"geometry_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8",newline="\n")
if not all(summary.values()):
    raise RuntimeError("HARD_STOP_GEOMETRY_CHECK_FAILED")
