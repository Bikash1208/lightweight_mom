"""Audit beta extrapolation/fixed-point residuals from existing tensors only."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import sub6_iter1_island_retune as R

groups = {}
for line in (R.OUT/"refine_nodes.jsonl").read_text().splitlines():
    row = json.loads(line)
    groups.setdefault((R.key(row["spec"]),row["mesh"]),[]).append(row)
result = []
for (ident,mesh),raw in groups.items():
    summary = R.summaries(raw)[0]
    for st in ("nominal|OFF","nominal|ON"):
        nodes = {}
        for row in raw:
            y = np.asarray(row["tensors"][st])
            nodes.setdefault(row["kt"],{})[row["direction"]] = {st:y[:,:,0]+1j*y[:,:,1]}
        b = summary["states"][st]["beta_over_k0"]
        yx,yy = R.S.interp_y(nodes,st,b)
        mean = .5*(yx[0,0].imag+yy[1,1].imag)
        b_next = R.S.find_beta_from_sheet_admittance(1j*mean,R.S.make_sub(R.SUB),5.8e9).beta_rad_per_m/R.free_space_wavenumber(5.8e9)
        result.append(dict(id=ident,mesh=mesh,state=st,beta_over_k0=b,
            beta_from_sheet_at_returned_beta_over_k0=b_next,
            fixed_point_absolute_residual=abs(b_next-b),
            fixed_point_relative_residual=abs(b_next-b)/abs(b),
            kt_node_interval=[min(nodes),max(nodes)],
            extrapolated=summary["states"][st]["extrapolated"],
            fixed_point_gate="RESIDUAL_LE_1E_MINUS5" if abs(b_next-b)<=1e-5 else "FIXED_POINT_NOT_CONVERGED"))
out = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED; post-stop analysis of existing tensors only",rows=result)
(R.OUT/"sheet_iteration_audit.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
for row in result:
    print(row["id"],row["mesh"],row["state"],"b",row["beta_over_k0"],"next",row["beta_from_sheet_at_returned_beta_over_k0"],row["fixed_point_gate"])
