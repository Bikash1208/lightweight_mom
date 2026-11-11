import json, re, math
from pathlib import Path
rows = Path(__file__).resolve().parents[2] / "results" / "sub6_gate" / "rows"
for p in sorted(rows.glob("*.json")):
    r = json.loads(p.read_text())
    t = r["mesh_raw"]["convergence_text"] or ""
    passes = re.search(r"Completed : (\d+)", t); conv = re.search(r"Converged : (\w+)", t)
    cur = re.search(r"Current   : ([\d.eE+-]+)", t)
    tab = re.findall(r"^\s*(\d+)\|\s*(\d+)\|\s*([\w./]+)\|", t, re.M)
    d = [x for x in r["per_frequency"] if abs(x["f_ghz"] - 5.8) < 1e-6][0]
    g = d.get("port_gamma_1_per_m")
    print(r["run_id"], "passes", passes and passes.group(1), "conv", conv and conv.group(1), "dS", cur and cur.group(1),
          "tets", tab[-1][1] if tab else None, "rb", r["readback_ok"], "gamma", g and [round(v, 4) for v in g],
          "t", round(r["solve_seconds"]))
