import sys
from pathlib import Path
import gate_spec as g
rows = Path(__file__).resolve().parents[2] / "results" / "sub6_gate" / "rows"
print(" ".join(r["run_id"] for r in g.cases() if not (rows / (r["run_id"] + ".json")).exists()))
