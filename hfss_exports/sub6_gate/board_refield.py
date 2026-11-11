"""Re-extract phasor field metrics from a SAVED, already-solved board project (no re-solve)."""
import json, sys, glob
from pathlib import Path
import board_hfss as BH
import gate_hfss as GH

run, fs, port = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
import re
cands = [x for x in glob.glob(str(GH.WORK / "projects" / f"{run}_*_a1")) if re.fullmatch(rf"{re.escape(run)}_\d{{10}}_a1", Path(x).name)]
d = sorted(cands)[-1]
proj = Path(d) / f"{run}.aedt"
row = BH.all_rows()[run]
h = GH.connect(proj, port)
try:
    if len(sys.argv) > 4 and sys.argv[4] == "flux":
        res = BH.flux_metrics(h, fs, f"{BH.SETUP} : LastAdaptive")
        res.update(run_id=run, solve_ghz=fs, project=str(proj.relative_to(GH.WORK)))
        GH.write_json(BH.OUT / "refield" / f"{run}_{fs:g}GHz_flux.json", res)
    else:
        BH.make_lines(h, row)
        res = BH.field_metrics(h, fs, f"{BH.SETUP} : LastAdaptive")
        res.update(run_id=run, solve_ghz=fs, project=str(proj.relative_to(GH.WORK)))
        GH.write_json(BH.OUT / "refield" / f"{run}_{fs:g}GHz.json", res)
    print("REFIELD_DONE", run, flush=True)
finally:
    try:
        h.release_desktop(close_projects=True, close_desktop=True)
    except Exception:
        pass
