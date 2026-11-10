"""9-cell sanity (model prediction, MODEL_ONLY): pattern of the 9-cell test board (BOARD_*_9cell.step cell set = 9 switched cells nearest the feed + static cells within 17 mm of them) from the full-board model currents,
OFF and ON, identity and rot180 (the HFSS no-FR4 full board matches rot180).  Reports, per state, the RHCP peak (theta, phi) over theta 6..60 and the RHCP level (dBi, theta <= 90 normalisation) at (25, 0) and (25, 180).
    SUB6_MESH=96 python jobs1_nine_cell_sanity.py <out.json>"""
import csv, json, math, os, sys
from pathlib import Path
os.environ.setdefault("SUB6_MESH", "96")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1]))
import numpy as np
import diag_hfss_vs_model as G
import sub6_final_design as FD
import sub6_board_model4 as R4
ROOT = HERE.parents[1]
cells = list(csv.DictReader(open(ROOT / "results/sub6_final_board/mesh96/fab/cell_table.csv")))
sw = sorted([c for c in cells if c["class"] == "SWITCHED"], key=lambda c: math.hypot(float(c["x_mm"]), float(c["y_mm"])))[:9]
keep = {int(c["site"]) for c in sw}
for c in cells:
    if c["class"] != "SWITCHED" and any(math.hypot(float(c["x_mm"]) - float(s["x_mm"]), float(c["y_mm"]) - float(s["y_mm"])) <= 17.0 for s in sw):
        keep.add(int(c["site"]))
fd = FD.final_design("cal"); D, sa, ph, S = fd["D"], fd["sa"], fd["ph"], fd["D"]["S"]
m = np.zeros(S.n, bool); m[list(keep)] = True
fl = fd["flips"].astype(float)
res = dict(n_cells=int(m.sum()), switched_sites=sorted(int(c["site"]) for c in sw))
class Sub:  # restrict the model field to the 9-cell cell set
    pass
sub = Sub(); sub.x, sub.y, sub.green = S.x[m], S.y[m], S.green
for st in ("off", "on"):
    jx, jy, _, _ = R4.currents(D, sa, ph, fd["direction"], fl, st)
    et, ep = G.model_field(sub, jx[m], jy[m]); e2 = (G.pad_to_180(et), G.pad_to_180(ep))
    for name, (a, b) in (("identity", e2), ("rot180", G.transform(*e2, "rot", 180))):
        rh = G.circ(a, b)[0]; pt = G.prad(a, b, 90); d = G.dir_db(rh, pt)
        pk = G.peaks(rh, 6, 60)
        res[f"{st}_{name}"] = dict(peaks=pk, D25_0=G.at(d, 25, 0), D25_180=G.at(d, 25, 180), Dmax=float(d[:31].max()))
        print(st, name, pk, "D(25,0) %.1f D(25,180) %.1f Dmax %.1f" % (res[f"{st}_{name}"]["D25_0"], res[f"{st}_{name}"]["D25_180"], res[f"{st}_{name}"]["Dmax"]))
json.dump(res, open(sys.argv[1], "w"), indent=1, default=str)
