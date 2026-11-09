"""SUB6 FINAL BOARD, 13-pad variant (KiCad + Freerouting route): independent connectivity / DRC gates from the routed tracks, fab tables, RF_UNCHANGED.  MODEL_ONLY / NOT_HFSS_VALIDATED.
    SUB6_MESH=96 SUB6_VARIANT=13pad python sub6_final_13pad_verify.py
Reads  kicad/tracks.json (net, layer, w, a, b), kicad/nets.json, results/sub6_final_board/mesh96/routing.json (RF data, unchanged)
Writes fab/{traces.json, cell_table.csv, drill.csv, drill.drl, string_table.csv, connector_pinout.csv, layout_summary.json, connectivity.json, drc.json, rf_unchanged.json}
"""
import csv
import json
import math
import os
import shutil
import sys
from pathlib import Path

os.environ["SUB6_VARIANT"] = "13pad"
os.environ.setdefault("SUB6_MESH", "96")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_final_layout as L
import sub6_final_13pad as M
import sub6_final_13pad_export as X

KD = L.OUT / "kicad"
FAB = L.OUT / "fab"
SRC = M.SRC
R3 = X.R3


def main():
    FAB.mkdir(parents=True, exist_ok=True)
    tracks = json.load(open(KD / "tracks.json"))
    nets_js = json.load(open(KD / "nets.json"))
    d = json.load(open(SRC / "routing.json"))
    vias = {int(k): (np.array(v["island"]), np.array(v["disc"])) for k, v in d["vias"].items()}
    pads = [np.array([p["x"], p["y"]]) for p in nets_js["pads"]]
    names = {k: p["net"] for k, p in enumerate(nets_js["pads"])}
    traces = [dict(kind="return" if t["net"] == "RETURN" else ("lead" if t["net"].startswith("S") else "link"), net=t["net"], string=0, pos="", w=t["w"], poly=[tuple(t["a"]), tuple(t["b"])]) for t in tracks if abs(t["a"][0] - t["b"][0]) + abs(t["a"][1] - t["b"][1]) > 1e-9]
    res = dict(d=d, vias=vias, pads=pads)
    conn, drc = X.checks(res, traces, names)
    for t in traces:
        pass
    print(conn["state"], conn["fails"][:5], flush=True)
    print(drc["state"], drc["min_clearances_mm"], "min width", drc["min_trace_width_mm"], "min RETURN width", drc["min_RETURN_width_mm"], drc["RETURN_widths_mm"], flush=True)
    # width statistics by net kind
    wk = {}
    for t in traces:
        wk.setdefault(t["kind"], {}).setdefault(str(round(t["w"], 3)), 0.0)
        wk[t["kind"]][str(round(t["w"], 3))] += math.hypot(t["poly"][1][0] - t["poly"][0][0], t["poly"][1][1] - t["poly"][0][1])
    drc["length_mm_by_kind_and_width"] = {k: {w: round(v, 1) for w, v in vv.items()} for k, vv in wk.items()}
    ang = {}
    for t in traces:
        dx, dy = t["poly"][1][0] - t["poly"][0][0], t["poly"][1][1] - t["poly"][0][1]
        a = round(math.degrees(math.atan2(dy, dx)) % 45, 1)
        ang[a if a < 44.95 else 0.0] = ang.get(a if a < 44.95 else 0.0, 0) + 1
    drc["segment_angle_mod_45_deg_counts"] = ang
    # SMA keep-out: square half-size = pad half 3.675 + 1.5 margin = 5.175 mm; trace edge (Chebyshev distance to the origin minus half width) must be >= 5.175
    kmin = 9e9
    for t in traces:
        a, b = t["poly"]
        n = max(2, int(math.hypot(b[0] - a[0], b[1] - a[1]) / 0.05))
        kmin = min(kmin, min(max(abs(a[0] + (b[0] - a[0]) * i / n), abs(a[1] + (b[1] - a[1]) * i / n)) for i in range(n + 1)) - t["w"] / 2)
    drc["SMA_keepout_half_size_mm"] = 5.175
    drc["min_trace_edge_chebyshev_distance_to_origin_mm"] = round(kmin, 3)
    drc["trace_edge_to_SMA_keepout_mm"] = round(kmin - 5.175, 3)
    if kmin < 5.175 - 1e-6:
        drc["state"] = "DRC_FAIL"
        drc.setdefault("violations", {})["SMA_keepout"] = round(kmin - 5.175, 3)
    (FAB / "connectivity.json").write_text(json.dumps(conn, indent=1), encoding="utf-8", newline="\n")
    (FAB / "drc.json").write_text(json.dumps(drc, indent=1), encoding="utf-8", newline="\n")
    (FAB / "traces.json").write_text(json.dumps([dict(net=t["net"], layer="B.Cu", w=t["w"], a=[R3(t["poly"][0][0]), R3(t["poly"][0][1])], b=[R3(t["poly"][1][0]), R3(t["poly"][1][1])]) for t in traces]), encoding="utf-8", newline="\n")
    for fn in ("cell_table.csv", "drill.csv", "drill.drl"):
        shutil.copyfile(SRC / "fab" / fn, FAB / fn)
    # string table / pinout
    pad_of = {}
    for k, p in enumerate(nets_js["pads"]):
        if p["net"].startswith("S"):
            pad_of[int(p["net"][1:])] = k
    srows = []
    for s in d["strings"]:
        sid = s["id"]
        k = pad_of[sid]
        srows.append(dict(string=sid, n_diodes=len(s["cells"]), drive_net=f"S{sid:02d}", drive_pad=nets_js["pads"][k]["net"], drive_pad_index_from_top=k + 1, return_pad="RETURN", return_pad_index_from_top=M.RET_K + 1,
                          link_trace_mm=R3(sum(math.hypot(t["poly"][1][0] - t["poly"][0][0], t["poly"][1][1] - t["poly"][0][1]) for t in traces if t["net"].startswith(f"N{sid:02d}_"))),
                          lead_trace_mm=R3(sum(math.hypot(t["poly"][1][0] - t["poly"][0][0], t["poly"][1][1] - t["poly"][0][1]) for t in traces if t["net"] == f"S{sid:02d}")),
                          V_forward_at_10mA=round(len(s["cells"]) * L.V_DIODE, 2), I_mA=L.I_DIODE * 1e3))
    with open(FAB / "string_table.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(srows[0]))
        w.writeheader()
        w.writerows(srows)
    with open(FAB / "connector_pinout.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pad_index_from_top", "kicad_pad", "net", "string", "x_mm", "y_mm", "voltage_at_10mA_V", "current_mA", "note"])
        for k, p in enumerate(nets_js["pads"]):
            if p["net"] == "RETURN":
                w.writerow([k + 1, k + 1, "RETURN", "all 12 (cathode ends)", R3(p["x"]), R3(p["y"]), "0 (return)", 120, "common return, middle pad"])
            else:
                w.writerow([k + 1, k + 1, p["net"], int(p["net"][1:]), R3(p["x"]), R3(p["y"]), 19.0, 10, "drive (anode end of string)"])
    rf = X.rf_unchanged()
    (FAB / "rf_unchanged.json").write_text(json.dumps(rf, indent=1), encoding="utf-8", newline="\n")
    hexp = d["outline"]
    meta = dict(label="FINAL_CANDIDATE_13PAD_KICAD / MODEL_ONLY / NOT_HFSS_VALIDATED (bias re-route by KiCad + Freerouting; RF unchanged)", gates=dict(connectivity=conn["state"], drc=drc["state"], rf=rf["state"]),
                verification=dict(min_clearances_mm=drc["min_clearances_mm"], requirement_mm=0.15, antipad_requirement_mm=0.172, all_ge_0p15=bool(drc["state"] == "DRC_OK"), antipad_crossings=drc["antipad_crossings"], min_trace_width_mm=drc["min_trace_width_mm"], min_RETURN_width_mm=drc["min_RETURN_width_mm"]),
                strings=srows, n_switched=240, n_vias=480, outline_hex_circumradius_mm=R3(d["r_out"]), connector_x_mm=R3(d["conn"][0]), n_pads=13, return_pad_index_from_top=M.RET_K + 1, pad_names={str(k): p["net"] for k, p in enumerate(nets_js["pads"])},
                note="13-pad bias layout, KiCad + Freerouting (B.Cu only). Diode anode on island, cathode on disc; S<string> drives the anode end of its string, all 12 cathode ends join RETURN.")
    (FAB / "layout_summary.json").write_text(json.dumps(meta, indent=1), encoding="utf-8", newline="\n")
    print("VERIFY_13PAD_DONE", conn["state"], drc["state"], rf["state"], flush=True)


if __name__ == "__main__":
    main()
