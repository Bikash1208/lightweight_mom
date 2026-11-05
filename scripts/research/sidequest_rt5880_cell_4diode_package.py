"""HFSS package for the independent-axes candidate C1: the exported rim-disc cell with BOTH halves of each slit bridged: x diodes at (0, +-0.8) on one control line, y diodes at (+-0.8, 0) on another (4 diodes per cell).
Export standard (hfss_exports/EXPORT_STANDARD.md): METAL_PEC, SUBSTRATE, AIR (ground bottom to lid), sheets DX_TOP, DX_BOTTOM, DY_RIGHT, DY_LEFT; checks before export and on the re-imported STEP. MODEL_ONLY predictions from independent_axes.json.
    python sidequest_rt5880_cell_4diode_package.py   -> hfss_exports/sidequest_rt5880_cell_4diode/"""
import itertools
import json
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import cadquery as cq
import numpy as np

import check_step_cell as CK
import sidequest_rt5880_cell_package as K
from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap

OUT = ROOT / "hfss_exports/sidequest_rt5880_cell_4diode"
POS, WD, G, LS = K.POS, K.WD, K.G, K.LS


def build():
    parts = K.build()
    parts.pop("DX_SHEET"); parts.pop("DY_SHEET")
    parts["DX_TOP"] = K.sheet(0, POS, G, WD)
    parts["DX_BOTTOM"] = K.sheet(0, -POS, G, WD)
    parts["DY_RIGHT"] = K.sheet(POS, 0, WD, G)
    parts["DY_LEFT"] = K.sheet(-POS, 0, WD, G)
    return parts


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    parts = build()
    comp = cq.Compound.makeCompound(list(parts.values()))
    pre = CK.check(comp, pos=POS, ls=LS, g=G, wd=WD, both=True)
    asm = cq.Assembly()
    cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "AIR": (0.7, 0.8, 1.0), "DX_TOP": (0.9, 0.1, 0.1), "DX_BOTTOM": (0.9, 0.1, 0.1), "DY_RIGHT": (0.1, 0.1, 0.9), "DY_LEFT": (0.1, 0.1, 0.9)}
    for nm, shp in parts.items():
        asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
    path = OUT / "cell_geometry_4diode.step"
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    snapped, stats = snap_step_text(raw, K.intended(), description="sidequest RT5880 coded cell, 4 diodes, mm, 0.001 mm grid")
    path.write_text(snapped, encoding="utf-8", newline="\n")
    post = CK.check_file(path, pos=POS, ls=LS, g=G, wd=WD, both=True)
    bad = inexact_literals(path.read_text(encoding="utf-8"), K.intended())
    vs = verify_snap(raw, snapped)
    exact = {"ok": len(bad) == 0 and vs.get("collapsed_edges", 1) == 0 and not vs.get("distinct_merge_violations"), "inexact_literals": len(bad), "n_snapped": stats["n_snapped"]}
    pre["checks"]["EXACT_DIMS"] = exact
    post["checks"]["EXACT_DIMS"] = exact
    pre["PASS"] &= exact["ok"]; post["PASS"] &= exact["ok"]
    print("PRE ", {k: v["ok"] for k, v in pre["checks"].items()}, pre["PASS"])
    print("POST", {k: v["ok"] for k, v in post["checks"].items()}, post["PASS"])
    (OUT / "geometry_checks.json").write_text(json.dumps({"before_export": pre, "re_imported_step": post}, indent=1, default=str), encoding="utf-8", newline="\n")
    if not (pre["PASS"] and post["PASS"]):
        sys.exit("geometry checks failed")
    with open(OUT / "states.csv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("state,x_diodes(DX_TOP+DX_BOTTOM),y_diodes(DY_RIGHT+DY_LEFT),note\n")
        for sx, sy in itertools.product(("OFF", "ON"), repeat=2):
            el = lambda s: "C 28.16 fF" if s == "OFF" else "R 2.44 ohm"  # noqa: E731
            fh.write(f"x{sx}_y{sy},{el(sx)},{el(sy)},R only = ON / C only = OFF; current line across the 0.4 mm gap; x diodes share one control line\n")
    d = json.loads((ROOT / "results/sidequest_hybrid/independent_axes.json").read_text(encoding="utf-8"))
    rows = [r for r in d["rows"] if r["candidate"] == "C1_rim_4diodes"]
    pred = {"label": "MODEL_ONLY (loaded Bloch tensor, mesh 64, lid-aware transverse resonance, +-0.15 GHz); OFF = open (HFSS C/C ~ open) or 28.16 fF load", "rows": []}
    for r in rows:
        for st_, f in r["f_x_run_GHz"].items():
            a, b = st_.split("/")
            pred["rows"].append({"OFF_model": r["OFF"], "state_x/y": st_, "f_x_run_GHz": f, "f_y_run_GHz": r["f_x_run_GHz"][f"{b}/{a}"], "Y_xx_mS_at_f_res": r["Y_at_f_res"][st_]["xx"]})
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1, default=float), encoding="utf-8", newline="\n")
    for r in pred["rows"]:
        print(r)


if __name__ == "__main__":
    main()
