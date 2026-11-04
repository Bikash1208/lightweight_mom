"""HFSS package: four ON-cell variants (+ OFF twins of the two most promising). STEP only. MODEL_ONLY predictions (calibrated on the HFSS baseline run).
    python scripts/research/sidequest_mushroom_variants_package.py -> hfss_exports/sidequest_mushroom_variants/"""
import csv
import json
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import cadquery as cq

import sidequest_mushroom_predict as pr

OUT = ROOT / "hfss_exports/sidequest_mushroom_variants"
P, H, T, PAD, LID = 6.4, 1.575, 0.035, 0.8, 30.0
TWINS = ("V1", "V2")  # most promising (highest predicted ON f(180))


def box(cx, cy, lx, ly, z0, z1):
    return cq.Workplane("XY").workplane(offset=z0).center(cx, cy).rect(lx, ly).extrude(z1 - z0)


def build(name, state, w, gd, r):
    hh = PAD / 2 + gd
    ring = [box(-(w / 2 + hh) / 2, 0, w / 2 - hh, w, 0, T), box((w / 2 + hh) / 2, 0, w / 2 - hh, w, 0, T),
            box(0, (w / 2 + hh) / 2, 2 * hh, w / 2 - hh, 0, T), box(0, -(w / 2 + hh) / 2, 2 * hh, w / 2 - hh, 0, T)]
    patch = ring[0]
    for q in ring[1:]:
        patch = patch.union(q)
    pad = box(0, 0, PAD, PAD, 0, T)
    via = cq.Workplane("XY").workplane(offset=-H).circle(r).extrude(H)
    sub = box(0, 0, P, P, -H, 0).cut(via)
    gnd = box(0, 0, P, P, -H - T, -H)
    air = box(0, 0, P, P, 0, LID).cut(patch).cut(pad)
    zs = T / 2
    x0, x1, y0, y1 = PAD / 2, hh, -PAD / 2, PAD / 2
    sheet = cq.Workplane("XY").add(cq.Face.makeFromWires(cq.Wire.makePolygon([cq.Vector(x0, y0, zs), cq.Vector(x1, y0, zs), cq.Vector(x1, y1, zs), cq.Vector(x0, y1, zs)], close=True)))
    asm = cq.Assembly()
    for nm, o, c in (("Substrate", sub, (0.25, 0.55, 0.25)), ("Ground", gnd, (0.6, 0.6, 0.65)), ("Via", via, (1, 1, 0)), ("Pad", pad, (0.8, 0.55, 0.1)), ("Patch", patch, (0.9, 0.7, 0.2)), ("Air", air, (0.7, 0.8, 1.0))):
        asm.add(o, name=nm, color=cq.Color(*c, 1))
    asm.add(sheet, name="D_ON" if state == "ON" else "D_OFF", color=cq.Color(0.9, 0.1, 0.1, 1))
    asm.export(str(OUT / f"{name}_{state}.step"), exportType="STEP")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    pred = json.loads((ROOT / "results/sidequest_hybrid/mushroom_variant_predictions.json").read_text(encoding="utf-8"))
    rows, files = [], []
    for k, v in pr.VARS.items():
        for st in ("ON",) + (("OFF",) if k in TWINS else ()):
            build(k, st, v["w"], v["gd"], v["r"])
            files.append(f"{k}_{st}.step")
        p = pred["variants"][k]
        rows.append({"variant": k, "patch_w_mm": v["w"], "patch_gap_mm": round(6.4 - v["w"], 2), "pad_mm": PAD, "diode_gap_mm": v["gd"], "via_radius_mm": v["r"], "files": ", ".join(f for f in files if f.startswith(k)),
                     "ON_f90_M2": p["S0"]["M2"]["ON_90"], "ON_f180_M2": p["S0"]["M2"]["ON_180"], "ON_f90_M1": p["S0"]["M1"]["ON_90"], "ON_f180_M1": p["S0"]["M1"]["ON_180"],
                     "ON_f180_M2_S1": p["S1"]["M2"]["ON_180"], "OFF_beta_M2": (p["OFF_beta_over_k0_10p54"]["M2"] or [None])[0], "OFF_beta_M1": (p["OFF_beta_over_k0_10p54"]["M1"] or [None])[0]})
    with open(OUT / "variants.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    (OUT / "predictions.json").write_text(json.dumps(pred | {"files": files, "twins": list(TWINS), "lid_mm": LID}, indent=1, default=float), encoding="utf-8", newline="\n")
    print(files)
    for r_ in rows:
        print({k: (round(x, 2) if isinstance(x, float) else x) for k, x in r_.items() if k != "files"})
