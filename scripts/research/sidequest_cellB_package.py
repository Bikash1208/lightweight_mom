"""HFSS unit-cell package, FAMILY B (ON raises Y): RT5880 disc r 3.0 with a cross slit cut THROUGH the disc (ls 3.1: four islands, no rim), the two diodes (both on the line) bridge the slit at 2.0 mm from the centre
(x sheet across the vertical slit at y = +2.0, y sheet across the horizontal slit at x = +2.0, 0.4 x 1.0). OFF = open (no boundary), ON = Impedance boundary Rs = 2.44 * 1.0 / 0.4 = 6.1 ohm/sq.
MoM prediction (MODEL_ONLY): OFF/ON eigen frequencies at 90 deg (transverse resonance on the loaded Bloch tensor, mesh 64, 30 mm PEC lid), with the HFSS-validated family-A cell (rim, ls 2.8) computed in the same run as the pipeline check
(HFSS 2026-11-05: A open 10.670, A R/R 11.006 GHz).     python sidequest_cellB_package.py   -> hfss_exports/sidequest_cellB_through_slit/  (+ sync to C:\\hmsa_sync\\data\\sidequest\\sidequest_cellB_through_slit)"""
import csv
import json
import shutil
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import cadquery as cq
import numpy as np

import sidequest_cell_search as C
import sidequest_k2_default_lib as L
import sidequest_rt5880_cell_package as M
from sidequest_board_k2_lp_package import read_names
from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap

OUT = ROOT / "hfss_exports/sidequest_cellB_through_slit"
SYNC = Path(r"C:\hmsa_sync\data\sidequest\sidequest_cellB_through_slit")
LS_B, POS_B = 3.1, 2.0
RS = 2.44 * 1.0 / 0.4


def _job(a):
    import warnings as w
    w.simplefilter("ignore")
    ls, pos, mode, f = a
    return C.bloch_loaded(L.cell(ls, pos, mode), f, 64, M.BETA, 0.0)


def predictions():
    cases = {"B_through_slit_pos2.0": (LS_B, POS_B), "A_rim_pos0.8 (HFSS-validated cell)": (2.8, 0.8)}
    jobs = [(ls, pos, mode, f) for ls, pos in cases.values() for mode in ("open", "ON") for f in M.FG]
    with ProcessPoolExecutor(12) as ex:
        res = list(ex.map(_job, jobs))
    it = iter(res)
    out = {"label": "MODEL_ONLY (loaded Bloch tensor, mesh 64, transverse resonance with the 30 mm PEC lid, +-0.15 GHz); OFF = open slit, ON = 2.44 ohm only", "bare_substrate_GHz_with_lid": M.f_bare(True), "rows": []}
    for name, (ls, pos) in cases.items():
        fr = {}
        for mode in ("open", "ON"):
            yt = np.array([next(it) for _ in M.FG])
            fr[mode] = M.resonance(yt[:, 0, 0])
        row = {"cell": name, "ls_mm": ls, "pos_mm": pos, "f_OFF_GHz": None if fr["open"] is None else round(fr["open"] / 1e9, 4), "f_ON_GHz": None if fr["ON"] is None else round(fr["ON"] / 1e9, 4)}
        if fr["open"] and fr["ON"]:
            row["shift_ON_minus_OFF_GHz"] = round((fr["ON"] - fr["open"]) / 1e9, 4)
        out["rows"].append(row)
    out["HFSS_reference_A"] = {"C/C (open)": 10.670, "R/R": 11.006, "note": "A pipeline check: model open slit 10.662, R/R 11.014"}
    out["falsifier"] = "ON must be BELOW OFF for family B (Y rises, resonance moves down); if ON is above OFF, or the shift is below 0.05 GHz, family B is not available in this geometry"
    return out


def checks(path, parts):
    res, ok_all = {}, True

    def rec(n, ok, **kw):
        nonlocal ok_all
        ok_all &= bool(ok)
        res[n] = {"ok": bool(ok), **kw}

    names = read_names(path)
    rec("bodies_and_names", set(names) == {"METAL_PEC", "SUBSTRATE", "AIR", "DX_SHEET", "DY_SHEET"}, names=sorted(names))
    solids = names["METAL_PEC"].Solids()
    isl = [s for s in solids if s.BoundingBox().zmin > -1e-6]
    gnd = [s for s in solids if s.BoundingBox().zmax < -1.5]
    vols = [s.Volume() for s in isl]
    rec("four_islands_symmetric", len(isl) == 4 and len(gnd) == 1 and max(vols) - min(vols) < 1e-9, islands=len(isl), island_volume_mm3=vols[0] if vols else None, ground_lumps=len(gnd), note="one named part METAL_PEC; lumps DC-isolated (4 islands + ground)")
    a = M.G / 2
    strip = 2 * (a * np.sqrt(M.R**2 - a**2) + M.R**2 * np.arcsin(a / M.R))  # area of one strip (width G) inside the disc
    disc_area = np.pi * M.R**2 - (2 * strip - 4 * a**2)  # two strips minus their crossing square
    rec("disc_minus_through_slit", abs(sum(vols) - disc_area * M.T) < 1e-6, total_volume=sum(vols), expected=disc_area * M.T, slit_width_mm=M.G)
    gap = min(isl[0].distance(s) for s in isl[1:])
    rec("island_gap", abs(gap - M.G) < 1e-9, min_distance_between_islands_mm=gap)
    sub = names["SUBSTRATE"]
    c1 = sub.intersect(cq.Compound.makeCompound(solids)).Volume()
    rec("no_intersections", c1 < 1e-9, substrate_x_metal=c1)
    dx, dy = names["DX_SHEET"].BoundingBox(), names["DY_SHEET"].BoundingBox()
    okx = abs(dx.xlen - 0.4) < 1e-9 and abs(dx.ylen - 1.0) < 1e-9 and abs(dx.ymin - 1.5) < 1e-9 and abs((dx.xmin + dx.xmax) / 2) < 1e-9 and abs(dx.zmin - 0.0175) < 1e-9
    oky = abs(dy.ylen - 0.4) < 1e-9 and abs(dy.xlen - 1.0) < 1e-9 and abs(dy.xmin - 1.5) < 1e-9 and abs((dy.ymin + dy.ymax) / 2) < 1e-9 and abs(dy.zmin - 0.0175) < 1e-9
    metal = cq.Compound.makeCompound(isl)
    corners = [cq.Vertex.makeVertex(sx * 0.2, y, 0.0175) for sx in (-1, 1) for y in (1.5, 2.5)] + [cq.Vertex.makeVertex(x, sy * 0.2, 0.0175) for sy in (-1, 1) for x in (1.5, 2.5)]
    dist = max(metal.distance(v) for v in corners)
    clear = float(np.hypot(max(dx.xmin - dy.xmax, dy.xmin - dx.xmax, 0), max(dx.ymin - dy.ymax, dy.ymin - dx.ymax, 0)))
    rec("sheets", okx and oky and dist < 1e-6 and clear >= 0.1, DX="0.4 x 1.0 at (0, 2.0)", DY="1.0 x 0.4 at (2.0, 0)", z=0.0175, corner_to_metal_max_mm=float(dist), DX_DY_clearance_mm=round(clear, 4), Rs_ohm_per_sq=2.44 * dy.xlen / dx.xlen)
    ab = names["AIR"].BoundingBox()
    rec("air_box_ground_to_lid_single_face_walls", len(names["AIR"].Faces()) == 6 and abs(ab.zmin + 1.61) < 1e-9 and abs(ab.zmax - M.LID) < 1e-9 and abs(ab.xlen - 6.4) < 1e-9 and abs(ab.ylen - 6.4) < 1e-9, faces=len(names["AIR"].Faces()), z=[ab.zmin, ab.zmax])
    return res, ok_all


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    parts = M.build(ls=LS_B, pos=POS_B, wd=1.0)
    asm = cq.Assembly()
    cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "AIR": (0.7, 0.8, 1.0), "DX_SHEET": (0.9, 0.1, 0.1), "DY_SHEET": (0.1, 0.1, 0.9)}
    for nm, shp in parts.items():
        asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
    path = OUT / "cellB_geometry.step"
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    it = M.intended(ls=LS_B, pos=POS_B, wd=1.0)
    snapped, stats = snap_step_text(raw, it, description="sidequest family-B cell (through slit), mm, 0.001 mm grid")
    path.write_text(snapped, encoding="utf-8", newline="\n")
    bad = inexact_literals(path.read_text(encoding="utf-8"), it)
    vs = verify_snap(raw, snapped)
    res, ok_all = checks(path, parts)
    exact = {"ok": len(bad) == 0 and vs.get("collapsed_edges", 1) == 0 and not vs.get("distinct_merge_violations"), "inexact_literals": len(bad), "n_snapped": stats["n_snapped"], "max_snap_mm": stats["max_snap_mm"], "grid_mm": 0.001}
    res["EXACT_DIMS"] = exact
    ok_all &= exact["ok"]
    (OUT / "geometry_checks.json").write_text(json.dumps({"label": "MODEL_ONLY geometry checks on the re-imported STEP", "all_pass": bool(ok_all), "checks": res}, indent=1, default=float), encoding="utf-8", newline="\n")
    print({k: v["ok"] for k, v in res.items()}, "ALL", ok_all, flush=True)
    if not ok_all:
        sys.exit("checks failed")
    pred = predictions()
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1, default=float), encoding="utf-8", newline="\n")
    print(json.dumps(pred["rows"], indent=0), pred["bare_substrate_GHz_with_lid"], flush=True)
    with open(OUT / "states.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["state", "DX_SHEET", "DY_SHEET", "note"])
        w.writerow(["OFF", "no boundary", "no boundary", "diode OFF = open (HFSS-confirmed model of the OFF diode)"])
        w.writerow(["ON", f"Impedance boundary Rs = {RS} ohm/sq", f"Impedance boundary Rs = {RS} ohm/sq", "Rs = R_on * width / gap = 2.44 * 1.0 / 0.4; Reactance 0; ONE boundary on both sheets"])
    with open(OUT / "results_template.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write("state,run,f_GHz,mode_is_TM0_below_light_line,Ez_fraction,notes\nOFF,x,,,,\nON,x,,,,\nOFF,y,,,,\nON,y,,,,\n")
    row = {r["cell"]: r for r in pred["rows"]}
    rb = row["B_through_slit_pos2.0"]
    setup = f"""# HFSS: family-B unit cell (through-slit RT5880 cell, ON RAISES Y), eigenmode OFF/ON at 90 deg. EXPORT STANDARD. Predictions MODEL_ONLY (`predictions.json`).

**1. Assign material `pec` to `METAL_PEC` (4 island patches + ground = ONE named part, DC-isolated lumps; the STEP imports as vacuum).** `SUBSTRATE`: RT5880 (eps_r 2.2, tan d 0.0009), z -1.575..0.
**2. `AIR`** = one vacuum box z = -1.610 ... +30 mm, 6.4 x 6.4, ONE face per side wall; it overlaps `SUBSTRATE` and `METAL_PEC` on purpose; first in the tree. **Lattice pairs on the four side faces of `AIR` only.** **3. PEC lid:** Perfect E on the top face of `AIR` (z = +30).
**Files:** `cellB_geometry.step` (`METAL_PEC`, `SUBSTRATE`, `AIR`, `DX_SHEET`, `DY_SHEET`), `states.csv`, `predictions.json`, `geometry_checks.json` (all PASS on the re-imported STEP), `results_template.csv`.
**Geometry** (mm): disc r 3.0 (35 um copper, z 0..0.035) cut by a cross slit 0.4 wide THROUGH the disc (four islands, no rim); ground 35 um (z -1.610..-1.575). Sheets at z = 0.0175: `DX_SHEET` 0.4 x 1.0 across the vertical slit at y = +1.5..2.5; `DY_SHEET` 1.0 x 0.4 across the horizontal slit at x = +1.5..2.5; edges on the island faces.
**4. States:** OFF = **no boundary** on the sheets. ON = **ONE Impedance boundary on `DX_SHEET` and `DY_SHEET`: Rs = {RS} ohm/sq (= 2.44 ohm x 1.0 / 0.4 from the real sheet size), reactance 0** (not an RLC).
**5. Eigenmode:** Lattice Pair 1 (x faces, Phase Delay `$ph_x`), Pair 2 (y faces, `$ph_y`); x run `$ph_x = 90deg`, `$ph_y = 0deg`; y run swapped (the cell is x/y symmetric: expect equal values). Min frequency 5 GHz, 8 modes, delta f 0.1 %; mesh: substrate 1 mm, islands 0.2 mm, sheets 0.1 mm. Clean Up Solutions before each run. The bound mode is the TM0 below the light line 11.711 GHz (track it by the Ez fraction in each state).
**Record** (`results_template.csv`): f of the TM0 mode for OFF and ON, x and y.
**Prediction (MODEL_ONLY, +-0.15 GHz):** family B OFF {rb['f_OFF_GHz']} GHz, ON {rb['f_ON_GHz']} GHz (shift {rb.get('shift_ON_minus_OFF_GHz')} GHz); pipeline check on the HFSS-validated family-A cell: model {row['A_rim_pos0.8 (HFSS-validated cell)']['f_OFF_GHz']} / {row['A_rim_pos0.8 (HFSS-validated cell)']['f_ON_GHz']} GHz against HFSS 10.670 / 11.006.
**Falsifier:** ON must lie BELOW OFF (Y rises); ON above OFF, or |shift| < 0.05 GHz, means family B is not available in this geometry. Report the shift; a board for design (ii) is built only after this cell passes.
"""
    (OUT / "HFSS_SETUP.md").write_text(setup, encoding="utf-8", newline="\n")
    SYNC.mkdir(parents=True, exist_ok=True)
    for f in OUT.iterdir():
        shutil.copy2(f, SYNC / f.name)
    print("synced", sorted(p.name for p in SYNC.iterdir()))


if __name__ == "__main__":
    main()
