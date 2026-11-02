"""Phase 30 item 3 (MODEL_ONLY): HFSS unit-cell EIGENMODE package for the approved cell family -- RO3010 0.635 mm, 3.2 mm
pitch, full-fill super-ellipse (n = 4) with a centre slot along x -- with the model's predictions written BEFORE any HFSS
run.  Re-uses the proven master/slave + PEC-lid construction of `phase26_eigenmode_package.py` (its template is patched,
not copied); the patch is a native polyline (72 vertices) with a rectangular slot subtracted.

    python scripts/research/phase30_eigenmode_package.py <refine_dir> [--mesh 48] [--check-mesh 64]

Cells E1..E4 (picked from `results/phase30/refine/refine_grid.json`): Bloch phase per cell chosen at the model's mean beta
at 10.54 GHz; the model predicts the eigenfrequency of the bound mode along x (slot axis) and along y.  S1: 7-cell
rotating-axis SUPERCELL of E1 with a radiation boundary on top (leaky eigenmode): the predicted complex frequency
(Q) gives alpha via alpha = Im(f) * 2 pi / v_g.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase26_eigenmode_package as p26  # noqa: E402

OUT = ROOT / "hfss_exports" / "phase30_cell_eigenmode"
P_MM, H_MM = 3.2, 0.635
EPS_R, TAND = 10.2, 0.0022
K0 = 2 * np.pi * 10.54e9 / 299792458.0
POLY_CODE = '''
def squircle_pts(a, b, n_exp, rot_deg, nv=72):
    pts = []
    for k in range(nv):
        t = 2.0 * math.pi * k / nv
        c, s = math.cos(t), math.sin(t)
        u = 0.5 * a * math.copysign(abs(c) ** (2.0 / n_exp), c)
        v = 0.5 * b * math.copysign(abs(s) ** (2.0 / n_exp), s)
        pts.append((u, v))
    r = math.radians(rot_deg); cr, sr = math.cos(r), math.sin(r)
    return [(u * cr - v * sr, u * sr + v * cr) for (u, v) in pts]

def poly(name, pts, z, cx=0.0, cy=0.0):
    oEditor.CreatePolyline(["NAME:PolylineParameters", "IsPolylineCovered:=", True, "IsPolylineClosed:=", True,
        ["NAME:PolylinePoints"] + [["NAME:PLPoint", "X:=", mm(cx + x), "Y:=", mm(cy + y), "Z:=", mm(z)] for (x, y) in pts + [pts[0]]],
        ["NAME:PolylineSegments"] + [["NAME:PLSegment", "SegmentType:=", "Line", "StartIndex:=", k, "NoOfPoints:=", 2] for k in range(len(pts))],
        ["NAME:PolylineXSection", "XSectionType:=", "None", "XSectionOrient:=", "Auto", "XSectionWidth:=", "0mm",
         "XSectionTopWidth:=", "0mm", "XSectionHeight:=", "0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=", "Corner"]],
        attrs(name, "vacuum", True))

def make_patch(name, a, b, n_exp, sl, sw, rot_deg, cx=0.0, cy=0.0):
    poly(name, squircle_pts(a, b, n_exp, rot_deg), H, cx, cy)
    if sl > 0.0 and sw > 0.0:
        r = math.radians(rot_deg); cr, sr = math.cos(r), math.sin(r)
        rc = [(-sl / 2, -sw / 2), (sl / 2, -sw / 2), (sl / 2, sw / 2), (-sl / 2, sw / 2)]
        poly(name + "_slot", [(u * cr - v * sr, u * sr + v * cr) for (u, v) in rc], H, cx, cy)
        oEditor.Subtract(["NAME:Selections", "Blank Parts:=", name, "Tool Parts:=", name + "_slot"],
                         ["NAME:SubtractParameters", "KeepOriginals:=", False])
'''


def patch_block(cells):
    """cells: list of (cx, cy, rot_deg); a, b, n, slot from the script constants."""
    lines = ["PATCHES = []"]
    for k, (cx, cy, rot) in enumerate(cells):
        lines.append('make_patch("Patch%d", MAJOR, MINOR, EXPO, SLOT_L, SLOT_W, %.4f, %.4f, %.4f); PATCHES.append("Patch%d")' % (k, rot, cx, cy, k))
    return "\n".join(lines)


def build_script(name, spec, n_cells=1, radiation=False, phase_x=0.0, eps=EPS_R):
    t = p26.TEMPLATE
    start = t.index('oEditor.CreateCircle(["NAME:CircleParameters"')
    end = t.index('oEditor.CreateRectangle(["NAME:RectangleParameters", "IsCovered:=", True, "XStart:=", mm(-P / 2), "YStart:=", mm(-P / 2),\n    "ZStart:=", "0mm"')
    if n_cells == 1:
        cells = [(0.0, 0.0, 0.0)]
    else:
        cells = [(-(n_cells * P_MM) / 2 + (i + 0.5) * P_MM, 0.0, 180.0 * i / n_cells) for i in range(n_cells)]
    t = t[:start] + POLY_CODE + "\n" + patch_block(cells) + "\n\n" + t[end:]
    if radiation:  # open top: no PEC lid sheet
        i0 = t.index('oEditor.CreateRectangle(["NAME:RectangleParameters", "IsCovered:=", True, "XStart:=", mm(-P / 2), "YStart:=", mm(-P / 2),' + chr(10) + '    "ZStart:=", mm(H + LID)')
        i1 = t.index(chr(10) + chr(10), i0)
        t = t[:i0] + t[i1:]
    t = t.replace("P = 6.4; H = 1.575; LID = @@LID@@; D = @@D@@", "P = %g; H = %g; LID = @@LID@@; D = 0.0\nPX_CELLS = %d\nMAJOR = %r; MINOR = %r; EXPO = 4.0" % (P_MM, H_MM, n_cells, spec["major_mm"], spec["major_mm"] * spec["ratio"]))
    t = t.replace("SLOT_L = @@SLOT_L@@; SLOT_W = @@SLOT_W@@; SLOT_ANGLE = @@ANGLE@@", "SLOT_L = @@SLOT_L@@; SLOT_W = @@SLOT_W@@")
    if n_cells > 1:
        # supercell: the box is N cells long in x
        t = t.replace("mm(-P / 2), \"YPosition:=\", mm(-P / 2)", "mm(-P * PX_CELLS / 2), \"YPosition:=\", mm(-P / 2)")
        t = t.replace('"XSize:=", mm(P), "YSize:=", mm(P)', '"XSize:=", mm(P * PX_CELLS), "YSize:=", mm(P)')
        t = t.replace('"XStart:=", mm(-P / 2), "YStart:=", mm(-P / 2)', '"XStart:=", mm(-P * PX_CELLS / 2), "YStart:=", mm(-P / 2)')
        t = t.replace('"Width:=", mm(P), "Height:=", mm(P)', '"Width:=", mm(P * PX_CELLS), "Height:=", mm(P)')
        t = t.replace("(\"X\", (-P / 2, 0.0), (P / 2, 0.0), PHASE_X, ((-P / 2, -P / 2), (-P / 2, P / 2)), ((P / 2, -P / 2), (P / 2, P / 2)))",
                      "(\"X\", (-P * PX_CELLS / 2, 0.0), (P * PX_CELLS / 2, 0.0), PHASE_X, ((-P * PX_CELLS / 2, -P / 2), (-P * PX_CELLS / 2, P / 2)), ((P * PX_CELLS / 2, -P / 2), (P * PX_CELLS / 2, P / 2)))")
        t = t.replace("(\"Y\", (0.0, -P / 2), (0.0, P / 2), PHASE_Y, ((-P / 2, -P / 2), (P / 2, -P / 2)), ((-P / 2, P / 2), (P / 2, P / 2)))",
                      "(\"Y\", (0.0, -P / 2), (0.0, P / 2), PHASE_Y, ((-P * PX_CELLS / 2, -P / 2), (P * PX_CELLS / 2, -P / 2)), ((-P * PX_CELLS / 2, P / 2), (P * PX_CELLS / 2, P / 2)))")
        t = t.replace("oBnd.AssignPerfectE([\"NAME:PerfE_\" + nm, \"Objects:=\", [nm], \"InfGroundPlane:=\", False])\n",
                      "oBnd.AssignPerfectE([\"NAME:PerfE_\" + nm, \"Objects:=\", [nm], \"InfGroundPlane:=\", False])\n")
    t = t.replace('for nm in ("Patch", "Ground", "Lid"):', 'for nm in PATCHES + ["Ground"] + ([] if %s else ["Lid"]):' % ("True" if radiation else "False"))
    if radiation:
        t = t.replace('"XSize:=", mm(P), "YSize:=", mm(P), "ZSize:=", mm(LID)]', '"XSize:=", mm(P), "YSize:=", mm(P), "ZSize:=", mm(LID)]')
        t = t.replace("oAna = oDesign.GetModule(\"AnalysisSetup\")", 'oBnd.AssignRadiation(["NAME:RadTop", "Faces:=", [int(oEditor.GetFaceByPosition(["NAME:FaceParameters", "BodyName:=", "Air", "XPosition:=", "0mm", "YPosition:=", "0mm", "ZPosition:=", mm(H + LID)]))], "IsFssReference:=", False, "IsForPML:=", False])\noAna = oDesign.GetModule("AnalysisSetup")', 1)
    header_old = [ln for ln in t.split("\n") if ln.startswith("# HFSS (AEDT) IronPython script")][0]
    t = t.replace(header_old, "# HFSS (AEDT) IronPython script -- unit-cell EIGENMODE of the Phase 30 cell \"@@NAME@@\" (generated by lwmom). NOT run in AEDT by the author. MODEL_ONLY predictions: see README / predictions.json.")
    t = t.replace("#   cell 6.4 x 6.4 mm, RT/duroid 5880", "#   cell 3.2 x 3.2 mm (x N cells for a supercell), Rogers RO3010")
    vals = {
        "NAME": name, "EPS": repr(eps), "TAND": repr(TAND), "LID": "60.0",
        "SLOT_L": repr(spec["slot_l_mm"]), "SLOT_W": repr(spec["slot_w_mm"]),
        "SLOTDESC": "super-ellipse n=4 %.3f x %.3f mm, slot %.2f x %.2f mm along x" % (spec["major_mm"], spec["major_mm"] * spec["ratio"], spec["slot_l_mm"], spec["slot_w_mm"]),
        "ANGLE": "0.0", "D": "0.0", "PX": "%.3f" % phase_x, "PD": "0",
    }
    for k, v in vals.items():
        t = t.replace("@@%s@@" % k, v)
    t = t.replace('PHASE_X = "%.3fdeg"' % phase_x, 'PHASE_X = "%.3fdeg"' % phase_x)
    return t


def main(refine_dir, mesh_n=48, check_mesh=64):
    OUT.mkdir(parents=True, exist_ok=True)
    stage = json.loads((Path(refine_dir) / "refine_stage.json").read_text())["picks"]
    preds = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "PREDICTED_BEFORE_HFSS"], "cells": {}}
    from src.basis_functions import MeshSettings
    from src.dispersion import tensor_sheet as ts
    import phase30_refine as R

    names = {0: "E1_strong_M0p64", 1: "E2_longthin_M0p62", 2: "E3_mid_M0p27", 3: "E4_weak_M0p12"}
    for pi, nm in names.items():
        pk = stage[pi]
        spec = pk["spec"]
        f_row = [r for r in pk["freq"] if abs(r["f_hz"] - 10.54e9) < 1e3][0]
        beta = f_row["beta_over_k0"] * K0
        phase = float(np.degrees(beta * P_MM * 1e-3))
        cs = {"major_mm": spec["major"] * 1e3, "ratio": spec["ratio"], "slot_l_mm": spec["slot_l"] * 1e3, "slot_w_mm": spec["slot_w"] * 1e3}
        (OUT / ("build_%s.py" % nm)).write_text(build_script(nm, cs, phase_x=phase), encoding="utf-8", newline="\n")
        entry = {"cell": cs, "bloch_phase_deg_per_cell": phase, "beta_over_k0_at_10p54": f_row["beta_over_k0"], "meshes": {}}
        cell = R.make_cell(spec["major"], spec["ratio"], spec["slot_l"], spec["slot_w"], 0.0)
        for mn in (mesh_n, check_mesh):
            mesh = MeshSettings(mn, mn, mn)
            kx = np.radians(phase) / (P_MM * 1e-3)
            rx = ts.bloch_frequencies(cell, mesh, kx, 0.0, f_lo=8.5e9, f_hi=12.5e9)
            ry = ts.bloch_frequencies(cell, mesh, 0.0, kx, f_lo=8.5e9, f_hi=12.5e9)

            def tm(rs):
                return sorted(rs, key=lambda r: (-round(r[1], 2), r[0]))[0][0] / 1e9 if rs else None

            entry["meshes"][str(mn)] = {"f_along_x_ghz": tm(rx), "f_along_y_ghz": tm(ry),
                                         "roots_x": [[round(f / 1e9, 4), round(t, 2)] for f, t in rx], "roots_y": [[round(f / 1e9, 4), round(t, 2)] for f, t in ry]}
            print(nm, mn, entry["meshes"][str(mn)]["f_along_x_ghz"], entry["meshes"][str(mn)]["f_along_y_ghz"], flush=True)
        preds["cells"][nm] = entry
        if pi == 0:
            n_cell = f_row.get("N_cells_period", 7)
            lam = n_cell * P_MM
            beta_super = f_row.get("beta_floquet_over_k0", f_row["beta_over_k0"]) * K0
            ph_super = float((np.degrees(beta_super * lam * 1e-3)) % 360.0)
            fr = sorted(pk["freq"], key=lambda r: r["f_hz"])
            dbeta = (fr[-1]["beta_over_k0"] * (fr[-1]["f_hz"] / 10.54e9) - fr[0]["beta_over_k0"] * (fr[0]["f_hz"] / 10.3e9)) * K0
            # group velocity v_g = d omega / d beta from the three-frequency model beta(f)
            b_lo = fr[0]["beta_over_k0"] * 2 * np.pi * fr[0]["f_hz"] / 299792458.0
            b_hi = fr[-1]["beta_over_k0"] * 2 * np.pi * fr[-1]["f_hz"] / 299792458.0
            vg = 2 * np.pi * (fr[-1]["f_hz"] - fr[0]["f_hz"]) / (b_hi - b_lo)
            a_f = f_row.get("alpha_floquet_over_k0")
            a1 = f_row["alpha_over_k0"]
            sup = {"N": n_cell, "period_mm": lam, "bloch_phase_deg_per_supercell": ph_super, "v_g_over_c": float(vg / 299792458.0),
                   "alpha_first_order_over_k0": a1, "alpha_floquet_over_k0": a_f}
            for key, a in (("first_order", a1), ("floquet", a_f)):
                if a is not None:
                    alpha_np = a * K0
                    imf = alpha_np * vg / (2 * np.pi)
                    sup["im_f_ghz_" + key] = float(imf / 1e9)
                    sup["Q_" + key] = float(10.54e9 / (2 * imf))
            preds["supercell_S1"] = sup
            sp = dict(cs)
            (OUT / "build_S1_supercell_leaky.py").write_text(build_script("S1_supercell_leaky", sp, n_cells=n_cell, radiation=True, phase_x=ph_super), encoding="utf-8", newline="\n")
    (OUT / "predictions.json").write_text(json.dumps(preds, indent=1))
    print("done")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("refine_dir")
    ap.add_argument("--mesh", type=int, default=48)
    ap.add_argument("--check-mesh", type=int, default=64)
    a = ap.parse_args()
    main(a.refine_dir, a.mesh, a.check_mesh)
