"""Phase 30 item 3b (MODEL_ONLY): single-feed TEST BOARD plan to separate cell leakage from monopole direct radiation and edge
diffraction.  RO3010 0.635 mm, 3.2 mm lattice, full-fill super-ellipse (n=4) cells with a centre slot, one coax-probe monopole
at the centre, 90 mm circular board, broadside RHCP beam, M-tapered (the slot length sets M, the slot axis the phase).

    python scripts/research/phase30_testboard.py <refine_grid.json> [--eta 0.8] [--taper 0.5]

Outputs `hfss_exports/phase30_testboard/`: `build_board_active.py` (tapered tensor-CP board), `build_board_control.py` (same
outline and feed, every cell the unslotted squircle = isotropic, no leakage), `layout.csv`, `predictions.json`, `README.md`.

Design chain = the Phase 26 ray march (`tensor_leaky_cp.march_rays` / `cell_design`) fed with a `SupercellTable` built from
the Phase 30 refined-grid cells (first-order alpha, beta of each cell; K/k0 dependence NOT modelled: the table is replicated
over K columns).  Far field: `_vector_far_field` with the RO3010 layered Green's function.  Launch fraction of the monopole
is NOT modelled (NA); gain is given for launch fraction 1 (an upper bound).
"""

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.frequency import free_space_wavenumber  # noqa: E402
from src.green_function import LayeredGreenFunction  # noqa: E402
from src.hologram.tensor_leaky_cp import SupercellTable, cell_design, march_rays  # noqa: E402
from src.materials import Substrate  # noqa: E402
from src.pattern.vector_aperture_pattern import _vector_far_field, circular_components, copol_component  # noqa: E402

OUT = ROOT / "hfss_exports" / "phase30_testboard"
F = float(os.environ.get("P30_F0", 10.54e9))  # default unchanged
K0 = free_space_wavenumber(F)
PITCH = 3.2e-3
MAJOR_MM = 3.05
EXPO = 4.0
H_MM = 0.635
R_BOARD = 0.045  # 90 mm circle
RHO_CLEAR = 4.0e-3  # keep-out for the coax opening (r_out 2.05 mm) + one cell
SUB = Substrate("ro3010", float(os.environ.get("P30_EPS", 10.2)), 0.0022, H_MM * 1e-3)
M_CAP = 0.66


def build_table(grid_json):
    rows = [r for r in json.loads(Path(grid_json).read_text())["rows"] if r["ok"] and abs(r["ratio"] - 1.0) < 1e-9]
    cand = [r for r in rows if r["M_t"] <= M_CAP and r["aniso_ratio"] < 10 and r["x_par_ohm"] > 0 and r["x_perp_ohm"] > 0]
    uniq = {}
    for r in cand:  # W = 0.1 / 0.2 mm give the same mesh-48 result: keep the thinner
        key = (round(r["slot_l"] * 1e3, 3), round(r["slot_w"] * 1e3, 3))
        uniq[key] = r
    keys = sorted(uniq, key=lambda k: uniq[k]["alpha_np_m"])
    k_grid = np.arange(1.0, 1.61, 0.05)
    alpha = np.array([[uniq[k]["alpha_np_m"]] * k_grid.size for k in keys])
    beta = np.array([[uniq[k]["beta_over_k0"]] * k_grid.size for k in keys])
    table = SupercellTable(tuple(keys), k_grid, alpha, beta, MAJOR_MM, PITCH * 1e3)
    return table, uniq, keys


def circle_poly(r, n=96):
    return np.array([[r * np.cos(2 * np.pi * i / n), r * np.sin(2 * np.pi * i / n)] for i in range(n)])


def lattice(r_board):
    n = int(np.ceil(r_board / PITCH)) + 2
    k = np.arange(-n, n + 1)
    xg, yg = np.meshgrid(PITCH * (k + 0.5), PITCH * (k + 0.5))  # feed at a lattice-cell corner: 4-fold symmetric lattice
    x, y = xg.ravel(), yg.ravel()
    rho = np.hypot(x, y)
    ok = (rho <= r_board - 0.5 * PITCH) & (rho >= RHO_CLEAR)
    return x[ok], y[ok]


def far_field(x, y, d, green, theta, phi):
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    th, ph = tg.ravel(), pg.ravel()
    et = np.zeros(th.size, complex)
    ep = np.zeros(th.size, complex)
    for lo in range(0, th.size, 3000):
        sl = slice(lo, lo + 3000)
        et[sl], ep[sl] = _vector_far_field(th[sl], ph[sl], x, y, d["jx"], d["jy"], F, green, PITCH, PITCH)
    return et.reshape(tg.shape), ep.reshape(tg.shape), tg, pg


def pattern_numbers(et, ep, tg, pg, theta, phi):
    u = np.abs(et) ** 2 + np.abs(ep) ** 2
    dth, dph = np.radians(theta[1] - theta[0]), np.radians(phi[1] - phi[0])
    w = np.sin(np.radians(tg)) * dth * dph
    p_tot = float(np.sum(u * w))
    rh = np.abs(copol_component(et, ep, "rhcp")) ** 2
    lh = np.abs(copol_component(et, ep, "lhcp")) ** 2
    i = np.unravel_index(np.argmax(rh), rh.shape)
    d_rh = float(4 * np.pi * rh[i] / p_tot)
    cop, cross = np.sqrt(rh[i]), np.sqrt(lh[i])
    ar = 20 * np.log10((cop + cross) / max(abs(cop - cross), 1e-12))
    dpk = np.degrees(np.arccos(np.clip(np.cos(np.radians(tg)) * np.cos(np.radians(tg[i])) + np.sin(np.radians(tg)) * np.sin(np.radians(tg[i])) * np.cos(np.radians(pg) - np.radians(pg[i])), -1, 1)))
    side = rh[dpk > 2.0 * 12.0]
    near = float(np.sum(u * w * (tg > 60.0)) / p_tot)
    cone = float(np.sum(u * w * (tg <= 12.0)) / p_tot)
    return {"directivity_rhcp_dbi": 10 * np.log10(d_rh), "ar_at_peak_db": float(ar), "peak_theta_deg": float(tg[i]), "peak_phi_deg": float(pg[i]),
            "rhcp_sll_db_outside_24deg": float(10 * np.log10(side.max() / rh[i])) if side.size else None,
            "radiated_fraction_theta_gt_60": near, "radiated_fraction_theta_le_12": cone,
            "lhcp_over_rhcp_at_peak_db": float(20 * np.log10(max(cross, 1e-12) / max(cop, 1e-12)))}


def aedt_board(cells, name, note):
    rows = ",\n    ".join("(%.3f, %.3f, %.3f, %.3f, %.4f)" % c for c in cells)
    return '''# -*- coding: utf-8 -*-
# HFSS (AEDT) IronPython script -- Phase 30 single-feed TEST BOARD "%(name)s" (generated by lwmom). MODEL_ONLY / NOT_HFSS_VALIDATED.
# Run inside AEDT: Tools > Run Script...  Native geometry in mm. %(note)s
#   Substrate Rogers RO3010 (eps_r 10.2, tan d 0.0022) 0.635 mm, 100 x 100 mm, z = 0 .. 0.635 mm; PEC ground with a coax opening
#   (r_out 2.05 mm); 90 mm circular patch field of 3.2 mm cells (full-fill super-ellipse n=4, major 3.05 mm, centre slot along the
#   cell axis, angle per cell); monopole at the board centre (pin r 0.635 mm, 6.0 mm above the substrate), lumped 50 ohm port.
import math
import ScriptEnv
ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oProject = oDesktop.NewProject()
oProject.InsertDesign("HFSS", "%(name)s", "DrivenModal", "")
oDesign = oProject.SetActiveDesign("%(name)s")
oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.SetModelUnits(["NAME:Units Parameter", "Units:=", "mm", "Rescale:=", False])
ADD_SWEEP = True
H = "%(h)gmm"
WX = "100mm"; WY = "100mm"
MAJOR = %(major)g; EXPO = %(expo)g; NV = 64
CELLS = [  # (x, y, slot length, slot width, slot angle deg) mm; slot centred, long axis at the angle from +x
    %(rows)s
]
oProject.GetDefinitionManager().AddMaterial(["NAME:HMSA_substrate", "CoordinateSystemType:=", "Cartesian",
    "BulkOrSurfaceType:=", 1, ["NAME:PhysicsTypes", "set:=", ["Electromagnetic"]],
    "permittivity:=", "10.2", "dielectric_loss_tangent:=", "0.0022"])
def var(name, value):
    oDesign.ChangeProperty(["NAME:AllTabs", ["NAME:LocalVariableTab", ["NAME:PropServers", "LocalVariables"],
        ["NAME:NewProps", ["NAME:" + name, "PropType:=", "VariableProp", "UserDef:=", True, "Value:=", value]]]])
for n, v in (("mono_h", "%(mono_h)gmm"), ("r_pin", "0.635mm"), ("r_out", "2.05mm"), ("pad_xy", "15mm"), ("pad_top", "30mm"), ("pad_bot", "15mm")):
    var(n, v)
def attrs(name, material, solve_inside, color="(255 128 0)", transp=0):
    return ["NAME:Attributes", "Name:=", name, "Flags:=", "", "Color:=", color, "Transparency:=", transp,
            "PartCoordinateSystem:=", "Global", "UDMId:=", "", "MaterialValue:=", '"' + material + '"',
            "SurfaceMaterialValue:=", '""', "SolveInside:=", solve_inside, "IsMaterialEditable:=", True,
            "UseMaterialAppearance:=", False, "IsLightweight:=", False]
def mm(v):
    return repr(float(v)) + "mm"
def rect(name, x, y, z, w, h, color="(255 128 0)"):
    oEditor.CreateRectangle(["NAME:RectangleParameters", "IsCovered:=", True, "XStart:=", x, "YStart:=", y, "ZStart:=", z,
        "Width:=", w, "Height:=", h, "WhichAxis:=", "Z"], attrs(name, "vacuum", True, color))
def circle(name, x, y, z, r, color="(255 128 0)"):
    oEditor.CreateCircle(["NAME:CircleParameters", "IsCovered:=", True, "XCenter:=", x, "YCenter:=", y, "ZCenter:=", z,
        "Radius:=", r, "WhichAxis:=", "Z", "NumSegments:=", "0"], attrs(name, "vacuum", True, color))
def subtract(blank, tools, keep):
    oEditor.Subtract(["NAME:Selections", "Blank Parts:=", blank, "Tool Parts:=", ",".join(tools)],
                     ["NAME:SubtractParameters", "KeepOriginals:=", keep])
def poly_sheet(name, pts, z):
    oEditor.CreatePolyline(["NAME:PolylineParameters", "IsPolylineCovered:=", True, "IsPolylineClosed:=", True,
        ["NAME:PolylinePoints"] + [["NAME:PLPoint", "X:=", mm(x), "Y:=", mm(y), "Z:=", z] for (x, y) in pts + [pts[0]]],
        ["NAME:PolylineSegments"] + [["NAME:PLSegment", "SegmentType:=", "Line", "StartIndex:=", k, "NoOfPoints:=", 2] for k in range(len(pts))],
        ["NAME:PolylineXSection", "XSectionType:=", "None", "XSectionOrient:=", "Auto", "XSectionWidth:=", "0mm",
         "XSectionTopWidth:=", "0mm", "XSectionHeight:=", "0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=", "Corner"]],
        attrs(name, "vacuum", True))
def squircle(cx, cy):
    pts = []
    for k in range(NV):
        t = 2.0 * math.pi * k / NV
        c, s = math.cos(t), math.sin(t)
        pts.append((round(cx + 0.5 * MAJOR * math.copysign(abs(c) ** (2.0 / EXPO), c), 4), round(cy + 0.5 * MAJOR * math.copysign(abs(s) ** (2.0 / EXPO), s), 4)))
    return pts
x_lo = "(-" + WX + "/2)"; y_lo = "(-" + WY + "/2)"
oEditor.CreateBox(["NAME:BoxParameters", "XPosition:=", x_lo, "YPosition:=", y_lo, "ZPosition:=", "0mm",
    "XSize:=", WX, "YSize:=", WY, "ZSize:=", H], attrs("Substrate", "HMSA_substrate", True, "(0 128 0)", 0.6))
rect("Ground", x_lo, y_lo, "0mm", WX, WY, "(128 128 128)")
oBnd = oDesign.GetModule("BoundarySetup")
oEditor.CreateCylinder(["NAME:CylinderParameters", "XCenter:=", "0mm", "YCenter:=", "0mm", "ZCenter:=", "0mm",
    "Radius:=", "r_pin", "Height:=", "mono_h", "WhichAxis:=", "Z", "NumSides:=", "0"], attrs("Monopole1", "pec", False, "(255 255 0)"))
circle("Hole1", "0mm", "0mm", "0mm", "r_out")
circle("PortSheet1", "0mm", "0mm", "0mm", "r_out")
circle("PinCut1", "0mm", "0mm", "0mm", "r_pin")
subtract("PortSheet1", ["PinCut1"], False)
subtract("Ground", ["Hole1"], False)
subtract("Substrate", ["Monopole1"], True)
oBnd.AssignLumpedPort(["NAME:Port1", "Objects:=", ["PortSheet1"], "DoDeembed:=", False, "RenormalizeAllTerminals:=", True,
    ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", True,
        ["NAME:IntLine", "Start:=", ["r_pin", "0mm", "0mm"], "End:=", ["r_out", "0mm", "0mm"]],
        "AlignmentGroup:=", 0, "CharImp:=", "Zpi", "RenormImp:=", "50ohm"]],
    "ShowReporterFilter:=", False, "ReporterFilter:=", [True], "Impedance:=", "50ohm"])
names = []
for k, (cx, cy, L, W, ang) in enumerate(CELLS):
    nm = "C%%d" %% k
    poly_sheet(nm, squircle(cx, cy), H)
    if L > 0.0 and W > 0.0:
        a = math.radians(ang); c = math.cos(a); s = math.sin(a)
        pts = [(round(cx + u * c - v * s, 4), round(cy + u * s + v * c, 4)) for (u, v) in ((-L / 2, -W / 2), (L / 2, -W / 2), (L / 2, W / 2), (-L / 2, W / 2))]
        poly_sheet("S%%d" %% k, pts, H)
        subtract(nm, ["S%%d" %% k], False)
    names.append(nm)
oBnd.AssignPerfectE(["NAME:PerfE_Ground", "Objects:=", ["Ground"], "InfGroundPlane:=", False])
oBnd.AssignPerfectE(["NAME:PerfE_Patches", "Objects:=", names, "InfGroundPlane:=", False])
oEditor.CreateBox(["NAME:BoxParameters", "XPosition:=", x_lo + "-pad_xy", "YPosition:=", y_lo + "-pad_xy", "ZPosition:=", "-pad_bot",
    "XSize:=", WX + "+2*pad_xy", "YSize:=", WY + "+2*pad_xy", "ZSize:=", "pad_bot+" + H + "+pad_top"],
    attrs("Air", "air", True, "(143 175 143)", 0.9))
oBnd.AssignRadiation(["NAME:Rad1", "Objects:=", ["Air"], "IsFssReference:=", False, "IsForPML:=", False])
oAna = oDesign.GetModule("AnalysisSetup")
oAna.InsertSetup("HfssDriven", ["NAME:Setup1", "Frequency:=", "10.54GHz", "MaxDeltaS:=", 0.02, "MaximumPasses:=", 15,
    "MinimumPasses:=", 2, "MinimumConvergedPasses:=", 2, "PercentRefinement:=", 30, "IsEnabled:=", True, "BasisOrder:=", 1,
    "DoLambdaRefine:=", True, "DoMaterialLambda:=", True, "SetLambdaTarget:=", False, "Target:=", 0.3333, "UseMaxTetIncrease:=", False])
if ADD_SWEEP:
    oAna.InsertFrequencySweep("Setup1", ["NAME:Sweep", "IsEnabled:=", True, "RangeType:=", "LinearStep", "RangeStart:=", "10.0GHz",
        "RangeEnd:=", "11.1GHz", "RangeStep:=", "0.02GHz", "Type:=", "Interpolating", "SaveFields:=", False, "SaveRadFields:=", False,
        "InterpTolerance:=", 0.5, "InterpMaxSolns:=", 250, "InterpMinSolns:=", 0, "InterpMinSubranges:=", 1])
oRad = oDesign.GetModule("RadField")
oRad.InsertFarFieldSphereSetup(["NAME:3D", "UseCustomRadiationSurface:=", False, "ThetaStart:=", "0deg", "ThetaStop:=", "180deg",
    "ThetaStep:=", "1deg", "PhiStart:=", "-180deg", "PhiStop:=", "180deg", "PhiStep:=", "1deg", "UseLocalCS:=", False])
oEditor.FitAll()
oProject.SaveAs(oDesktop.GetProjectDirectory() + "/%(name)s.aedt", True)
oDesktop.AddMessage(oProject.GetName(), "%(name)s", 0, "Built: %%d cells, 1 monopole, 1 lumped port. Validate, then Analyze All." %% len(CELLS))
''' % {"name": name, "note": note, "h": H_MM, "major": MAJOR_MM, "expo": EXPO, "rows": rows, "mono_h": H_MM + 6.0}


def main(grid_json, eta=0.8, taper=0.5):
    OUT.mkdir(parents=True, exist_ok=True)
    table, uniq, keys = build_table(grid_json)
    x, y = lattice(R_BOARD)
    poly = circle_poly(R_BOARD)
    b0 = float(np.min(table.beta[:, 0]))
    m = march_rays(table, poly, (0.0, 0.0), K0, 0.0, 0.0, eta, taper_c=taper, n_psi=360, rho_min_m=RHO_CLEAR, beta0_over_k0=b0)
    d = cell_design(m, table, x, y, (0.0, 0.0), K0, 0.0, 0.0, hand=+1)
    theta = np.arange(0.0, 90.5, 1.0)
    phi = np.arange(0.0, 360.0, 5.0)
    green = LayeredGreenFunction(SUB, F, backing="pec")
    et, ep, tg, pg = far_field(x, y, d, green, theta, phi)
    pn = pattern_numbers(et, ep, tg, pg, theta, phi)
    eta_conv = float(m["eta_conv"])
    mt = [uniq[k]["M_t"] for k in keys]
    used = np.asarray(d["slot_index"])
    m_used = np.array([mt[i] for i in used])
    cells = [(float(x[i] * 1e3), float(y[i] * 1e3), float(d["slot_L_mm"][i]), float(d["slot_W_mm"][i]), float(np.degrees(d["slot_angle_rad"][i]))) for i in range(x.size)]
    ctrl = [(c[0], c[1], 0.0, 0.0, 0.0) for c in cells]
    (OUT / "build_board_active.py").write_text(aedt_board(cells, "phase30_testboard_active", "ACTIVE: tapered tensor-CP board."), encoding="utf-8", newline="\n")
    (OUT / "build_board_control.py").write_text(aedt_board(ctrl, "phase30_testboard_control", "CONTROL: same outline/feed, every cell the unslotted squircle (isotropic, no leakage)."), encoding="utf-8", newline="\n")
    with open(OUT / "layout.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["x_mm", "y_mm", "slot_L_mm", "slot_W_mm", "slot_angle_deg", "alpha_np_per_m", "beta_over_k0", "rho_mm", "M_t_cell"])
        for i, c in enumerate(cells):
            w.writerow(["%.3f" % c[0], "%.3f" % c[1], "%.3f" % c[2], "%.3f" % c[3], "%.3f" % c[4], "%.3f" % d["alpha"][i], "%.4f" % d["beta_over_k0"][i], "%.2f" % (d["rho"][i] * 1e3), "%.3f" % m_used[i]])
    rho = d["rho"] * 1e3
    prof = {}
    for lo, hi in ((4, 12), (12, 20), (20, 28), (28, 36), (36, 45)):
        sel = (rho >= lo) & (rho < hi)
        if sel.any():
            prof["%d-%d mm" % (lo, hi)] = {"alpha_np_m_mean": float(np.mean(d["alpha"][sel])), "M_t_mean": float(np.mean(m_used[sel])), "beta_over_k0_mean": float(np.mean(d["beta_over_k0"][sel]))}
    pred = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "LAUNCH_FRACTION_NOT_MODELLED"],
        "board": {"shape": "circle", "diameter_mm": 2 * R_BOARD * 1e3, "cells": int(x.size), "substrate": "RO3010 0.635 mm", "pitch_mm": 3.2, "feed": "coax probe monopole, centre, 6.0 mm above substrate", "beam": "broadside RHCP"},
        "design": {"eta_design": eta, "taper_c": taper, "alpha_max_available_np_m": float(np.nanmax(table.alpha)), "M_t_cap": M_CAP, "n_table_cells": len(keys)},
        "active": {"eta_conv_surface_wave_to_aperture": eta_conv, "edge_fraction_surface_wave": 1.0 - eta_conv,
                   "pattern_of_aperture_current": pn, "gain_rhcp_upper_bound_dbi_launch1": pn["directivity_rhcp_dbi"] + 10 * np.log10(max(eta_conv, 1e-9)),
                   "by_radius": prof},
        "control": {"edge_fraction_surface_wave": 1.0, "note": "isotropic cells: first-order alpha = 0, only dielectric/conductor loss (alpha_d ~0.2, alpha_c ~0.4 Np/m); expected near-zero RHCP broadside gain; the radiation that remains is monopole direct radiation + edge diffraction (NOT modelled: NA)"},
        "NA": ["launch fraction of the monopole into the surface wave", "monopole direct radiation pattern", "edge diffraction pattern", "HFSS gain"],
    }
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1, default=float))
    print(json.dumps({k: pred["active"][k] for k in ("eta_conv_surface_wave_to_aperture", "edge_fraction_surface_wave", "gain_rhcp_upper_bound_dbi_launch1")}), pn)
    print(prof)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("grid_json")
    ap.add_argument("--eta", type=float, default=0.8)
    ap.add_argument("--taper", type=float, default=0.5)
    a = ap.parse_args()
    main(a.grid_json, a.eta, a.taper)
