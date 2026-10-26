"""Phase 20 R5-6 (INFRASTRUCTURE): HFSS export of a coded finite RIS board of CP-2 cells -- geometry, diode table, STEP (`Load_<k>` sheets), AEDT build script, README.

A board is `ny x nx` cells at the cell pitch (15 mm), each in one of the four one-hot states (`grid[j][i]` = 0..3, library state = index + 1; (i, j) = column, row; the board is centred on
the origin, x along columns, y along rows -- the same conventions as `forward.py` / `coding.py` / `supercell.build_supercell`). Every cell has four radial stubs (0/45/90/135 deg), each with a diode
span that is an EMPTY GAP in the exported metal bridged by one sheet `Load_<k>`, k = 4 * (j * nx + i) + stub + 1 (the order of `supercell.build_supercell(...).loads`). The diode of the
stub whose index equals the cell state is ON (R_s), the other three are OFF (C_T): one single-type lumped RLC boundary per sheet (no parallel/series question).

Outputs (`export_board`): `layout.csv` (cells), `diodes.csv` (every Load sheet: position, size, current-line angle, state), `hfss_build_board.py` (native AEDT script: zero-thickness metal sheets
at z = h, Perfect E, lumped RLC per diode, radiation box, plane-wave excitation, far-field sphere), `board.step` (CadQuery, gitignored), `README.md`. The prediction table of the README is a
placeholder unless `predictions` is passed: this module never invents numbers. MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.

The AEDT script reuses the call patterns of the Phase 21/23 board scripts that were run in AEDT (polylines, Perfect E, radiation, setup, far-field sphere). The two calls that are new here --
`AssignLumpedRLC` and the incident-wave assignment -- are verified in AEDT (3x3 smoke build, 2026-10-24): the build completed without errors, `Load_1`..`Load_36_RLC` carry R = 5.2 ohm (ON) / C = 0.025 pF (OFF) and the incident wave was assigned. The current-line direction of the boundaries is for the user to confirm. The circular-polarisation incidence (two plane waves 90 deg apart) remains a MANUAL step as documented in the README. Each call stays wrapped so a failure prints what to assign by hand.
"""

import csv
import io
import math
from dataclasses import dataclass
from pathlib import Path

LABELS = (
    "MODEL_ONLY",
    "NOT_HFSS_VALIDATED",
    "PROVISIONAL_DATASHEET",
    "VIA_OMITTED_APPROXIMATION",
    "AEDT_CALLS_VERIFIED_3X3_SMOKE_2026-10-24",
)
MAX_CELLS = 12  # HFSS size guard: 12 x 12 cells = 576 diodes; larger boards are refused until feasibility is agreed


@dataclass(frozen=True)
class Diode:
    k: int  # 1-based Load_<k>
    i: int
    j: int
    stub: int  # 0..3
    x_mm: float
    y_mm: float
    span_mm: float
    width_mm: float
    angle_deg: float
    on: bool


@dataclass(frozen=True)
class RisBoard:
    cell: object  # Cp2Cell (duck-typed: period, circumradius, apothem, overlap, pad, diode_span, line_at, width_at, stub_angles_deg)
    kit: object  # DiodeKit (on.R_s, off.C_T)
    grid: tuple  # ((state index 0..3, ...), ...) rows j, columns i
    name: str = "ris_board"
    f_ghz: float = 10.5
    theta_i_deg: float = 0.0
    phi_i_deg: float = 0.0
    handedness: str = "rhcp"

    @property
    def nx(self):
        return len(self.grid[0])

    @property
    def ny(self):
        return len(self.grid)

    @property
    def pitch_mm(self):
        return self.cell.period * 1e3

    @property
    def width_mm(self):
        return self.nx * self.pitch_mm

    @property
    def height_mm(self):
        return self.ny * self.pitch_mm

    @property
    def h_mm(self):
        return 1.575  # RT/duroid 5880, the characterised substrate


def make_board(cell, kit, grid, **kw):
    g = tuple(tuple(int(v) for v in row) for row in grid)
    if not g or not g[0] or any(len(r) != len(g[0]) for r in g):
        raise ValueError("grid must be a non-empty rectangular array")
    if any(v not in (0, 1, 2, 3) for r in g for v in r):
        raise ValueError("grid entries are state indices 0..3 (library state = index + 1)")
    if max(len(g), len(g[0])) > MAX_CELLS:
        raise ValueError(f"boards above {MAX_CELLS} x {MAX_CELLS} cells are refused (HFSS size not agreed)")
    return RisBoard(cell, kit, g, **kw)


def _centre_mm(b, i, j):
    return (i - (b.nx - 1) / 2) * b.pitch_mm, (j - (b.ny - 1) / 2) * b.pitch_mm


def _rect(cx, cy, length, width, ang):
    u, v = (math.cos(ang), math.sin(ang)), (-math.sin(ang), math.cos(ang))
    return [
        (cx + sl * length / 2 * u[0] + sw * width / 2 * v[0], cy + sl * length / 2 * u[1] + sw * width / 2 * v[1])
        for sl, sw in ((-1, -1), (1, -1), (1, 1), (-1, 1))
    ]


def cell_pieces_mm(cell):
    """Metal polygons of ONE cell centred at the origin, in mm: the octagon and, per stub, the inner piece (overlap + pad) and the outer piece (pad + line); the diode span between them is
    the gap. Returns (pieces, diode_sites) with pieces = [(label, [(x, y), ...])] and sites = [(x, y, span, width, angle_deg)] per stub."""
    c = cell
    a = c.apothem
    pieces = [
        (
            "Oct",
            [
                (
                    c.circumradius * 1e3 * math.cos(math.pi / 8 + 2 * math.pi * k / 8),
                    c.circumradius * 1e3 * math.sin(math.pi / 8 + 2 * math.pi * k / 8),
                )
                for k in range(8)
            ],
        )
    ]
    sites = []
    for s, ang_deg in enumerate(c.stub_angles_deg):
        ang = math.radians(ang_deg)
        r0, rd0 = a - c.overlap, a + c.pad
        rd1 = rd0 + c.diode_span
        r_end = rd1 + c.pad + c.line_at(ang_deg)
        w = c.width_at(ang_deg)
        mid = lambda r0_, r1_, a_=ang: (0.5 * (r0_ + r1_) * math.cos(a_) * 1e3, 0.5 * (r0_ + r1_) * math.sin(a_) * 1e3)  # noqa: E731
        pieces.append((f"S{s}in", _rect(*mid(r0, rd0), (rd0 - r0) * 1e3, w * 1e3, ang)))
        pieces.append((f"S{s}out", _rect(*mid(rd1, r_end), (r_end - rd1) * 1e3, w * 1e3, ang)))
        x, y = mid(rd0, rd1)
        sites.append((x, y, c.diode_span * 1e3, w * 1e3, ang_deg))
    return pieces, sites


def diodes(b):
    _, sites = cell_pieces_mm(b.cell)
    out, k = [], 0
    for j in range(b.ny):
        for i in range(b.nx):
            cx, cy = _centre_mm(b, i, j)
            for s, (x, y, span, w, ang) in enumerate(sites):
                k += 1
                out.append(Diode(k, i, j, s, round(cx + x, 6), round(cy + y, 6), span, w, ang, s == b.grid[j][i]))
    return out


def board_pieces_mm(b):
    """Every metal polygon of the board with a unique name `C<j>_<i>_<label>`, translated to the cell centres."""
    pieces, _ = cell_pieces_mm(b.cell)
    out = []
    for j in range(b.ny):
        for i in range(b.nx):
            cx, cy = _centre_mm(b, i, j)
            out += [(f"C{j}_{i}_{lab}", [(round(cx + x, 6), round(cy + y, 6)) for x, y in pts]) for lab, pts in pieces]
    return out


# ------------------------------------------------------------------------------------------------------------------------------------- files
def layout_csv_text(b):
    buf = io.StringIO(newline="")
    w = csv.writer(buf)
    w.writerow(["row_j", "col_i", "x_mm", "y_mm", "state_1to4", "diode_on_stub", "load_first", "load_last"])
    for j in range(b.ny):
        for i in range(b.nx):
            cx, cy = _centre_mm(b, i, j)
            k0 = 4 * (j * b.nx + i) + 1
            w.writerow([j, i, f"{cx:.6g}", f"{cy:.6g}", b.grid[j][i] + 1, b.grid[j][i], f"Load_{k0}", f"Load_{k0 + 3}"])
    return buf.getvalue()


def diodes_csv_text(b):
    buf = io.StringIO(newline="")
    w = csv.writer(buf)
    w.writerow(
        [
            "load",
            "row_j",
            "col_i",
            "stub",
            "x_mm",
            "y_mm",
            "span_mm",
            "width_mm",
            "current_line_angle_deg",
            "diode",
            "element",
        ]
    )
    for d in diodes(b):
        elem = f"R = {b.kit.on.R_s:g} ohm" if d.on else "C = C_off"
        w.writerow(
            [
                f"Load_{d.k}",
                d.j,
                d.i,
                d.stub,
                d.x_mm,
                d.y_mm,
                f"{d.span_mm:.6g}",
                f"{d.width_mm:.6g}",
                d.angle_deg,
                "ON" if d.on else "OFF",
                elem,
            ]
        )
    return buf.getvalue()


def export_step(b, path):
    """STEP of the board as ONE periodic-cell-style assembly (`Substrate`, `Metal` with the diode gaps, `Ground`, `Load_<k>`) via the existing CAD export of the supercell (CadQuery; slow and large for 10 x 10)."""
    from ..cad_export import export_unit_cell_step
    from ..ris.supercell import build_supercell

    sc = build_supercell(b.cell, [[b.cell.states()[b.grid[j][i]] for i in range(b.nx)] for j in range(b.ny)], b.kit)
    return export_unit_cell_step(sc, path)


_AEDT = r"""
# -*- coding: utf-8 -*-
# HFSS (AEDT) IronPython script -- RIS board "@@NAME@@": @@NX@@ x @@NY@@ CP-2 cells, @@NDIODE@@ diodes (generated by lwmom, Phase 20 R5-6 infrastructure).
# Run inside AEDT: Tools > Run Script...  MODEL_ONLY / NOT_HFSS_VALIDATED. The two calls marked VERIFIED ran in AEDT in the 3x3 smoke build of 2026-10-24 (current-line direction: to be confirmed
# by the user); if one fails the script prints what to assign by hand and continues.
import ScriptEnv
ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oProject = oDesktop.NewProject()
oProject.InsertDesign("HFSS", "@@NAME@@", "DrivenModal", "")
oDesign = oProject.SetActiveDesign("@@NAME@@")
oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.SetModelUnits(["NAME:Units Parameter", "Units:=", "mm", "Rescale:=", False])

H = "@@H@@mm"
WX = "@@W@@mm"; WY = "@@HT@@mm"
F_GHZ = @@F_GHZ@@
THETA_I = @@THETA_I@@; PHI_I = @@PHI_I@@          # direction of the SPECULAR beam = tangential wave vector of the incident wave (see README)
R_ON = "@@R_ON@@ohm"; C_OFF = "@@C_OFF@@pF"
# pieces: (name, [(x, y), ...] mm)
PIECES = [
@@PIECES@@
]
# diodes: (load name, state ON?, [x0, y0, x1, y1] current line start/end mm, [(x, y) x4] sheet corners mm)
LOADS = [
@@LOADS@@
]

oProject.GetDefinitionManager().AddMaterial(["NAME:RT5880_board", "CoordinateSystemType:=", "Cartesian",
    "BulkOrSurfaceType:=", 1, ["NAME:PhysicsTypes", "set:=", ["Electromagnetic"]],
    "permittivity:=", "@@EPS_R@@", "dielectric_loss_tangent:=", "@@TAND@@"])

def attrs(name, material, solve_inside, color="(255 128 0)", transp=0):
    return ["NAME:Attributes", "Name:=", name, "Flags:=", "", "Color:=", color, "Transparency:=", transp,
            "PartCoordinateSystem:=", "Global", "UDMId:=", "", "MaterialValue:=", '"' + material + '"',
            "SurfaceMaterialValue:=", '""', "SolveInside:=", solve_inside, "IsMaterialEditable:=", True,
            "UseMaterialAppearance:=", False, "IsLightweight:=", False]

def mm(v):
    return repr(float(v)) + "mm"

def poly_sheet(name, pts, z, color="(255 128 0)"):
    pl = ["NAME:PolylineParameters", "IsPolylineCovered:=", True, "IsPolylineClosed:=", True,
          ["NAME:PolylinePoints"] + [["NAME:PLPoint", "X:=", mm(x), "Y:=", mm(y), "Z:=", z] for (x, y) in pts + [pts[0]]],
          ["NAME:PolylineSegments"] + [["NAME:PLSegment", "SegmentType:=", "Line", "StartIndex:=", k, "NoOfPoints:=", 2] for k in range(len(pts))],
          ["NAME:PolylineXSection", "XSectionType:=", "None", "XSectionOrient:=", "Auto", "XSectionWidth:=", "0mm",
           "XSectionTopWidth:=", "0mm", "XSectionHeight:=", "0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=", "Corner"]]
    oEditor.CreatePolyline(pl, attrs(name, "vacuum", True, color))

def unite(names):
    oEditor.Unite(["NAME:Selections", "Selections:=", ",".join(names)], ["NAME:UniteParameters", "KeepOriginals:=", False])

# substrate and ground (centred on the origin)
oEditor.CreateBox(["NAME:BoxParameters", "XPosition:=", "-" + WX + "/2", "YPosition:=", "-" + WY + "/2", "ZPosition:=", "0mm",
    "XSize:=", WX, "YSize:=", WY, "ZSize:=", H], attrs("Substrate", "RT5880_board", True, "(0 128 0)", 0.6))
oEditor.CreateRectangle(["NAME:RectangleParameters", "IsCovered:=", True, "XStart:=", "-" + WX + "/2", "YStart:=", "-" + WY + "/2",
    "ZStart:=", "0mm", "Width:=", WX, "Height:=", WY, "WhichAxis:=", "Z"], attrs("Ground", "vacuum", True, "(128 128 128)"))

# metal: every piece a zero-thickness sheet at z = h, united in blocks of 100 into "Metal"
names = []
for nm, pts in PIECES:
    poly_sheet(nm, pts, H)
    names.append(nm)
blocks = []
for b in range(0, len(names), 100):
    grp = names[b:b + 100]
    if len(grp) > 1:
        unite(grp)
    blocks.append(grp[0])
if len(blocks) > 1:
    unite(blocks)
oEditor.ChangeProperty(["NAME:AllTabs", ["NAME:Geometry3DAttributeTab", ["NAME:PropServers", names[0]],
    ["NAME:ChangedProps", ["NAME:Name", "Value:=", "Metal"]]]])

# diode sheets Load_<k> bridging the gaps
for nm, on, line, corners in LOADS:
    poly_sheet(nm, corners, H, "(255 0 0)")

oBnd = oDesign.GetModule("BoundarySetup")
oBnd.AssignPerfectE(["NAME:PerfE_Ground", "Objects:=", ["Ground"], "InfGroundPlane:=", False])
oBnd.AssignPerfectE(["NAME:PerfE_Metal", "Objects:=", ["Metal"], "InfGroundPlane:=", False])

# VERIFIED in AEDT (3x3 smoke build, 2026-10-24): one single-type Lumped RLC boundary per diode sheet, current line along the stub
failed = []
for nm, on, line, corners in LOADS:
    try:
        rlc = ["NAME:" + nm + "_RLC", "Objects:=", [nm],
               ["NAME:CurrentLine", "Start:=", [mm(line[0]), mm(line[1]), H], "End:=", [mm(line[2]), mm(line[3]), H]],
               "RLC Type:=", "Parallel",
               "UseResist:=", bool(on), "Resistance:=", R_ON if on else "0ohm",
               "UseInduct:=", False, "Inductance:=", "0nH",
               "UseCap:=", (not on), "Capacitance:=", "0pF" if on else C_OFF]
        oBnd.AssignLumpedRLC(rlc)
    except Exception as e:
        failed.append(nm)
if failed:
    oDesktop.AddMessage(oProject.GetName(), "@@NAME@@", 2, "Lumped RLC assignment failed for %d sheets (first: %s). Assign by hand: ON sheets (diodes.csv) R = %s, OFF sheets C = %s, one element type each, current line along the stub (diodes.csv, current_line_angle_deg)." % (len(failed), failed[0], R_ON, C_OFF))

# air box + radiation
PAD = "@@PAD@@mm"
oEditor.CreateBox(["NAME:BoxParameters", "XPosition:=", "-" + WX + "/2-" + PAD, "YPosition:=", "-" + WY + "/2-" + PAD, "ZPosition:=", "-" + PAD,
    "XSize:=", WX + "+2*" + PAD, "YSize:=", WY + "+2*" + PAD, "ZSize:=", PAD + "+" + H + "+" + PAD + "+@@PAD_TOP@@mm"],
    attrs("Air", "air", True, "(143, 175, 143)", 0.9))
oBnd.AssignRadiation(["NAME:Rad1", "Objects:=", ["Air"], "IsFssReference:=", False, "IsForPML:=", False])

# VERIFIED in AEDT (3x3 smoke build, 2026-10-24): incident plane wave. @@POL_NOTE@@
try:
    oBnd.AssignPlaneWave(["NAME:IncWave1", "IncludePortPostProcessing:=", False, "SpecifyAsPlaneWave:=", True,
        "IncomingDirection:=", ["PhiAngle:=", str(PHI_I + 180.0) + "deg", "ThetaAngle:=", str(THETA_I) + "deg"],
        "Polarization:=", "Vertical"])
except Exception as e:
    oDesktop.AddMessage(oProject.GetName(), "@@NAME@@", 2, "Incident wave not assigned by script (%s): add it by hand as in the README (circular polarisation = two plane waves)." % str(e))

oAna = oDesign.GetModule("AnalysisSetup")
oAna.InsertSetup("HfssDriven", ["NAME:Setup1", "Frequency:=", str(F_GHZ) + "GHz", "MaxDeltaS:=", 0.02,
    "MaximumPasses:=", 12, "MinimumPasses:=", 2, "MinimumConvergedPasses:=", 2, "PercentRefinement:=", 30,
    "IsEnabled:=", True, "BasisOrder:=", 1, "DoLambdaRefine:=", True, "DoMaterialLambda:=", True,
    "SetLambdaTarget:=", False, "Target:=", 0.3333, "UseMaxTetIncrease:=", False])
oRad = oDesign.GetModule("RadField")
oRad.InsertFarFieldSphereSetup(["NAME:3D", "UseCustomRadiationSurface:=", False,
    "ThetaStart:=", "0deg", "ThetaStop:=", "90deg", "ThetaStep:=", "0.5deg",
    "PhiStart:=", "0deg", "PhiStop:=", "360deg", "PhiStep:=", "1deg", "UseLocalCS:=", False])
oEditor.FitAll()
oProject.SaveAs(oDesktop.GetProjectDirectory() + "/@@NAME@@.aedt", True)
oDesktop.AddMessage(oProject.GetName(), "@@NAME@@", 0, "Built: @@NX@@ x @@NY@@ cells, @@NDIODE@@ diode sheets. Check the RLC boundaries and the incident wave, validate, then Analyze.")
""".lstrip("\n")


def aedt_script_text(b):
    sub = {"eps_r": 2.2, "tan_d": 0.0009}
    pieces = ",\n".join(
        f"  ({nm!r}, {[(round(x, 6), round(y, 6)) for x, y in pts]!r})" for nm, pts in board_pieces_mm(b)
    )
    loads = []
    for d in diodes(b):
        a = math.radians(d.angle_deg)
        ux, uy = math.cos(a), math.sin(a)
        line = [
            d.x_mm - ux * d.span_mm / 2,
            d.y_mm - uy * d.span_mm / 2,
            d.x_mm + ux * d.span_mm / 2,
            d.y_mm + uy * d.span_mm / 2,
        ]
        corners = [
            (round(p[0] + d.x_mm, 6), round(p[1] + d.y_mm, 6)) for p in _rect(0.0, 0.0, d.span_mm, d.width_mm, a)
        ]
        loads.append(f"  ('Load_{d.k}', {d.on!r}, {[round(v, 6) for v in line]!r}, {corners!r})")
    on = b.kit.on.R_s
    c_off = b.kit.off.C_T * 1e12
    values = {
        "NAME": "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in b.name),
        "NX": str(b.nx),
        "NY": str(b.ny),
        "NDIODE": str(4 * b.nx * b.ny),
        "H": repr(b.h_mm),
        "W": repr(round(b.width_mm, 6)),
        "HT": repr(round(b.height_mm, 6)),
        "F_GHZ": repr(b.f_ghz),
        "THETA_I": repr(float(b.theta_i_deg)),
        "PHI_I": repr(float(b.phi_i_deg)),
        "R_ON": f"{on:g}",
        "C_OFF": f"{c_off:g}",
        "PIECES": pieces,
        "LOADS": ",\n".join(loads),
        "EPS_R": repr(sub["eps_r"]),
        "TAND": repr(sub["tan_d"]),
        "PAD": "15",
        "PAD_TOP": "30",
        "POL_NOTE": "RHCP in the propagation-relative sense = two plane waves, theta-polarised and phi-polarised, 90 deg apart in time phase (README); the script assigns ONE linear wave as a placeholder.",
    }
    text = _AEDT
    for k, v in values.items():
        text = text.replace(f"@@{k}@@", v)
    return text


def readme_text(b, predictions=None):
    n = b.nx * b.ny
    ons = [sum(1 for row in b.grid for v in row if v == s) for s in range(4)]
    pred = (
        "**PREDICTIONS_NOT_YET_GENERATED.** The forward-model predictions for this coding (main-beam gain, pointing, SLL, opposite-handed specular lobe, cuts) are written **before** the HFSS run in a separate commit and compared with "
        "`phase20_r5_hfss_compare.py` (comparison only). Until then this package has no `mom_reference.csv` and no number to compare with."
        if predictions is None
        else predictions
    )
    return f"""# {b.name} -- {b.nx} x {b.ny} CP-2 RIS board for HFSS (Phase 20 R5-6)

**State: INFRASTRUCTURE / MODEL_ONLY / NOT_HFSS_VALIDATED.** Labels: {", ".join(LABELS)}. {pred}

Board: {b.width_mm:g} x {b.height_mm:g} mm, {n} cells ({b.pitch_mm:g} mm pitch), {4 * n} diodes (MADP-000907-14020, PROVISIONAL_DATASHEET: ON R_s = {b.kit.on.R_s:g} ohm, OFF C_T = {b.kit.off.C_T * 1e12:g} pF, L_s = 0 assumed, flip-chip pads not modelled),
RT/duroid 5880 eps_r 2.2, tan d 0.0009, h = {b.h_mm:g} mm, PEC ground, via omitted. State counts (states 1-4): {ons[0]} / {ons[1]} / {ons[2]} / {ons[3]}. Design frequency {b.f_ghz:g} GHz.

## Files
`hfss_build_board.py` (AEDT script: Tools > Run Script), `layout.csv` (one row per cell: position, state 1-4), `diodes.csv` (one row per `Load_<k>`: position, size, current-line angle, ON/OFF, element), `board.step` (STEP with `Substrate`, `Metal`, `Ground`, `Load_<k>`; gitignored, regenerate with the exporter).

## Conventions
* x along columns, y along rows, board centred on the origin; cell (row j, column i) has state s = layout.csv `state_1to4`; its stub `s-1` (angle 45 (s-1) deg) has the diode ON, the other three OFF.
* `Load_<k>`: k = 4 (j nx + i) + stub + 1. Each is a sheet bridging the 0.75 mm diode gap; the current line runs ALONG the stub (angle in `diodes.csv`). **One single-type Lumped RLC boundary per sheet**: ON = R only, OFF = C only (L_s = 0).
* Metal is zero-thickness (sheets at z = h, Perfect E), like the MoM; the `board.step` metal is the 35 um slab (state which you used).
* Incident wave: propagation-relative **RHCP** at the tangential wave vector of the specular direction (theta_i, phi_i) = ({b.theta_i_deg:g}, {b.phi_i_deg:g}) deg. In HFSS: two plane waves of equal amplitude, theta- and phi-polarised, 90 deg apart in phase
  with the sign that gives right-handed rotation about the propagation direction (check with the bare grounded slab: the reflected wave must be left-handed there, `docs/theory/CONVENTIONS.md`). The script assigns one linear wave as a placeholder.
* Radiation box padded by 15 mm; far-field sphere theta 0-90 deg, phi 0-360 deg. Report the scattered far field: co-handed (RHCP) and opposite-handed (LHCP) gain, normalised to the power incident on the board area.

## Honest limits
Model uncertainty ~10 deg on the cell phases (staircased octagon in the MoM), the local-response assumption of the forward model (checked at array level by the R5-4 server runs), no via, no bias network, lumped diode without pad geometry, finite 150 mm board edge effects in HFSS that the aperture model does not have.
The two new AEDT calls (`AssignLumpedRLC`, `AssignPlaneWave`) are verified in AEDT (3x3 smoke build, 2026-10-24); the current-line direction is for the user to confirm, and the circular-polarisation incidence (two waves) stays a manual step.
"""


def export_board(out_dir, b, with_step=True, predictions=None):
    """Write the package; returns the list of file names. No prediction file is written unless `predictions` text is given (it only changes the README)."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "layout.csv").write_text(layout_csv_text(b), encoding="utf-8", newline="")
    (out / "diodes.csv").write_text(diodes_csv_text(b), encoding="utf-8", newline="")
    (out / "hfss_build_board.py").write_text(aedt_script_text(b), encoding="utf-8", newline="\n")
    (out / "README.md").write_text(readme_text(b, predictions), encoding="utf-8", newline="\n")
    names = ["layout.csv", "diodes.csv", "hfss_build_board.py", "README.md"]
    if with_step:
        export_step(b, out / "board.step")
        names.append("board.step")
    return names
