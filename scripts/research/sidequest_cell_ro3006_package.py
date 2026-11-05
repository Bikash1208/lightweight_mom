"""HFSS unit-cell package: RO3006 1.27 mm, disc r 2.2 (rim 0.2, cross slit 0.4 x 4.0), two diode sheets (DX at (0, +0.8), DY at (+0.8, 0), 0.4 x 1.0, both diodes on the line = scalar cell), lattice 6.4 mm.
OFF = diode sheets unassigned (open), ON = ONE Impedance boundary on both sheets, Rs = R_on * W / L = 2.44 * 1.0 / 0.4 = 6.1 ohm/sq, reactance 0.  Eigenmode at 90 deg (phase along the wave), PEC lid 30 mm.
Predictions (MODEL_ONLY): loaded Bloch tensor (MoM mesh 64) + transverse resonance with the lid, as for the RT5880 cell (HFSS-checked there: model 10.662 / 11.007 against HFSS 10.670 / 11.006 GHz).
    python sidequest_cell_ro3006_package.py [geom|pred|all]  ->  hfss_exports/sidequest_cell_RO3006_r2p2/  (+ sync C:\\hmsa_sync\\data\\sidequest\\cell_RO3006_r2p2)"""
import csv
import json
import shutil
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import cadquery as cq
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from sidequest_board_k2_lp_package import read_names
from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap

OUT = ROOT / "hfss_exports/sidequest_cell_RO3006_r2p2"
SYNC = Path(r"C:\hmsa_sync\data\sidequest\cell_RO3006_r2p2")
P, H, T, GT, LID = 6.4, 1.27, 0.035, 0.035, 30.0
R, G, WD, POS = 2.2, 0.4, 1.0, 0.8
LS = R - 0.2
ER, TAND = 6.15, 0.002
RS = 2.44 * WD / G
C0, EPS0 = 299792458.0, 8.8541878128e-12
BETA = np.radians(90.0) / (P * 1e-3)
FG = np.array([7.6, 8.0, 8.4, 8.8, 9.2, 9.6, 10.0, 10.4, 10.8, 11.2]) * 1e9


def box(cx, cy, lx, ly, z0, z1):
    return cq.Workplane("XY").workplane(offset=z0).center(cx, cy).rect(lx, ly).extrude(z1 - z0)


def sheet(cx, cy, lx, ly):
    z = T / 2
    pts = [cq.Vector(cx - lx / 2, cy - ly / 2, z), cq.Vector(cx + lx / 2, cy - ly / 2, z), cq.Vector(cx + lx / 2, cy + ly / 2, z), cq.Vector(cx - lx / 2, cy + ly / 2, z)]
    return cq.Face.makeFromWires(cq.Wire.makePolygon(pts, close=True))


def build():
    patch = cq.Workplane("XY").circle(R).extrude(T).cut(box(0, 0, G, 2 * LS, 0, T)).cut(box(0, 0, 2 * LS, G, 0, T))
    ground = box(0, 0, P, P, -H - GT, -H)
    metal = cq.Compound.makeCompound([patch.val(), ground.val()])
    return {"METAL_PEC": metal, "SUBSTRATE": box(0, 0, P, P, -H, 0).val(), "AIR": box(0, 0, P, P, -H - GT, LID).val(), "DX_SHEET": sheet(0, POS, G, WD), "DY_SHEET": sheet(POS, 0, WD, G)}


def intended():
    xy = {0.0, P / 2, R, G / 2, LS, POS - WD / 2, POS + WD / 2, WD / 2}
    xy |= {-v for v in xy}
    z = {0.0, T, T / 2, -H, -H - GT, LID}
    return dict(xy=np.array(sorted(round(v, 6) for v in xy)), z=np.array(sorted(round(v, 6) for v in z)), dir=np.array([-1.0, 0.0, 1.0]))


def checks(path):
    res, ok_all = {}, True

    def rec(n, ok, **kw):
        nonlocal ok_all
        ok_all &= bool(ok)
        res[n] = {"ok": bool(ok), **kw}

    names = read_names(path)
    rec("bodies_and_names", set(names) == {"METAL_PEC", "SUBSTRATE", "AIR", "DX_SHEET", "DY_SHEET"}, names=sorted(names))
    solids = names["METAL_PEC"].Solids()
    patch = [s for s in solids if s.BoundingBox().zmin > -1e-6]
    gnd = [s for s in solids if s.BoundingBox().zmax < -1.2]
    rec("metal_one_part_two_lumps", len(patch) == 1 and len(gnd) == 1, patch_lumps=len(patch), ground_lumps=len(gnd))
    exp_area = np.pi * R**2 - (2 * (2 * LS) * G - G * G)
    pb = patch[0].BoundingBox()
    rec("patch_disc_minus_cross_slit", abs(patch[0].Volume() - exp_area * T) < 2e-4 and abs(pb.xlen - 2 * R) < 1e-6, volume=patch[0].Volume(), expected=exp_area * T, radius_mm=pb.xlen / 2, slit_width_mm=G, slit_length_mm=2 * LS, rim_mm=R - LS)
    gb = gnd[0].BoundingBox()
    rec("ground", abs(gb.zmin - (-H - GT)) < 1e-9 and abs(gb.zmax + H) < 1e-9 and abs(gb.xlen - P) < 1e-9, z=[gb.zmin, gb.zmax], xlen=gb.xlen)
    sub = names["SUBSTRATE"]
    sb = sub.BoundingBox()
    rec("substrate", abs(sb.zmin + H) < 1e-9 and abs(sb.zmax) < 1e-9 and abs(sb.xlen - P) < 1e-9 and abs(sb.ylen - P) < 1e-9, thickness_mm=sb.zlen)
    c1 = sub.intersect(cq.Compound.makeCompound(solids)).Volume()
    rec("no_intersections", c1 < 1e-9, substrate_x_metal=c1)
    dx, dy = names["DX_SHEET"].BoundingBox(), names["DY_SHEET"].BoundingBox()
    okx = abs(dx.xlen - G) < 1e-9 and abs(dx.ylen - WD) < 1e-9 and abs(dx.ymin - (POS - WD / 2)) < 1e-9 and abs((dx.xmin + dx.xmax) / 2) < 1e-9 and abs(dx.zmin - T / 2) < 1e-9
    oky = abs(dy.ylen - G) < 1e-9 and abs(dy.xlen - WD) < 1e-9 and abs(dy.xmin - (POS - WD / 2)) < 1e-9 and abs((dy.ymin + dy.ymax) / 2) < 1e-9 and abs(dy.zmin - T / 2) < 1e-9
    metal = patch[0]
    corners = [cq.Vertex.makeVertex(sx * G / 2, y, T / 2) for sx in (-1, 1) for y in (POS - WD / 2, POS + WD / 2)] + [cq.Vertex.makeVertex(x, sy * G / 2, T / 2) for sy in (-1, 1) for x in (POS - WD / 2, POS + WD / 2)]
    dist = max(metal.distance(v) for v in corners)
    clear = float(np.hypot(max(dx.xmin - dy.xmax, dy.xmin - dx.xmax, 0), max(dx.ymin - dy.ymax, dy.ymin - dx.ymax, 0)))
    rec("sheets", okx and oky and dist < 1e-6 and clear >= 0.1, DX="0.4 x 1.0 at (0, 0.8)", DY="1.0 x 0.4 at (0.8, 0)", z=T / 2, corner_to_metal_max_mm=float(dist), DX_DY_clearance_mm=round(clear, 4), W_over_L=WD / G, Rs_ohm_per_sq=2.44 * dy.xlen / dx.xlen)
    ab = names["AIR"].BoundingBox()
    rec("air_box_ground_to_lid_single_face_walls", len(names["AIR"].Faces()) == 6 and abs(ab.zmin + H + GT) < 1e-9 and abs(ab.zmax - LID) < 1e-9 and abs(ab.xlen - P) < 1e-9 and abs(ab.ylen - P) < 1e-9, faces=len(names["AIR"].Faces()), z=[ab.zmin, ab.zmax])
    return res, ok_all


def geom():
    OUT.mkdir(parents=True, exist_ok=True)
    parts = build()
    asm = cq.Assembly()
    cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "AIR": (0.7, 0.8, 1.0), "DX_SHEET": (0.9, 0.1, 0.1), "DY_SHEET": (0.1, 0.1, 0.9)}
    for nm, shp in parts.items():
        asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
    path = OUT / "cell_RO3006_r2p2.step"
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    it = intended()
    snapped, stats = snap_step_text(raw, it, description="sidequest RO3006 r2.2 diode cell, mm, 0.001 mm grid")
    path.write_text(snapped, encoding="utf-8", newline="\n")
    bad = inexact_literals(path.read_text(encoding="utf-8"), it)
    vs = verify_snap(raw, snapped)
    res, ok_all = checks(path)
    exact = {"ok": len(bad) == 0 and vs.get("collapsed_edges", 1) == 0 and not vs.get("distinct_merge_violations"), "inexact_literals": len(bad), "n_snapped": stats["n_snapped"], "max_snap_mm": stats["max_snap_mm"], "grid_mm": 0.001}
    res["EXACT_DIMS"] = exact
    ok_all &= exact["ok"]
    (OUT / "geometry_checks.json").write_text(json.dumps({"label": "MODEL_ONLY geometry checks on the re-imported STEP", "all_pass": bool(ok_all), "checks": res}, indent=1, default=float), encoding="utf-8", newline="\n")
    print({k: v["ok"] for k, v in res.items()}, "ALL", ok_all, flush=True)
    if not ok_all:
        sys.exit("checks failed")


# ---------------------------------------------------------------- predictions
def y_air_lid(beta, om, k0, L=LID * 1e-3):
    kz = np.sqrt(k0**2 - beta**2 + 0j)
    if kz.imag > 0:
        kz = -kz
    return -1j * om * EPS0 / (kz * np.tan(kz * L))


def stack_b(beta, f, er, h):
    from src.dispersion import sheet_model as sm
    from src.materials import Substrate
    sub = Substrate("s", er, TAND, h * 1e-3)
    om, k0 = 2 * np.pi * f, 2 * np.pi * f / C0
    return (y_air_lid(beta, om, k0) + sm._y_slab(beta, om, sub.eps_complex, k0, sub.thickness)).imag


def pjob(a):
    import warnings as w
    w.simplefilter("ignore")
    import sidequest_cell_search as C
    import sidequest_gate_cell as gc
    import sidequest_k2_default_lib as L
    from src.materials import Substrate
    r, er, mode, f = a
    gc.SUB = Substrate("ro3006", er, TAND, H * 1e-3)
    L.R = r
    return C.bloch_loaded(L.cell(round(r - 0.2, 3), POS, mode), f, 64, BETA, 0.0)


def predictions():
    cases = {"nominal": (R, ER), "r-0.05": (R - 0.05, ER), "r+0.05": (R + 0.05, ER), "er-0.15": (R, ER - 0.15), "er+0.15": (R, ER + 0.15)}
    jobs = [(r, er, m, f) for r, er in cases.values() for m in ("open", "ON") for f in FG]
    with ProcessPoolExecutor(30) as ex:
        res = list(ex.map(pjob, jobs))
    it = iter(res)
    out = {"label": "MODEL_ONLY: loaded Bloch tensor (mesh 64), transverse resonance with the 30 mm PEC lid, 90 deg phase along x (beta = pi/2 / 6.4 mm); OFF = open slit, ON = 2.44 ohm only (Rs 6.1 ohm/sq on the 0.4 x 1.0 sheets)", "phase_deg": 90, "rows": []}
    for name, (r, er) in cases.items():
        fr = {}
        for mode in ("open", "ON"):
            yt = np.array([next(it) for _ in FG])
            cs = CubicSpline(FG, yt[:, 0, 0].imag)
            try:
                fr[mode] = brentq(lambda f, cs=cs, er=er: stack_b(BETA, f, er, H) + float(cs(f)), FG[0], FG[-1], xtol=1e3)
            except ValueError:
                fr[mode] = None
        out["rows"].append({"case": name, "r_mm": r, "eps_r": er, "f_OFF_GHz": None if fr["open"] is None else round(fr["open"] / 1e9, 4), "f_ON_GHz": None if fr["ON"] is None else round(fr["ON"] / 1e9, 4)})
    for rw in out["rows"]:
        if rw["f_OFF_GHz"] and rw["f_ON_GHz"]:
            rw["shift_ON_minus_OFF_GHz"] = round(rw["f_ON_GHz"] - rw["f_OFF_GHz"], 4)
    fs_ = np.linspace(3e9, 11.7e9, 400); vs_ = np.array([stack_b(BETA, f, ER, H) for f in fs_]); ix = np.where(np.diff(np.sign(vs_)) != 0)[0]
    out["bare_substrate_TM0_GHz_with_lid"] = round(brentq(lambda f: stack_b(BETA, f, ER, H), fs_[ix[0]], fs_[ix[0] + 1]) / 1e9, 4) if len(ix) else None
    out["note"] = "the RT5880 version of this pipeline matched HFSS to 0.01 GHz (model 10.662 / 11.007, HFSS 10.670 / 11.006); the RO3006 cell is not HFSS-checked. +-0.15 GHz is the stated model uncertainty; nominal-case sensitivities are in the other rows."
    return out


def docs(pred):
    row = pred["rows"][0]
    sens = "; ".join("%s: %s / %s" % (r["case"], r["f_OFF_GHz"], r["f_ON_GHz"]) for r in pred["rows"][1:])
    setup = f"""# HFSS: RO3006 1.27 mm diode cell, disc r 2.2 mm, eigenmode OFF/ON at 90 deg. EXPORT STANDARD. Predictions MODEL_ONLY (`predictions.json`).

**1. Assign material `pec` to `METAL_PEC` (patch + ground = ONE named part, 2 DC-isolated lumps; the STEP imports as vacuum).** `SUBSTRATE`: RO3006, eps_r {ER}, tan d {TAND}, z -{H}..0 (assign by name).
**2. `AIR`** = ONE vacuum box z = {-H-GT:.3f} ... +{LID:g} mm (ground bottom to the lid), {P} x {P}, ONE face per side wall; it overlaps `SUBSTRATE` and `METAL_PEC` on purpose and is FIRST in the tree. **Primary/Secondary (Lattice Pair) boundaries on the four side faces of `AIR`** (they cover the substrate region too). **3. PEC lid:** Perfect E on the top face of `AIR` (z = +{LID:g}).
**Files:** `cell_RO3006_r2p2.step` (`METAL_PEC`, `SUBSTRATE`, `AIR`, `DX_SHEET`, `DY_SHEET`), `states.csv`, `predictions.json`, `geometry_checks.json` (all PASS on the re-imported STEP, EXACT_DIMS PASS), `results_template.csv`.
**Geometry** (mm): disc r {R} (35 um copper, z 0..0.035) with a cross slit {G} wide and {2*LS} long (rim {R-LS:.1f} at the four ends); ground 35 um (z {-H-GT:.3f}..{-H}); `DX_SHEET` {G} x {WD} at (0, +{POS}), `DY_SHEET` {WD} x {G} at (+{POS}, 0), at z = {T/2}; the sheet edges lie on the metal faces across the slit (current along the {G} mm gap, width {WD} mm), the two sheets are 0.14 mm apart.
**4. States:** OFF = **no boundary on the sheets** (diode OFF = open). ON = **ONE Impedance boundary on `DX_SHEET` and `DY_SHEET`: Rs = {RS:g} ohm/sq, Reactance 0** (Rs = R_on x W / L = 2.44 ohm x {WD} mm / {G} mm; **W/L = {WD/G:g}**). No RLC.
**5. Eigenmode:** Lattice Pair 1 (x faces, Phase Delay `$ph_x`), Pair 2 (y faces, `$ph_y`); x run `$ph_x = 90deg`, `$ph_y = 0deg`; y run swapped (the cell is x/y symmetric: expect the same values). Min frequency 5 GHz, 8 modes, delta f 0.1 %; mesh: substrate 0.8 mm, patch 0.2 mm, sheets 0.1 mm. Clean Up Solutions before each run. The bound mode is the lowest TM0-like mode below the light line (11.711 GHz at 90 deg); identify it by the Ez fraction.
**Record** (`results_template.csv`): f of that mode for OFF and ON, x and y.
**Prediction at 90 deg (MODEL_ONLY, +-0.15 GHz):** OFF {row['f_OFF_GHz']} GHz, ON {row['f_ON_GHz']} GHz (shift {row.get('shift_ON_minus_OFF_GHz')} GHz); bare substrate TM0 with the lid {pred['bare_substrate_TM0_GHz_with_lid']} GHz. Sensitivity OFF / ON [GHz]: {sens}.
**Falsifier:** ON must lie ABOVE OFF (the diode shorts the slit and lowers Y, as in the RT5880 cell: 10.670 -> 11.006 GHz); a shift below 0.1 GHz or ON below OFF means the cell does not switch as modelled; report both values. A board for this cell is built only after this check passes.
"""
    (OUT / "HFSS_SETUP.md").write_text(setup, encoding="utf-8", newline="\n")
    with open(OUT / "states.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["state", "DX_SHEET", "DY_SHEET", "note"])
        w.writerow(["OFF", "no boundary", "no boundary", "diode OFF = open"])
        w.writerow(["ON", f"Impedance boundary Rs = {RS:g} ohm/sq, X = 0", f"same boundary (ONE boundary on both sheets)", f"Rs = R_on * W / L = 2.44 * {WD} / {G}; W/L = {WD/G:g}"])
    with open(OUT / "results_template.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write("state,run,f_GHz,mode_is_TM0_below_light_line,Ez_fraction,notes\nOFF,x,,,,\nON,x,,,,\nOFF,y,,,,\nON,y,,,,\n")
    SYNC.mkdir(parents=True, exist_ok=True)
    for f in OUT.iterdir():
        shutil.copy2(f, SYNC / f.name)
    print("synced", sorted(p.name for p in SYNC.iterdir()))


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("geom", "all"):
        geom()
    if what in ("pred", "all"):
        pr = predictions()
        (OUT / "predictions.json").write_text(json.dumps(pr, indent=1, default=float), encoding="utf-8", newline="\n")
        print(json.dumps(pr["rows"], indent=0), pr["bare_substrate_TM0_GHz_with_lid"], flush=True)
        docs(pr)
