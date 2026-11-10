"""STEP + manual eigenmode packages; no AEDT scripts. MODEL_ONLY."""
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(Path(__file__).resolve().parent))
import cadquery as cq
from OCP.BRepCheck import BRepCheck_Analyzer
from sidequest_board_k2_lp_package import read_names
from sub6_step_exact import exactify_step_text,audit_step_text
from src.transmitarray.export import verify_snap

OUT = ROOT/"hfss_exports/sub6_island_baseline_g1"
P,H,T,LID = 11.6,3.175,.035,150.
VIAS = [(-1.780,-.530),(-3.630,-3.630)]
CELLS = {
    "BASELINE":dict(sl=4.348,sw=1.450,il=2.174,bot=-.850,top=.450),
    "G1":dict(sl=4.592,sw=1.329,il=1.933,bot=-.725,top=.329),
}


def box(x0,x1,y0,y1,z0,z1):
    return cq.Workplane("XY").workplane(offset=z0).center((x0+x1)/2,(y0+y1)/2).rect(x1-x0,y1-y0).extrude(z1-z0).val()


def cylinder(x,y,r,z0,z1):
    return cq.Workplane("XY").workplane(offset=z0).center(x,y).circle(r).extrude(z1-z0).val()


def rotated(shape):
    return shape.rotate((0,0,0),(0,0,1),90.)


def build(cell,via):
    c = CELLS[cell]
    patch = cylinder(0,0,5.435,0,T).cut(box(-c["sl"],c["sl"],-c["sw"],c["sw"],-.01,T+.01))
    island = box(-c["il"],c["il"],c["bot"],c["top"],0,T)
    ground = box(-P/2,P/2,-P/2,P/2,-H-T,-H)
    sub = box(-P/2,P/2,-P/2,P/2,-H,0)
    fr4 = box(-P/2,P/2,-P/2,P/2,-4.010,-3.210) if via else None
    if via:
        metals = [island,patch]
        for index,(x,y) in enumerate(VIAS):
            barrel = cylinder(x,y,.150,-4.045,0.)
            pad = cylinder(x,y,.300,-4.045,-4.010)
            stub = box(x,x+2.,y-.150,y+.150,-4.045,-4.010)
            metals[index] = metals[index].fuse(barrel,pad,stub).clean()
            ground = ground.cut(cylinder(x,y,.300,-3.220,-3.165))
            sub = sub.cut(cylinder(x,y,.150,-3.185,.010))
            fr4 = fr4.cut(cylinder(x,y,.150,-4.020,-3.200))
        island,patch = metals
    bottom = -34.045 if via else -33.210
    air = box(-P/2,P/2,-P/2,P/2,bottom,LID)
    vertices = [cq.Vector(-.175,c["top"],T/2),cq.Vector(.175,c["top"],T/2),
                cq.Vector(.175,c["sw"],T/2),cq.Vector(-.175,c["sw"],T/2)]
    sheet = cq.Face.makeFromWires(cq.Wire.makePolygon(vertices,close=True))
    parts = dict(METAL_PEC=cq.Compound.makeCompound([patch,island,ground]),
                 SUBSTRATE=sub,AIR=air,DIODE_SHEET=sheet)
    if via:
        parts["SUBSTRATE_FR4"] = fr4
    return {k:rotated(v) for k,v in parts.items()},bottom


def checks(parts,cell,via,bottom):
    c = CELLS[cell]
    checks = {}
    def rec(name,ok,**data):
        checks[name] = dict(ok=bool(ok),**data)
    expected = {"METAL_PEC","SUBSTRATE","AIR","DIODE_SHEET"} | ({"SUBSTRATE_FR4"} if via else set())
    rec("NAMES",set(parts)==expected,names=sorted(parts))
    rec("BREP_VALID",all(BRepCheck_Analyzer(v.wrapped).IsValid() for v in parts.values()))
    solids = parts["METAL_PEC"].Solids()
    rec("ONE_METAL_OBJECT_THREE_SOLIDS",len(solids)==3,count=len(solids))
    solid_parts = [(k,v) for k,v in parts.items() if k not in ("AIR","DIODE_SHEET")]
    overlaps = {}
    for i,(n,a) in enumerate(solid_parts):
        for m,b in solid_parts[i+1:]:
            overlaps[n+"/"+m] = abs(a.intersect(b).Volume())
    for i,a in enumerate(solids):
        for j,b in enumerate(solids[i+1:],i+1):
            overlaps[f"metal{i}/metal{j}"] = abs(a.intersect(b).Volume())
    rec("NO_SOLID_OVERLAP",all(v<1e-8 for v in overlaps.values()),intersection_mm3=overlaps)
    air = parts["AIR"]
    bb = air.BoundingBox()
    rec("AIR_SINGLE_BOX",len(air.Faces())==6 and abs(bb.zmin-bottom)<1e-6 and abs(bb.zmax-LID)<1e-6,faces=len(air.Faces()),z=[bb.zmin,bb.zmax])
    contained = all(v.BoundingBox().zmin>=bb.zmin-1e-6 and v.BoundingBox().zmax<=bb.zmax+1e-6 for v in parts.values())
    rec("AIR_CONTAINS_ALL",contained)
    sheet = parts["DIODE_SHEET"]
    sb = sheet.BoundingBox()
    rec("ONE_SHEET_1_BY_0P35",len(sheet.Faces())==1 and len(sheet.Solids())==0 and abs(sheet.Area()-.35)<1e-7 and abs(sb.xlen-1.)<1e-6 and abs(sb.ylen-.35)<1e-6,area_mm2=sheet.Area(),bbox=[sb.xmin,sb.xmax,sb.ymin,sb.ymax,sb.zmin,sb.zmax])
    edges = []
    for x in (-c["top"],-c["sw"]):
        edge = cq.Edge.makeLine(cq.Vector(x,-.175,T/2),cq.Vector(x,.175,T/2))
        distances = sorted(s.distance(edge) for s in solids)
        edges.append(distances)
    rec("SHEET_EDGE_COINCIDENCE",all(d[0]<1e-7 for d in edges),distances_mm=edges)
    rec("SHEET_NO_INTERIOR_METAL",abs(parts["METAL_PEC"].intersect(sheet).Area())<1e-8)
    # CAD in the rotated calibration frame; verify the top copper away from boundaries.
    bad = total = symmetry_bad = 0
    for ix in range(-52,53,2):
        for iy in range(-52,53,2):
            x,y = ix/10,iy/10
            if min(abs(math.hypot(x,y)-5.435),abs(abs(x)-c["sl"]),abs(abs(y)-c["sw"]),abs(abs(x)-c["il"]),abs(y-c["bot"]),abs(y-c["top"]))<.02:
                continue
            expected_metal = (math.hypot(x,y)<5.435 and not(abs(x)<c["sl"] and abs(y)<c["sw"])) or (abs(x)<c["il"] and c["bot"]<y<c["top"])
            actual = any(s.isInside(cq.Vector(-y,x,T/2),1e-7) for s in solids)
            reflected = any(s.isInside(cq.Vector(-y,-x,T/2),1e-7) for s in solids)
            symmetry_bad+=actual!=reflected
            total+=1;bad+=actual!=expected_metal
    rec("PHYSICAL_TOP_MASK",bad==0,mismatches=bad,points=total)
    rec("TOP_MIRROR_SYMMETRY",symmetry_bad==0,mismatches=symmetry_bad)
    volume = (math.pi*5.435**2-4*c["sl"]*c["sw"]+2*c["il"]*(c["top"]-c["bot"]))*T
    top = box(-P/2,P/2,-P/2,P/2,0,T)
    actual_volume = parts["METAL_PEC"].intersect(top).Volume()
    rec("TOP_VOLUME",abs(actual_volume-volume)<1e-7,actual_mm3=actual_volume,expected_mm3=volume)
    ground_solids=[s for s in solids if abs(s.BoundingBox().zmin+3.210)<1e-6 and abs(s.BoundingBox().zmax+3.175)<1e-6]
    ground_expected=(P*P-(2*math.pi*.300**2 if via else 0))*T
    rec("GROUND_AND_ANTI_PAD_VOLUME",len(ground_solids)==1 and abs(ground_solids[0].Volume()-ground_expected)<1e-7,expected_mm3=ground_expected,actual_mm3=ground_solids[0].Volume() if ground_solids else None,anti_pad_diameter_mm=.6 if via else None)
    slab_area=P*P-(2*math.pi*.150**2 if via else 0)
    rec("SUBSTRATE_BARREL_HOLES",abs(parts["SUBSTRATE"].Volume()-slab_area*H)<1e-7)
    if via:
        rec("FR4_THICKNESS_AND_HOLES",abs(parts["SUBSTRATE_FR4"].Volume()-slab_area*.8)<1e-7)
        pad_stub_overlap=.150*math.sqrt(.300**2-.150**2)+.300**2*math.asin(.150/.300)
        bottom_area=math.pi*.300**2+2.*.3-pad_stub_overlap
        expected_metal=volume+ground_expected+2*(math.pi*.150**2*4.010+bottom_area*T)
        rec("VIA_PAD_STUB_VOLUME",abs(parts["METAL_PEC"].Volume()-expected_metal)<1e-7,expected_mm3=expected_metal,actual_mm3=parts["METAL_PEC"].Volume())
        margin = c["il"]-abs(VIAS[0][0])-.150
        rec("ISLAND_BARREL_FITS",margin>0 and VIAS[0][1]-.150>c["bot"] and VIAS[0][1]+.150<c["top"],minimum_x_clearance_mm=margin)
    return dict(all_pass=all(v["ok"] for v in checks.values()),checks=checks)


def setup(cell,via,bottom):
    c = CELLS[cell]
    hardware = "2 vias (dia 0.3), pads (dia 0.6), anti-pads (dia 0.6), two +y-directed 2.0 x 0.3 mm bottom stubs; FR4 0.8 mm" if via else "no vias, FR4, pads or stubs"
    return f"""# {cell} {'VIAS_FR4_STUBS' if via else 'NO_VIAS'} — HFSS eigenmode (one page)
MODEL_ONLY / NOT_HFSS_VALIDATED. One geometry STEP; OFF/ON differ only in sheet impedance.

**1. Import STEP in mm; assign material pec to METAL_PEC (all three metal solids). STEP imports as vacuum.** SUBSTRATE: RT5880 eps_r=2.2, tanD=0.0009; {'SUBSTRATE_FR4: eps_r=4.4, tanD=0.02.' if via else ''} AIR: vacuum. Copper thickness 0.035 mm. {hardware}.

2. AIR intentionally overlaps all objects. Create project variables $ph_x, $ph_y in degrees. Lattice pairs on its FOUR complete side faces ONLY: x=-5.8/+5.8, phase $ph_x; y=-5.8/+5.8, phase $ph_y. TOP z=150: Perfect E. BOTTOM z={bottom:.3f}: Radiation, **no PEC below the stubs**. Single six-face AIR box retains the calibration cavity fix.

3. Frame is the calibration's +90-degree rotation: slot along y. ACROSS: ($ph_x,$ph_y)=(phi,0); ALONG: (0,phi). p=11.6 mm. ONE Impedance boundary on DIODE_SHEET (1.00 x 0.35 mm): OFF Rs=5.424, Xs=-79.320; ON Rs=1.05, Xs=+7.65 ohm/sq. Current line: (-{c['top']:.3f},0,0.0175) to (-{c['sw']:.3f},0,0.0175). Do not multiply sheet values by W/L again. Record actual boundary values; these are the 5.8-GHz values, replacing the old parallel-RLC model.

4. Eigenmode: min f=4.5 GHz, 8 modes, delta f=0.1%, max 15 passes. Initial mesh: substrate 1 mm; top metal 0.5 mm; diode/gaps/vias/anti-pads 0.1 mm; FR4 0.5 mm. Clean Up Solutions for every state/phase. Start phi=80,85,90,95,100 deg; extend until the SAME top-surface TM branch brackets 5.8 GHz. Track Ez above the patch and decay upward; record all nearby modes/Q. Reject underside/air/via modes. ON mode 1 need not be the wanted branch.

5. Read f (Hz), Q and phase for >=3 nearby points per direction/state; refine phase spacing locally if necessary. Interpolate phi_58 and Q_58. beta=phi_58*pi/(180*0.0116) rad/m; v_g=2*pi*0.0116*df/d(phi radians); alpha_x=2*pi*5.8e9/(2*Q_58*v_g) Np/m. Return results_template.csv and mode-Ez identification. This alpha_x is Q-derived attenuation, NOT the MoM holographic leakage alpha. beta/k0=beta/(2*pi*5.8e9/299792458).

6. Sanity: bare RT5880/ground, 90 deg: historical model 6.291 GHz (150-mm lid); identify the top bound mode. No-via/via comparisons require the same state/phase/branch. Predictions and their limitations are in predictions.json; via+FR4 attenuation has no 3D MoM prediction. Method source: sub6_final_cal.py and results/sub6_final_board/hfss_data.json; archived setup at commit 21fe444, hfss_exports/sub6_island_choke/HFSS_SETUP.md.
"""


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    all_checks = {}
    for cell in CELLS:
        for via in (False,True):
            variant = "vias_fr4_stubs" if via else "novia"
            folder = OUT/(cell+"_"+variant)
            folder.mkdir(exist_ok=True)
            parts,bottom = build(cell,via)
            pre = checks(parts,cell,via,bottom)
            assembly = cq.Assembly()
            for name,shape in parts.items():
                assembly.add(shape,name=name)
            path = folder/(cell+"_"+variant+".step")
            assembly.export(str(path),exportType="STEP")
            raw = path.read_text(encoding="utf-8")
            fixed,stats = exactify_step_text(raw)
            path.write_text(fixed,encoding="utf-8",newline="\n")
            post = checks(read_names(path),cell,via,bottom)
            bad = audit_step_text(fixed)
            vs = verify_snap(raw,fixed)
            exact_ok = len(bad)==0 and vs.get("collapsed_edges",1)==0 and not vs.get("distinct_merge_violations")
            post["checks"]["EXACT_DIMS"] = dict(ok=bool(exact_ok),inexact_literals=len(bad),snap=stats,grid_mm=.001,mid_copper_z_mm=.0175,collapsed_edges=vs.get("collapsed_edges"),distinct_merge_violations=str(vs.get("distinct_merge_violations")))
            post["all_pass"] = post["all_pass"] and exact_ok
            data = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED",cell=cell,variant=variant,physical_dimensions_mm=CELLS[cell],frame="+90deg; ACROSS=x; ALONG=y",air_bottom_mm=bottom,pre=pre,post=post,step_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            (folder/"geometry_checks.json").write_text(json.dumps(data,indent=2),encoding="utf-8",newline="\n")
            all_checks[cell+"_"+variant] = pre["all_pass"] and post["all_pass"]
            (folder/"HFSS_SETUP.md").write_text(setup(cell,via,bottom),encoding="utf-8",newline="\n")
            with (folder/"states.csv").open("w",encoding="utf-8",newline="") as stream:
                writer=csv.writer(stream,lineterminator="\n")
                writer.writerow(["state","geometry","sheet","Rs_ohm_per_sq","Xs_ohm_per_sq","frequency_GHz"])
                writer.writerows([["OFF",path.name,"DIODE_SHEET",5.424,-79.320,5.8],["ON",path.name,"DIODE_SHEET",1.05,7.65,5.8]])
            with (folder/"results_template.csv").open("w",encoding="utf-8",newline="") as stream:
                writer=csv.writer(stream,lineterminator="\n")
                writer.writerow(["cell","variant","state","direction","phase_deg","mode_index","f_GHz","Q","Ez_top_TM_identified","notes"])
                writer.writerows([cell,variant,st,d,ph,mode,"","","",""] for st in ("OFF","ON") for d in ("across_x","along_y") for ph in (80,85,90,95,100) for mode in range(1,9))
            print(cell,variant,"PRE",pre["all_pass"],"POST",post["all_pass"],flush=True)
    (OUT/"geometry_summary.json").write_text(json.dumps(all_checks,indent=2),encoding="utf-8",newline="\n")
    if not all(all_checks.values()):
        raise RuntimeError("HARD_STOP_HFSS_GEOMETRY_GATE_FAILED")


if __name__ == "__main__":
    main()
