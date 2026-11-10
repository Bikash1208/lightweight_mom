"""Isolated BAR50-02V retune; MODEL_ONLY / NOT_HFSS_VALIDATED.

No solver changes. Geometry extends the established island mask with diode x position.
All dimensions recorded in mm; one 1.00 x 0.35 mm lumped sheet. Mesh gap is
8 cells at 96 and 11 at 128, as accepted by the user. Pool uses all logical CPUs.
"""
import argparse
import itertools
import json
import math
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path
from types import SimpleNamespace

# One BLAS thread per process; worker count is not capped.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.linalg import lu_factor, lu_solve
import sub6_slotdisc_lib as S
from src.basis_functions import RooftopBasis
from src.dispersion.leaky_supercell_mom import assemble_z_complex
from src.frequency import free_space_wavenumber
from src.geometry import UnitCell
from src.green_function import LayeredGreenFunction
from src.loads import LumpedLoad, load_matrix
from src.materials import Metal
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.periodic_bc import RectangularLattice
from src.solver import solve

OUT = ROOT / "results/sub6_iter1_island_retune"
SUB = ("rt5880", 2.2, 0.0009, 3.175)
BASE = dict(slot_l=8.7, slot_w=2.9, island_l=4.35, clear=5*11.6/96, diode_x=0.0)
CASES = [("nominal", .11, 4000)] + [
    (f"C{c}_Rp{rp}", c, rp) for c, rp in itertools.product((.105, .12), (3000, 5000))]


class Island:
    kind = "patch_circle_slot_island"
    polarity = "metal"
    rotation = 0.0
    radius = .005435

    def __init__(self, spec, n):
        d = 11.6/n
        self.dx = d
        self.sl = round(spec["slot_l"]/2/d)*d
        self.sw = round(spec["slot_w"]/2/d)*d
        self.il = round(spec["island_l"]/2/d)*d
        self.cb = round(spec["clear"]/d)*d
        self.gap = round(1.0/d)*d  # nearest cells to physical 1.00 mm; unchanged at 96/128
        self.x = round(spec["diode_x"]/d)*d
        self.y = self.sw-self.gap/2
        self.bot = -self.sw+self.cb
        self.top = self.sw-self.gap
        if not (self.top > self.bot and self.sl > self.il > abs(self.x)+.175):
            raise ValueError("Invalid disconnected island / diode geometry")

    @property
    def feature_area(self):
        return math.pi*self.radius**2

    def corners(self):
        r = self.radius
        return np.array([[r,r],[-r,r],[-r,-r],[r,-r]])

    def in_feature(self, x, y):
        x, y = np.asarray(x)*1e3, np.asarray(y)*1e3
        slot = (abs(x) <= self.sl+1e-6) & (abs(y) <= self.sw+1e-6)
        disk = np.hypot(x,y) <= 5.435
        island = (abs(x) <= self.il+1e-6) & (y >= self.bot-1e-6) & (y <= self.top+1e-6)
        sheet = (abs(x-self.x) <= .175) & (abs(y-self.y) <= .5)
        return (disk & ~slot) | (slot & (island | sheet))

    def dimensions(self):
        return dict(slot_l=2*self.sl, slot_w=2*self.sw, island_l=2*self.il,
                    island_w=self.top-self.bot, clearance=self.cb,
                    end_clearance=self.sl-self.il, diode_x=self.x,
                    diode_y=self.y, gap=self.gap, physical_gap=1.0,
                    sheet_width=.35, radius=5.435, period=11.6)


def impedance(f, c, rp, state):
    w = 2*math.pi*f
    return (3+1j*w*.6e-9 if state == "ON" else
            1j*w*.6e-9+1/(1/rp+1j*w*c*1e-12))


def low_rank_currents(basis,geo,green,z0,v,pairs,f):
    """Exact Woodbury identity, using the unchanged load_matrix integrals.

    One nominal-OFF LU; remaining states differ only on sheet-touching
    y-rooftops. No fitted physics, dropped matrix entries or approximations.
    """
    x,y = geo.x*.001,geo.y*.001
    selected = np.flatnonzero((basis.direction == 1) &
        (basis.center_x+basis.dx/2 > x-.000175) &
        (basis.center_x-basis.dx/2 < x+.000175) &
        (basis.center_y+basis.dy > y-.0005) &
        (basis.center_y-basis.dy < y+.0005))
    view = SimpleNamespace(n_unknowns=len(selected),n_x=0,dx=basis.dx,dy=basis.dy,
        center_x=basis.center_x[selected],center_y=basis.center_y[selected])
    unit = LumpedLoad(x,y,length=.001,width=.00035,direction="y",R=1.)
    block = load_matrix(view,[unit],green.omega)
    reference = impedance(f,.11,4000,"OFF")
    z0[np.ix_(selected,selected)] += reference*block
    rhs = np.zeros((basis.n_unknowns,2+len(selected)),complex)
    rhs[:,:2] = v
    rhs[selected,2+np.arange(len(selected))] = 1
    # A.T is an F-contiguous view of C-contiguous A: factor in place, then
    # solve its transpose. This avoids a second dense matrix allocation.
    trans = 1 if z0.flags.c_contiguous else 0
    lu = lu_factor(z0.T if trans else z0,overwrite_a=True,check_finite=False)
    answer = lu_solve(lu,rhs,trans=trans,check_finite=False)
    base, response = answer[:,:2],answer[:,2:]
    out = []
    for name,c,rp,st in pairs:
        difference = (impedance(f,c,rp,st)-reference)*block
        small = np.eye(len(selected))+difference@response[selected]
        correction = np.linalg.solve(small,difference@base[selected])
        cur = base-response@correction
        if not np.isfinite(cur).all():
            raise RuntimeError("HARD_STOP_NONFINITE_WOODBURY_CURRENT")
        out.append(cur)
    return out,len(selected)


def node(task):
    spec, n, f, kt, direction, order, sensitive = task
    t0 = time.time()
    geo = Island(spec, n)
    cell = UnitCell(RectangularLattice(.0116,.0116), geo, S.make_sub(SUB), Metal(), backing="pec")
    green = LayeredGreenFunction(cell.substrate, f, backing="pec")
    basis = RooftopBasis(cell, n, n)
    k0 = free_space_wavenumber(f)
    kx, ky = (kt*k0,0) if direction == "x" else (0,kt*k0)
    fast = os.environ.get("RETUNE_FFT", "0") == "1"
    z0 = (assemble_z_complex(basis,cell.substrate,f,kx,ky,order,order,chunk=int(os.environ.get("RETUNE_FFT_CHUNK","2048"))) if fast
          else assemble_impedance_matrix(basis, green, kx, ky, order))
    v = np.stack([assemble_excitation_vector(basis,e,kx,ky) for e in np.eye(2)],axis=1)
    g00 = np.asarray(green.spectral_dyad(kx,ky), complex)
    result = {}
    pairs = [(name,c,rp,"OFF") for name,c,rp in (CASES if sensitive else CASES[:1])]
    pairs += [("nominal",.11,4000,"ON")]
    woodbury = os.environ.get("RETUNE_WOODBURY","0") == "1"
    currents,rank = low_rank_currents(basis,geo,green,z0,v,pairs,f) if woodbury else (None,None)
    if woodbury:
        del z0
    for index,(name,c,rp,st) in enumerate(pairs):
        zvalue = impedance(f,c,rp,st)
        ld = LumpedLoad(geo.x*1e-3,geo.y*1e-3,length=.001,width=.00035,
                        direction="y",z_of_f=lambda unused: zvalue)
        cur = currents[index] if woodbury else solve(z0+load_matrix(basis,[ld],green.omega),v)
        j = np.stack([basis.current_coefficients(cur[:,i],kx,ky) for i in range(2)],axis=1)
        yt = j @ np.linalg.inv(np.eye(2)+g00@j)
        if not np.isfinite(yt).all():
            raise RuntimeError("HARD_STOP: nonfinite tensor")
        result[name+"|"+st] = [[[float(a.real),float(a.imag)] for a in row] for row in yt]
    return dict(spec=spec,mesh=n,frequency_hz=f,kt=kt,direction=direction,
                floquet_order=order,dimensions_mm=geo.dimensions(),tensors=result,
                unknowns=basis.n_unknowns,seconds=time.time()-t0,
                assembly="existing_fft" if fast else "production_sum",
                solve_method="exact_woodbury" if woodbury else "direct",load_update_size=rank)


def key(spec):
    return "_".join(f"{k}{v:.5f}" for k,v in spec.items())


def summaries(rows):
    grouped = {}
    for row in rows:
        g = grouped.setdefault((key(row["spec"]),row["mesh"],row["frequency_hz"],row["floquet_order"]),[])
        g.append(row)
    out = []
    for (ident,n,f,order), group in grouped.items():
        states = sorted(set(k for row in group for k in row["tensors"]))
        per = {}
        for st in states:
            nodes = {}
            for row in group:
                if st in row["tensors"]:
                    y = np.array(row["tensors"][st])
                    nodes.setdefault(row["kt"],{})[row["direction"]] = {st:y[:,:,0]+1j*y[:,:,1]}
            if any(set(v) != {"x","y"} for v in nodes.values()):
                continue
            r = S.state_summary(nodes,st,SUB,f=f,guess=1.06)
            per[st] = {k:(float(v) if isinstance(v,(float,np.floating)) else v) for k,v in r.items() if not isinstance(v,complex)}
            per[st]["Delta_mS"] = [r["Delta_mS"].real,r["Delta_mS"].imag]
            per[st]["beta_rad_per_m"] = r["beta_over_k0"]*free_space_wavenumber(f)
        if "nominal|OFF" not in per or "nominal|ON" not in per:
            continue
        on = per["nominal|ON"]
        for st,r in per.items():
            if st.endswith("|OFF"):
                r["suppression_dB"] = 10*math.log10(r["alpha_np_per_m"]/on["alpha_np_per_m"])
                r["beta_ON_minus_OFF_over_k0"] = on["beta_over_k0"]-r["beta_over_k0"]
        out.append(dict(id=ident,spec=group[0]["spec"],mesh=n,frequency_hz=f,floquet_order=order,
                        dimensions_mm=group[0]["dimensions_mm"],states=per))
    return out


def candidates():
    specs = [BASE.copy()]
    for knob,values in dict(slot_l=(7.733333333,8.216666667,9.183333333,9.666666667),
                            slot_w=(2.416666667,2.658333333,3.141666667,3.383333333),
                            island_l=(3.383333333,3.866666667,4.833333333,5.316666667),
                            clear=(.241666667,.3625,.845833333,1.0875),diode_x=(-1.0875,-.483333333,.483333333,1.0875)).items():
        for value in values:
            spec = BASE.copy(); spec[knob]=value
            try: Island(spec,96)
            except ValueError: continue
            specs.append(spec)
    return specs


def combined_candidates():
    """Second grid follows the completed single-variable screen, not a surrogate."""
    return [dict(slot_l=l,slot_w=w,island_l=i,clear=c,diode_x=0.)
            for l,w,i,c in itertools.product((8.7,8.941666667,9.183333333),
                (2.658333333,2.9),(3.866666667,4.35),(.241666667,.604166666667))]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage",choices=("baseline","search","search2","refine","frequency"))
    args = parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if args.stage == "baseline":
        specs, meshes, fs, sensitive = [BASE],[96],[5.8e9],False
    elif args.stage == "search":
        specs, meshes, fs, sensitive = candidates(),[96],[5.8e9],False
    elif args.stage == "search2":
        specs, meshes, fs, sensitive = combined_candidates(),[96],[5.8e9],False
    else:
        selected = json.loads((OUT/"selected.json").read_text())
        specs = [r["spec"] for r in selected]
        meshes = [96,128] if args.stage == "refine" else [96]
        fs = [5.8e9] if args.stage == "refine" else [v*1e9 for v in np.arange(4.0,7.61,.2)]
        sensitive = args.stage == "refine"
    tasks = [(s,n,float(f),kt,d,128 if args.stage == "baseline" else 2*n,sensitive)
             for s,n,f,kt,d in itertools.product(specs,meshes,fs,(1.017,1.097),"xy")]
    output_tag = args.stage+os.environ.get("RETUNE_OUTPUT_SUFFIX","")
    path = OUT/(output_tag+"_nodes.jsonl")
    done = []
    if path.exists():
        done = [json.loads(line) for line in path.read_text().splitlines() if line]
    def identity(t): return (key(t[0]),t[1],t[2],t[3],t[4],t[5])
    prior = {(key(r["spec"]),r["mesh"],r["frequency_hz"],r["kt"],r["direction"],r["floquet_order"]) for r in done}
    tasks = [t for t in tasks if identity(t) not in prior]
    print(json.dumps(dict(stage=args.stage,workers=os.cpu_count(),tasks=len(tasks),label="MODEL_ONLY / NOT_HFSS_VALIDATED")),flush=True)
    with path.open("a",encoding="utf-8") as stream, Pool() as pool:
        for row in pool.imap_unordered(node,tasks):
            done.append(row)
            stream.write(json.dumps(row)+"\n"); stream.flush()
            print(f"done {len(done)}: mesh {row['mesh']} f {row['frequency_hz']/1e9:.2f} kt {row['kt']} {row['direction']} {row['seconds']:.1f}s",flush=True)
    result = summaries(done)
    (OUT/(output_tag+".json")).write_text(json.dumps(dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED",rows=result),indent=2),encoding="utf-8")
    for r in result:
        a,b = r["states"]["nominal|OFF"],r["states"]["nominal|ON"]
        print(r["id"],r["mesh"],f"alpha {a['alpha_np_per_m']:.5g}/{b['alpha_np_per_m']:.5g} supp {a['suppression_dB']:.3f}",flush=True)


if __name__ == "__main__":
    main()
