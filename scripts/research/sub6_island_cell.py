"""SUB6 island cell scan (MODEL_ONLY / NOT_HFSS_VALIDATED).  jobs/sub6-island-cell.md.

Cell: closed slotted disc (anchor r 5.435, p 11.6, RT5880 3.175, slot L x W along x, as SlotDiscNotch without a notch) + a rectangular metal ISLAND centred in the slot
+ ONE BAR50 lumped sheet (1.00 mm across y x 0.35 mm along x) bridging the 1.00 mm clearance between the island's +y long edge and the slot's +y edge, at x = 0.
Dimension conventions (this script, flagged):
  * clearance 'c' (0.3 / 0.6 mm) = the clearance between the island's -y long edge and the slot's -y edge ("other side"); the island length Li = (Li/L) L
    sets the clearance to the two short slot ends ((L - Li)/2); only the -y side takes the swept clearance.
  * mesh 64 on p = 11.6 mm: ONE MESH CELL = 0.18125 mm.  The mask is sampled at cell centres, so every edge is snapped to a multiple of the mesh cell:
    slot half length / half width = floor(L/2/dx + 1/2) cells; island half length = round(Li/2/dx) cells; clearance c -> round(c/dx) cells (0.3 -> 2 cells = 0.3625 mm,
    0.6 -> 3 cells = 0.544 mm); the 1.00 mm +y gap -> 5 cells = 0.906 mm (the same convention as the 1.00 mm rim notch of the cell scan, where the sheet length is 1.00 mm
    over a 5-cell = 0.906 mm gap).  Realised dimensions are written to every row.
  * ISLAND_ABSENT (reference) = the closed slotted disc (NOCUT), one per (W, L/2r).
States OFF_BAR50 ((0.15 pF || 5 k) + 0.6 nH), ON_BAR50 (3.0 ohm + 0.6 nH) from sub6_slotdisc_lib (unchanged).  Self-consistent beta, Delta = Y_par - Y_perp, A_perp, alpha
(leakage_alpha_first_order_tensor at 25 deg), rho = Delta_ON / Delta_OFF, leakage class (threshold 1 Np/m), resonance proxy (sign of Im Y_xx/Y_yy at kt = 0, 5.0 vs 6.6 GHz).
Gate OK: alpha_ON or alpha_OFF >= 1 Np/m AND |rho - 1| >= 1.
Scheduling: all (row, direction) node assemblies of one round go to one pool (task-level parallel); nodes added adaptively as in sub6_slotdisc_grid (max 4).
    python sub6_island_cell.py list | run [--workers N] | report
"""
import argparse
import json
import math
import os
import sys
import time
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_slotdisc_grid as G
import sub6_slotdisc_lib as S
from src.basis_functions import MeshSettings, RooftopBasis
from src.frequency import free_space_wavenumber
from src.geometry import UnitCell
from src.green_function import LayeredGreenFunction
from src.loads import LumpedLoad, load_matrix
from src.materials import Metal
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.periodic_bc import RectangularLattice
from src.solver import solve

ROOT = G.ROOT
OUT = ROOT / "results" / os.environ.get("ISL_OUT", "sub6_island_cell")  # env overrides: mesh-check runs (defaults = the original mesh-64 scan)
MM = 1e-3
SUB = ("rt5880", 2.2, 0.0009, 3.175)
P, R_ANCHOR = 11.6, 5.435
N_MESH = int(os.environ.get("ISL_MESH", "64"))
FORDER = 2 * N_MESH
DX = P / N_MESH  # mm
WS = (2.9, 3.6, 4.4)
RATIOS = (0.7, 0.8, 0.9)
LI_OVER_L = (0.3, 0.5, 0.7)
CLEAR = (0.3, 0.6)
if os.environ.get("ISL_GRID") == "refine":  # job sub6-island-refine: grid around ISL_W2.9_q0.8_Li0.5_c0.6 (diode gap stays 1.00 mm = GAP_CELLS)
    WS, RATIOS, LI_OVER_L, CLEAR = (2.6, 2.9, 3.2), (0.75, 0.8, 0.85), (0.4, 0.5, 0.6), (0.5, 0.6, 0.7)
GAP_CELLS = int(os.environ.get("ISL_GAP_CELLS", "5"))  # 5 cells = 0.906 mm at mesh 64; mesh 96: 8 cells = 0.967 mm
CORNER_MIN, RIM_MIN = 0.2, 0.35
ISL_STATES = ["OFF_BAR50", "ON_BAR50"]


@dataclass(frozen=True)
class SlotDiscIsland:
    """Duck-types geometry.Geometry (like SlotDiscNotch).  Integer mesh-cell counts (see the module docstring); `island` False = closed slotted disc."""
    radius: float  # m
    n_l: int  # slot half length in cells
    n_w: int  # slot half width in cells
    n_i: int  # island half length in cells
    cb: int  # clearance island -y edge to slot -y edge, cells
    island: bool = True
    kind = "patch_circle_slot_island"
    polarity = "metal"
    rotation = 0.0

    @property
    def dxm(self):
        return DX * MM

    @property
    def y_top_island(self):
        return (self.n_w - GAP_CELLS) * self.dxm

    @property
    def y_bot_island(self):
        return (-self.n_w + self.cb) * self.dxm

    @property
    def diode_y(self):
        return 0.5 * (self.y_top_island + self.n_w * self.dxm)

    @property
    def feature_area(self):
        return math.pi * self.radius**2

    def corners(self):
        r = self.radius
        return np.array([[r, r], [-r, r], [-r, -r], [r, -r]])

    def in_feature(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        d = self.dxm
        inside = np.hypot(x, y) <= self.radius
        slot = (np.abs(x) <= self.n_l * d + 1e-9) & (np.abs(y) <= self.n_w * d + 1e-9)
        inside &= ~slot
        if self.island:
            isl = (np.abs(x) <= self.n_i * d + 1e-9) & (y >= self.y_bot_island - 1e-9) & (y <= self.y_top_island + 1e-9)
            sheet = (np.abs(x) <= S.DIODE_ALONG / 2) & (np.abs(y - self.diode_y) <= S.G_C / 2)
            inside = inside | (slot & (isl | sheet))
        return inside


def cells_of(L, W, li_over_l, c):
    """Integer mesh-cell geometry (and the realised mm dimensions) for a nominal row."""
    n_l = int(math.floor(L / 2 / DX + 0.5))
    n_w = int(math.floor(W / 2 / DX + 0.5))
    n_i = max(1, int(round(li_over_l * L / 2 / DX)))
    cb = max(1, int(round(c / DX)))
    return n_l, n_w, n_i, cb


def realised(n_l, n_w, n_i, cb):
    return dict(slot_L_mm=2 * n_l * DX, slot_W_mm=2 * n_w * DX, island_Li_mm=2 * n_i * DX, island_Wi_mm=(2 * n_w - GAP_CELLS - cb) * DX,
                clear_minus_y_mm=cb * DX, gap_plus_y_mm=GAP_CELLS * DX, clear_short_ends_mm=(n_l - n_i) * DX, mesh_cell_mm=DX)


def rows():
    out, skipped = [], []
    r = R_ANCHOR
    for W in WS:
        for q in RATIOS:
            L = 2 * r * q
            if r - L / 2 < RIM_MIN - 1e-9 or r - math.hypot(L / 2, W / 2) < CORNER_MIN:
                skipped.append(f"W{W:g}_q{q:g}")
                continue
            for li in LI_OVER_L:
                for c in CLEAR:
                    n_l, n_w, n_i, cb = cells_of(L, W, li, c)
                    out.append(dict(id=f"ISL_W{W:g}_q{q:g}_Li{li:g}_c{c:g}", W=W, ratio=q, L=L, r=r, li_over_l=li, clear=c, cells=[n_l, n_w, n_i, cb], real=realised(n_l, n_w, n_i, cb)))
    refs = sorted({(W, q) for W in WS for q in RATIOS if f"W{W:g}_q{q:g}" not in skipped})
    return out, refs, skipped


def make_island_cell(r_mm, cells, state, island=True):
    n_l, n_w, n_i, cb = cells
    geo = SlotDiscIsland(r_mm * MM, n_l, n_w, n_i, cb, island=island)
    loads = ()
    if island and state != "NOCUT":
        ld = S.STATES[state][1]
        loads = (LumpedLoad(0.0, geo.diode_y, length=S.G_C, width=S.DIODE_ALONG, direction="y", **ld),)
    return UnitCell(RectangularLattice(P * MM, P * MM), geo, S.make_sub(SUB), Metal(), backing="pec", loads=loads)


def node_task(a):
    """(kind 'ISL'|'REF', r, cells, direction, kt, freq) -> {state: Y}.  One Z assembly; ISL: OFF/ON share it."""
    kind, r, cells, direction, kt, f = a
    t0 = time.time()
    island = kind == "ISL"
    states = ISL_STATES if island else ["NOCUT"]
    cell0 = make_island_cell(r, cells, states[0], island)
    k0 = free_space_wavenumber(f)
    kx, ky = (kt * k0, 0.0) if direction == "x" else (0.0, kt * k0)
    green = LayeredGreenFunction(cell0.substrate, f, backing="pec")
    basis = RooftopBasis(cell0, N_MESH, N_MESH)
    z0 = assemble_impedance_matrix(basis, green, kx, ky, FORDER)
    v = np.stack([assemble_excitation_vector(basis, e, kx, ky) for e in np.eye(2)], axis=1)
    g00 = np.asarray(green.spectral_dyad(kx, ky), dtype=complex)
    out = {}
    for st in states:
        c = make_island_cell(r, cells, st, island)
        z = z0 + load_matrix(basis, list(c.loads), green.omega) if c.loads else z0
        cur = solve(z, v)
        j = np.stack([basis.current_coefficients(cur[:, k], kx, ky) for k in range(2)], axis=1)
        out[st] = j @ np.linalg.inv(np.eye(2) + g00 @ j)
    return dict(kind=kind, direction=direction, kt=kt, f=f, Y=out, seconds=time.time() - t0, unknowns=int(basis.n_unknowns))


def ref_cells(W, q):
    L = 2 * R_ANCHOR * q
    n_l, n_w, _, _ = cells_of(L, W, 0.5, 0.3)
    return [n_l, n_w, 1, 1]


def ckpt_path(key):
    return OUT / "ckpt" / f"{key}.json"


def next_nodes_fix(nodes, states, tol=0.004, cap=9):
    """New nodes at each state's own self-consistent beta when it is farther than tol from every existing node (max `cap` nodes)."""
    out = []
    for st in states:
        try:
            r = S.state_summary(nodes, st, SUB, guess=1.06)
        except Exception:
            continue
        b = r["beta_over_k0"]
        if r["nearest_node_dist"] > tol and not any(abs(b - k) < 3e-4 for k in list(nodes) + out) and len(nodes) + len(out) < cap:
            out.append(float(b))
    return out


def run(workers):
    OUT.mkdir(parents=True, exist_ok=True)
    rs, refs, skipped = rows()
    if os.environ.get("ISL_ROWS"):  # mesh-check subset: only these row ids, no reference geometries
        keep = set(os.environ["ISL_ROWS"].split(","))
        rs, refs = [r for r in rs if r["id"] in keep], []
    jobs = {r["id"]: dict(kind="ISL", r=R_ANCHOR, cells=r["cells"], states=ISL_STATES) for r in rs}
    for W, q in refs:
        jobs[f"REF_W{W:g}_q{q:g}"] = dict(kind="REF", r=R_ANCHOR, cells=ref_cells(W, q), states=["NOCUT"])
    nodes = {k: G.load_nodes(ckpt_path(k)) for k in jobs}
    t0 = time.time()
    with Pool(workers) as pool:
        for rnd in range(8 if os.environ.get("ISL_NODEFIX") else 5):
            todo = []
            for k, j in jobs.items():
                if os.environ.get("ISL_NODEFIX"):  # one node at EACH state's own beta (the median rule of S.next_node lands between OFF and ON betas -> extrapolation)
                    kts = next_nodes_fix(nodes[k], j["states"]) if nodes[k] else [float(x) for x in os.environ.get("ISL_KT", "1.06").split(",")]
                else:
                    kt = 1.06 if not nodes[k] else S.next_node(nodes[k], j["states"], SUB)
                    kts = [] if kt is None else [kt]
                for kt in kts:
                    for d in "xy":
                        todo.append((k, (j["kind"], j["r"], j["cells"], d, kt, S.F)))
            if not todo:
                print("nodes converged", flush=True)
                break
            print(f"round {rnd}: {len(todo)} assemblies, {workers} workers, t={time.time() - t0:.0f}s", flush=True)
            for (k, a), res in zip(todo, pool.map(node_task, [a for _, a in todo], chunksize=1)):
                nodes[k].setdefault(res["kt"], {}).setdefault(res["direction"], {}).update(res["Y"])
            for k in {k for k, _ in todo}:
                G.save_nodes(ckpt_path(k), nodes[k])
        # resonance proxy: normal-incidence Y_xx, Y_yy at 5.0 / 6.6 GHz (x-direction assembly, kt = 0 gives both diagonals)
        res_tasks = [(k, (j["kind"], j["r"], j["cells"], "x", 0.0, fr)) for k, j in jobs.items() for fr in (5.0e9, 6.6e9)]
        print(f"resonance proxy: {len(res_tasks)} assemblies t={time.time() - t0:.0f}s", flush=True)
        resv = {}
        for (k, a), res in zip(res_tasks, pool.map(node_task, [a for _, a in res_tasks], chunksize=1)):
            resv.setdefault(k, {})[a[5]] = res["Y"]
    out_rows, out_refs = [], {}
    for k, j in jobs.items():
        states = {st: S.state_summary(nodes[k], st, SUB, guess=1.06) for st in j["states"]}
        flip = []
        for st in j["states"]:
            for i in (0, 1):
                if np.sign(resv[k][5.0e9][st][i, i].imag) != np.sign(resv[k][6.6e9][st][i, i].imag):
                    flip.append(f"{st}:{'xy'[i]}")
        if j["kind"] == "REF":
            out_refs[k] = dict(states=G.jc(states), nodes=sorted(nodes[k]), resonance_flip=flip)
            continue
        row = [r for r in rs if r["id"] == k][0]
        ron, roff = states["ON_BAR50"], states["OFF_BAR50"]
        rho = ron["Delta_mS"] / roff["Delta_mS"]
        a_on, a_off = ron["alpha_np_per_m"], roff["alpha_np_per_m"]
        cls = G.leakage_class(a_on, a_off)
        gate = "OK" if (a_on >= 1.0 or a_off >= 1.0) and abs(rho - 1.0) >= 1.0 else "NOT_OK"
        out_rows.append(dict(row=row, states=G.jc(states), nodes=sorted(nodes[k]), resonance_flip=flip, resonance_flag=bool(flip), rho=G.jc(rho), abs_rho_minus_1=float(abs(rho - 1.0)),
                             alpha_on=a_on, alpha_off=a_off, leakage=cls, gate=gate))
    json.dump(dict(rows=out_rows, refs=out_refs, skipped=skipped), open(OUT / "island_rows.json", "w"), indent=1, default=float)
    prov = G.provenance(dict(stage="island", workers=workers, seconds=time.time() - t0, rows=len(out_rows), refs=len(out_refs), mesh_cell_mm=DX))
    (OUT / "provenance_run.json").write_text(json.dumps(prov, indent=1, default=str))
    print(f"done {time.time() - t0:.0f}s: rows {len(out_rows)} OK={sum(r['gate'] == 'OK' for r in out_rows)}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["list", "run", "selftest"])
    ap.add_argument("--workers", type=int, default=64)
    a = ap.parse_args()
    if a.stage == "list":
        rs, refs, sk = rows()
        print(len(rs), "island rows,", len(refs), "reference geometries, skipped (corner/rim):", sk)
        for r in rs[:6]:
            print(r["id"], r["cells"], {k: round(v, 3) for k, v in r["real"].items()})
    elif a.stage == "selftest":
        rs, refs, sk = rows()
        r = rs[0]
        c = make_island_cell(r["r"], r["cells"], "ON_BAR50")
        b = RooftopBasis(c, N_MESH, N_MESH)
        m = b.metal_mask
        print(r["id"], "unknowns", b.n_unknowns, "metal cells", int(m.sum()))
        j = N_MESH // 2
        # print the central mask rows around the slot (x horizontal)
        n_l, n_w, n_i, cb = r["cells"]
        for jj in range(j + n_w + 1, j - n_w - 2, -1):
            print("".join("#" if m[ii, jj] else "." for ii in range(j - n_l - 2, j + n_l + 2)))
    else:
        run(a.workers)
