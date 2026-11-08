"""SUB6 island cell: mesh 96 vs 128 convergence check (MODEL_ONLY).  Built on sub6_island_cell.py (compute/sub6-island-cell).

Rows ISL_W2.9_q0.8_Li0.5_c0.6 (board row) and ISL_W2.9_q0.85_Li0.5_c0.5 (refine row), states OFF_BAR50 / ON_BAR50.
Mesh cell dx = p/N.  Same convention as the base script, but with the integer cell counts re-derived per mesh and the 1.00 mm diode gap = round(1.00/dx) cells:
slot half length/width = floor(L/2/dx + 1/2), island half length = round(Li/2/dx), clearance(-y) = round(c/dx), +y gap = round(1.00/dx).
N is even, so the x = 0 mirror plane is a mesh-cell edge: slot, island and diode sheet are mirror-symmetric in x (2 n cells).  Floquet order F = 2N (as before).
    python sub6_island_mesh128.py list | run [--workers N] | collect
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

import sub6_island_cell as IC
import sub6_slotdisc_grid as G
import sub6_slotdisc_lib as S
from src.basis_functions import RooftopBasis
from src.frequency import free_space_wavenumber
from src.geometry import UnitCell
from src.green_function import LayeredGreenFunction
from src.loads import LumpedLoad, load_matrix
from src.materials import Metal
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.periodic_bc import RectangularLattice
from src.solver import solve

OUT = IC.ROOT / "results" / "sub6_island_mesh128"
MM, P, R, SUB = 1e-3, IC.P, IC.R_ANCHOR, IC.SUB
MESHES = (96, 128)
ROWS = [dict(id="ISL_W2.9_q0.8_Li0.5_c0.6", W=2.9, q=0.8, li=0.5, c=0.6), dict(id="ISL_W2.9_q0.85_Li0.5_c0.5", W=2.9, q=0.85, li=0.5, c=0.5)]
STATES = IC.ISL_STATES
SEEDS = (1.02, 1.08, 1.16)
RES_F = (5.0e9, 5.8e9, 6.6e9)


@dataclass(frozen=True)
class Isl(IC.SlotDiscIsland):
    dx_mm: float = 0.0
    gap: int = 5

    @property
    def dxm(self):
        return self.dx_mm * MM

    @property
    def y_top_island(self):
        return (self.n_w - self.gap) * self.dxm


def cells(N, row):
    dx = P / N
    L = 2 * R * row["q"]
    n_l = int(math.floor(L / 2 / dx + 0.5))
    n_w = int(math.floor(row["W"] / 2 / dx + 0.5))
    n_i = max(1, int(round(row["li"] * L / 2 / dx)))
    cb = max(1, int(round(row["c"] / dx)))
    gap = max(1, int(round(1.0 / dx)))
    return n_l, n_w, n_i, cb, gap


def real(N, row):
    dx = P / N
    n_l, n_w, n_i, cb, gap = cells(N, row)
    return dict(mesh=N, dx_mm=dx, slot_L_mm=2 * n_l * dx, slot_W_mm=2 * n_w * dx, island_Li_mm=2 * n_i * dx, island_Wi_mm=(2 * n_w - gap - cb) * dx,
                clear_minus_y_mm=cb * dx, gap_plus_y_mm=gap * dx, clear_short_ends_mm=(n_l - n_i) * dx, cells=[n_l, n_w, n_i, cb, gap],
                nominal=dict(slot_L_mm=2 * R * row["q"], slot_W_mm=row["W"], island_Li_mm=row["li"] * 2 * R * row["q"], clear_mm=row["c"], gap_mm=1.0))


def make_cell(N, c, state):
    n_l, n_w, n_i, cb, gap = c
    geo = Isl(R * MM, n_l, n_w, n_i, cb, True, dx_mm=P / N, gap=gap)
    ld = S.STATES[state][1]
    return UnitCell(RectangularLattice(P * MM, P * MM), geo, S.make_sub(SUB), Metal(), backing="pec",
                    loads=(LumpedLoad(0.0, geo.diode_y, length=S.G_C, width=S.DIODE_ALONG, direction="y", **ld),))


def task(a):
    N, c, direction, kt, f = a
    t0 = time.time()
    cell0 = make_cell(N, c, STATES[0])
    k0 = free_space_wavenumber(f)
    kx, ky = (kt * k0, 0.0) if direction == "x" else (0.0, kt * k0)
    green = LayeredGreenFunction(cell0.substrate, f, backing="pec")
    basis = RooftopBasis(cell0, N, N)
    z0 = assemble_impedance_matrix(basis, green, kx, ky, 2 * N)
    v = np.stack([assemble_excitation_vector(basis, e, kx, ky) for e in np.eye(2)], axis=1)
    g00 = np.asarray(green.spectral_dyad(kx, ky), dtype=complex)
    out = {}
    for st in STATES:
        cl = make_cell(N, c, st)
        cur = solve(z0 + load_matrix(basis, list(cl.loads), green.omega), v)
        j = np.stack([basis.current_coefficients(cur[:, k], kx, ky) for k in range(2)], axis=1)
        out[st] = j @ np.linalg.inv(np.eye(2) + g00 @ j)
    return dict(direction=direction, kt=kt, f=f, Y=out, seconds=time.time() - t0, unknowns=int(basis.n_unknowns))


def ck(N, rid):
    return OUT / "ckpt" / f"m{N}_{rid}.json"


def run(workers):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "ckpt").mkdir(exist_ok=True)
    jobs = {(N, r["id"]): cells(N, r) for N in MESHES for r in ROWS}
    nodes = {k: G.load_nodes(ck(*k)) for k in jobs}
    resf = OUT / "resonance.json"
    resv = json.loads(resf.read_text()) if resf.exists() else {}
    t0 = time.time()

    def log(m):
        print(f"[{time.time() - t0:6.0f}s] {m}", flush=True)

    with Pool(workers) as pool:
        for rnd in range(3):
            todo = []
            for k, c in jobs.items():
                if not nodes[k]:
                    kts = SEEDS
                else:
                    kt = S.next_node(nodes[k], STATES, SUB)
                    kts = (kt,) if kt is not None else ()
                for kt in kts:
                    for d in "xy":
                        todo.append((k, (k[0], c, d, kt, S.F)))
            res_t = []
            if rnd == 0:
                res_t = [(k, (k[0], c, "x", 0.0, f)) for k, c in jobs.items() for f in RES_F if f"{k[0]}|{k[1]}|{f:g}" not in resv]
            allt = sorted(todo + res_t, key=lambda t: -t[0][0])  # mesh 128 first
            if not allt:
                log("nodes converged")
                break
            log(f"round {rnd}: {len(todo)} node assemblies + {len(res_t)} resonance, {workers} workers")
            for i, (kk, res) in enumerate(zip([t[0] for t in allt], pool.imap(task, [t[1] for t in allt], chunksize=1))):
                a = [t for t in allt if t[0] == kk][0]
                if res["kt"] == 0.0 if "kt" in res else False:
                    pass
                if res["kt"] == 0.0:
                    resv[f"{kk[0]}|{kk[1]}|{res['f']:g}"] = {st: np.asarray(y).tolist() for st, y in res["Y"].items()}
                    resv = {k2: ({s: [[complex(v) for v in row] for row in m] for s, m in v2.items()} if False else v2) for k2, v2 in resv.items()}
                else:
                    nodes[kk].setdefault(res["kt"], {}).setdefault(res["direction"], {}).update(res["Y"])
                log(f"{i + 1}/{len(allt)} mesh {kk[0]} {kk[1]} {res['direction']} kt={res['kt']:.4f} f={res['f'] / 1e9:.1f} {res['seconds']:.0f}s unknowns {res['unknowns']}")
                if (i + 1) % 4 == 0 or i + 1 == len(allt):
                    for k2 in jobs:
                        G.save_nodes(ck(*k2), nodes[k2])
                    json.dump({k2: {s: [[[z.real, z.imag] if isinstance(z, complex) else z for z in row] for row in np.asarray(m)] for s, m in v2.items()} if not isinstance(next(iter(v2.values())), list) else v2 for k2, v2 in resv.items()} if False else {k2: {s: np.asarray(m).tolist() if not isinstance(m, list) else m for s, m in v2.items()} for k2, v2 in resv.items()}, open(resf, "w"), default=lambda o: [o.real, o.imag] if isinstance(o, complex) else float(o))
    json.dump(dict(workers=workers, seconds=time.time() - t0, meshes=MESHES, realised={f"{N}|{r['id']}": real(N, r) for N in MESHES for r in ROWS}), open(OUT / "run_meta.json", "w"), indent=1)
    print("RUN_DONE", flush=True)


def collect():
    prov_rows = []
    out = {}
    for r in ROWS:
        out[r["id"]] = dict(realised={N: real(N, r) for N in MESHES}, mesh={})
        for N in MESHES:
            nd = G.load_nodes(ck(N, r["id"]))
            out[r["id"]]["mesh"][N] = {st: S.state_summary(nd, st, SUB, guess=1.06) for st in STATES}
            out[r["id"]]["mesh"][N]["nodes"] = sorted(nd)
    json.dump(out, open(OUT / "mesh128_check.json", "w"), indent=1, default=lambda o: G.jc(o) if True else float(o))
    for rid, d in out.items():
        print(rid)
        for N in MESHES:
            rl = d["realised"][N]
            print(f"  mesh {N}: dx {rl['dx_mm']:.4f}  slot {rl['slot_L_mm']:.3f}x{rl['slot_W_mm']:.3f}  Li {rl['island_Li_mm']:.3f}  clear {rl['clear_minus_y_mm']:.3f}  gap {rl['gap_plus_y_mm']:.3f}  ends {rl['clear_short_ends_mm']:.3f}  nodes {d['mesh'][N]['nodes']}")
        for st in STATES:
            a, b = d["mesh"][96][st], d["mesh"][128][st]
            for key, nm in (("alpha_np_per_m", "alpha"), ("beta_over_k0", "beta/k0")):
                print(f"    {st:9s} {nm:8s} {a[key]:.4f} -> {b[key]:.4f}  {100 * (b[key] / a[key] - 1):+.1f}%")
            da, db = abs(complex(a["Delta_mS"])), abs(complex(b["Delta_mS"]))
            print(f"    {st:9s} |Delta|  {da:.3f} -> {db:.3f}  {100 * (db / da - 1):+.1f}%   extrapolated {a['extrapolated']}/{b['extrapolated']}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["list", "run", "collect"])
    ap.add_argument("--workers", type=int, default=60)
    a = ap.parse_args()
    if a.stage == "list":
        for r in ROWS:
            for N in MESHES:
                print(r["id"], json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in real(N, r).items() if k != "nominal"}))
    elif a.stage == "run":
        run(a.workers)
    else:
        collect()
