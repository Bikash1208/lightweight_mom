"""SUB6 edge-diode cell (MODEL_ONLY / IDEAL_BIAS / NOT_HFSS_VALIDATED).  Job sub6-edge-diode-cell.

Closed slotted disc (rim intact) on a p = 11.6 mm square lattice, RT5880 3.175 mm on PEC, 5.8 GHz.  The disc radius is set by the inter-cell gap:
r = (p - 1.00) / 2 = 5.30 mm.  The BAR50 sheet (1.00 along x x 0.35 across, direction "x") bridges the gap to the +x neighbour, centred on the cell
edge.  In the periodic model EVERY cell has a diode to its +x edge neighbour: a line of switched cells (a row of cells is connected through the
sheets), not an isolated diode.

Implementation: the unit cell is translated by p/2 in x so the sheet sits at the cell centre (0, 0) (an interior LumpedLoad; the load integrals do not wrap
across the periodic boundary) and the disc centre is at the cell edge (x = +-p/2, wrapped mask).  A translation of the origin does not change the sheet
tensor Y(kt) (checked in `selftest` against the unshifted slotted disc).
Slot angle theta relative to the diode axis (+x) = 0 / 45 / 90 deg: theta = 0 is the diode ON the slot axis, 90 the diode at the slot end-wall side
(the diode at the +y edge is the same lattice problem, rotated).  A diode at a 45 deg cell CORNER has no 1.00 mm neighbour gap (the diagonal disc-to-disc
gap is p*sqrt(2) - 2r = 5.8 mm), so 45 deg is realised as the slot at 45 deg to the diode axis.
Y is computed for kt along the slot axis u (Y_par) and across it v (Y_perp) and rotated into the slot frame (Y_uu, Y_vv, A, and the off-diagonal Y_uv, which is
non-zero at 45 deg because the mirror about the diode axis is broken).
Mesh 64 (h = 0.181 mm), F = 2n = 128.  Gap resolution: 6 cell columns (1.087 mm) x 2 cell rows (0.363 mm) for the 1.00 x 0.35 sheet.
States: OFF_BAR50 = (0.15 pF || 5 kohm) + 0.6 nH, ON_BAR50 = 3.0 ohm + 0.6 nH, NO_DIODE = empty gap (no sheet).

    python sub6_edge_diode_cell.py selftest | run | res | collect
"""
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
from src._validate import finite_real
from src.basis_functions import RooftopBasis
from src.frequency import free_space_wavenumber
from src.geometry import UnitCell
from src.green_function import LayeredGreenFunction
from src.loads import LumpedLoad, load_matrix
from src.materials import Metal
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.periodic_bc import RectangularLattice
from src.solver import solve

MM = 1e-3
OUT = G.ROOT / "results/sub6_edge_diode"
N, FO = 64, 128
P = 11.6
GAP = 1.00
R = (P - GAP) / 2
SUB = ("rt5880", 2.2, 0.0009, 3.175)
SHEET_L, SHEET_W = 1.00 * MM, 0.35 * MM
THETAS = (0, 45, 90)
QS = (0.7, 0.8, 0.9)
WS = (2.4, 2.9, 3.6)
CLASS_STATES = {"NODIODE": ["NO_DIODE"], "BRIDGE": ["OFF_BAR50", "ON_BAR50"]}
ALL_STATES = ["NO_DIODE", "OFF_BAR50", "ON_BAR50"]
LOADS = {"OFF_BAR50": dict(z_of_f=S.z_off_bar50), "ON_BAR50": dict(R=3.0, L=0.6e-9)}
GUESS = 1.06
WORKERS = 61  # decision 2026-11-08: RAM is not a limit, always the full pool


@dataclass(frozen=True)
class EdgeDisc:
    radius: float
    slot_length: float
    slot_width: float
    slot_angle: float  # rad from +x (the diode axis)
    period: float
    bridge: bool = False  # True: the sheet region (1.00 x 0.35 mm at the cell centre) is metal
    kind = "patch_circle_rect_slot_edge"
    polarity = "metal"
    rotation = 0.0

    def __post_init__(self):
        for n in ("radius", "slot_length", "slot_width", "slot_angle", "period"):
            object.__setattr__(self, n, finite_real(n, getattr(self, n)))

    @property
    def feature_area(self):
        return math.pi * self.radius**2 - self.slot_length * self.slot_width + (SHEET_L * SHEET_W if self.bridge else 0.0)

    def corners(self):
        h = self.period / 2
        return np.array([[h, h], [-h, h], [-h, -h], [h, -h]])

    def in_feature(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        dx = np.mod(x, self.period) - self.period / 2  # disc centre at the cell edge x = +-p/2 (periodic mask)
        c, s = math.cos(self.slot_angle), math.sin(self.slot_angle)
        u, v = dx * c + y * s, -dx * s + y * c
        inside = np.hypot(dx, y) <= self.radius
        inside &= ~((np.abs(u) <= self.slot_length / 2) & (np.abs(v) <= self.slot_width / 2))
        if self.bridge:
            inside |= (np.abs(x) <= SHEET_L / 2) & (np.abs(y) <= SHEET_W / 2)
        return inside


def make_cell(L_mm, W_mm, theta_deg, state):
    bridge = state != "NO_DIODE"
    geo = EdgeDisc(R * MM, L_mm * MM, W_mm * MM, math.radians(theta_deg), P * MM, bridge=bridge)
    loads = (LumpedLoad(0.0, 0.0, length=SHEET_L, width=SHEET_W, direction="x", **LOADS[state]),) if bridge else ()
    return UnitCell(RectangularLattice(P * MM, P * MM), geo, S.make_sub(SUB), Metal(), backing="pec", loads=loads)


def rows():
    out = []
    for th in THETAS:
        for q in QS:
            L = 2 * R * q
            for w in WS:
                out.append(dict(id=f"EDGE_th{th}_q{q:g}_W{w:g}", theta=th, q=q, L=L, W=w, r=R, sub=SUB))
    return out


def slot_frame(Y, th):
    c, s = math.cos(math.radians(th)), math.sin(math.radians(th))
    Rm = np.array([[c, -s], [s, c]])
    return Rm.T @ Y @ Rm


def node_task(a):
    """(row, class, direction 'x'|'y' = slot axis u | across v, kt, n, fo, f) -> {state: 2x2 Y in the slot frame}."""
    row, cls, d, kt, n, fo, f = a
    th = row["theta"]
    cell = make_cell(row["L"], row["W"], th, CLASS_STATES[cls][0])
    k0 = free_space_wavenumber(f)
    ang = math.radians(th) + (0.0 if d == "x" else math.pi / 2)
    kx, ky = kt * k0 * math.cos(ang), kt * k0 * math.sin(ang)
    green = LayeredGreenFunction(cell.substrate, f, backing="pec")
    basis = RooftopBasis(cell, n, n)
    t0 = time.time()
    z0 = assemble_impedance_matrix(basis, green, kx, ky, fo)
    v = np.stack([assemble_excitation_vector(basis, e, kx, ky) for e in np.eye(2)], axis=1)
    g00 = np.asarray(green.spectral_dyad(kx, ky), dtype=complex)
    out = {}
    for st in CLASS_STATES[cls]:
        z = z0
        if st != "NO_DIODE":
            ld = LumpedLoad(0.0, 0.0, length=SHEET_L, width=SHEET_W, direction="x", **LOADS[st])
            z = z0 + load_matrix(basis, [ld], green.omega)
        cur = solve(z, v)
        j = np.stack([basis.current_coefficients(cur[:, k], kx, ky) for k in range(2)], axis=1)
        out[st] = slot_frame(j @ np.linalg.inv(np.eye(2) + g00 @ j), th)
    return dict(row=row["id"], cls=cls, direction=d, kt=kt, f=f, Y=out, seconds=time.time() - t0, unknowns=int(basis.n_unknowns))


def ckpt(rid):
    return OUT / "ckpt" / f"{rid}.json"


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    rs = rows()
    state = {r["id"]: dict(nodes=G.load_nodes(ckpt(r["id"])), pending=None) for r in rs}
    for r in rs:
        state[r["id"]]["pending"] = [] if state[r["id"]]["nodes"] else [GUESS]
        if state[r["id"]]["nodes"]:
            nk = S.next_node(state[r["id"]]["nodes"], ALL_STATES, SUB)
            state[r["id"]]["pending"] = [nk] if nk is not None else []
    t0 = time.time()
    with Pool(WORKERS, maxtasksperchild=1) as pool:
        for rnd in range(5):
            tasks = [(r, c, d, float(kt), N, FO, S.F) for r in rs for kt in state[r["id"]]["pending"] for c in CLASS_STATES for d in "xy"]
            if not tasks:
                break
            print(f"[{time.time() - t0:6.0f}s] round {rnd}: {len(tasks)} assemblies, {WORKERS} workers", flush=True)
            for k, o in enumerate(pool.imap_unordered(node_task, tasks)):
                s = state[o["row"]]
                s["nodes"].setdefault(o["kt"], {}).setdefault(o["direction"], {}).update(o["Y"])
                G.save_nodes(ckpt(o["row"]), s["nodes"])
                print(f"[{time.time() - t0:6.0f}s] {k + 1}/{len(tasks)} {o['row']} {o['cls']} {o['direction']} kt={o['kt']:.4f} {o['seconds']:.0f}s {o['unknowns']} unk", flush=True)
            for r in rs:
                s = state[r["id"]]
                nk = S.next_node(s["nodes"], ALL_STATES, SUB)
                s["pending"] = [nk] if nk is not None else []
    print("RUN_DONE", flush=True)


def res():
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / "resonance.json"
    done = json.load(open(p)) if p.exists() else {}
    tasks = [(r, c, "x", 0.0, N, FO, fr) for r in rows() for c in CLASS_STATES for fr in (5.0e9, 6.6e9) if f"{r['id']}|{c}|{fr:g}" not in done]
    with Pool(WORKERS, maxtasksperchild=1) as pool:
        for o in pool.imap_unordered(node_task, tasks):
            done[f"{o['row']}|{o['cls']}|{o['f']:g}"] = G.jc(o["Y"])
            p.write_text(json.dumps(done))
    print("RES_DONE", flush=True)


def collect():
    resj = json.load(open(OUT / "resonance.json")) if (OUT / "resonance.json").exists() else {}
    out, lines = [], ["id,theta_deg,q,W,alpha_NODIODE,alpha_OFF,alpha_ON,beta_OFF,beta_ON,absDelta_OFF_mS,absDelta_ON_mS,A_perp_OFF_mS,A_perp_ON_mS,Yuv_over_Delta_ON,rho_ON_over_OFF,rho_ON_over_NODIODE,LEAKAGE,resonance_flip"]
    for r in rows():
        nodes = G.load_nodes(ckpt(r["id"]))
        if not nodes:
            continue
        sm = {st: S.state_summary(nodes, st, SUB, guess=GUESS) for st in ALL_STATES}
        flip = None
        if f"{r['id']}|BRIDGE|5e+09" in resj:
            flip = []
            for cls in CLASS_STATES:
                y5 = {st: G.to_c(v) for st, v in resj[f"{r['id']}|{cls}|5e+09"].items()}
                y6 = {st: G.to_c(v) for st, v in resj[f"{r['id']}|{cls}|6.6e+09"].items()}
                for st in y5:
                    flip += [f"{st}:{'uv'[i]}" for i in (0, 1) if np.sign(y5[st][i, i].imag) != np.sign(y6[st][i, i].imag)]
        yx, yy = S.interp_y(nodes, "ON_BAR50", sm["ON_BAR50"]["beta_over_k0"])
        d_on, d_off, d_nd = sm["ON_BAR50"]["Delta_mS"], sm["OFF_BAR50"]["Delta_mS"], sm["NO_DIODE"]["Delta_mS"]
        yuv = float(max(abs(yx[0, 1]), abs(yy[1, 0])) * 1e3 / abs(d_on)) if abs(d_on) > 0 else float("nan")
        a_on, a_off = sm["ON_BAR50"]["alpha_np_per_m"], sm["OFF_BAR50"]["alpha_np_per_m"]
        cls_ = G.leakage_class(a_on, a_off)
        rec = dict(id=r["id"], theta=r["theta"], q=r["q"], W=r["W"], alpha_NODIODE=sm["NO_DIODE"]["alpha_np_per_m"], alpha_OFF=a_off, alpha_ON=a_on, beta_OFF=sm["OFF_BAR50"]["beta_over_k0"], beta_ON=sm["ON_BAR50"]["beta_over_k0"],
                   Delta_mS={k: v["Delta_mS"] for k, v in sm.items()}, A_perp_mS={k: v["A_perp_mS"] for k, v in sm.items()}, Yuv_over_Delta_ON=yuv, rho_ON_over_OFF=d_on / d_off, rho_ON_over_NODIODE=d_on / d_nd,
                   LEAKAGE=cls_, resonance_flip=flip, nodes=sorted(nodes), extrapolated={k: v["extrapolated"] for k, v in sm.items()})
        out.append(rec)
        lines.append(f"{r['id']},{r['theta']},{r['q']:g},{r['W']:g},{rec['alpha_NODIODE']:.4g},{a_off:.4g},{a_on:.4g},{rec['beta_OFF']:.4g},{rec['beta_ON']:.4g},{abs(d_off):.4g},{abs(d_on):.4g},{rec['A_perp_mS']['OFF_BAR50']:.4g},{rec['A_perp_mS']['ON_BAR50']:.4g},{yuv:.3g},{abs(rec['rho_ON_over_OFF']):.4g},{abs(rec['rho_ON_over_NODIODE']):.4g},{cls_},{'NA' if flip is None else ('|'.join(flip) or 'none')}")
    prov = G.provenance(dict(stage="edge_diode", r_mm=R, gap_mm=GAP))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "edge_diode.json").write_text(json.dumps(G.jc(dict(provenance=prov, results=out)), indent=1, default=str))
    (OUT / "edge_diode.csv").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def selftest():
    """(1) wrapped-mask disc (theta 0, no sheet) == the unshifted slotted disc translated by p/2 (n = 32: exact 16-cell shift): same Y tensor.
    (2) geometry facts: sheet cells, disc-edge gap, closed rim."""
    from src.dispersion.slotted_cells import SlottedCircularPatch
    n, fo = 32, 64
    f = S.F
    k0 = free_space_wavenumber(f)
    L, W = 8.48, 2.9
    sub = S.make_sub(SUB)
    green = LayeredGreenFunction(sub, f, backing="pec")
    ys = []
    for geo in (EdgeDisc(R * MM, L * MM, W * MM, 0.0, P * MM, bridge=False), SlottedCircularPatch(2 * R * MM, L * MM, W * MM)):
        cell = UnitCell(RectangularLattice(P * MM, P * MM), geo, sub, Metal(), backing="pec")
        b = RooftopBasis(cell, n, n)
        kx, ky = 1.06 * k0, 0.0
        z = assemble_impedance_matrix(b, green, kx, ky, fo)
        v = np.stack([assemble_excitation_vector(b, e, kx, ky) for e in np.eye(2)], axis=1)
        j = np.stack([b.current_coefficients(solve(z, v)[:, k], kx, ky) for k in range(2)], axis=1)
        ys.append(j @ np.linalg.inv(np.eye(2) + np.asarray(green.spectral_dyad(kx, ky), dtype=complex) @ j))
    err = np.abs(ys[0] - ys[1]).max() / np.abs(ys[1]).max()
    print("shifted-vs-centred Y relative difference", err)
    assert err < 1e-8, "origin shift changed Y"
    h = P / N
    xc = -P / 2 + (np.arange(N) + 0.5) * h
    X, Y = np.meshgrid(xc * MM, xc * MM, indexing="ij")
    for th in THETAS:
        geo = EdgeDisc(R * MM, 8.48 * MM, 2.9 * MM, math.radians(th), P * MM, bridge=True)
        m = geo.in_feature(X, Y)
        strip = (np.abs(X) <= SHEET_L / 2) & (np.abs(Y) <= SHEET_W / 2)
        gapcols = int(np.sum((m & strip)[:, N // 2]))
        print(f"theta {th}: sheet cells {int(strip.sum())} (6 x 2 expected), metal-cell count {int(m.sum())}")
    from scipy.ndimage import label
    geo = EdgeDisc(R * MM, 8.48 * MM, 2.9 * MM, 0.0, P * MM, bridge=False)
    _, nl = label(geo.in_feature(X, Y))
    print("components without the sheet (periodic wrap not applied in label):", nl)
    print("selftest ok")


if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else ""
    {"selftest": selftest, "run": run, "res": res, "collect": collect}[st]()
