"""SUB6_SPLITDISC_1BIT stage 1 (job `jobs/sub6-splitdisc-cell-scan.md`): split-disc cell scan at 5.8 GHz. MODEL_ONLY, IDEAL_BIAS, NOT_HFSS_VALIDATED.

Cell: a disc of radius r (variant S) or radii r1 / r2 = q r1 (variant T) centred in a square cell of period p, cut by a straight gap g through the centre;
one lumped diode load bridges the gap at the centre (footprint g across x `foot` along the gap = the footprint of the HFSS RLC sheet). Grounded RT5880 slab.

GEOMETRY PATH: `PolyComposite` (src/polygon_geometry.py), halves as arc polygons + the footprint strip as metal. `SlottedCircularPatch` cannot give variant T (two radii).
CONVENTION: a half is the part of the disc of radius r (centred at the cell centre) beyond the gap, |u| >= g/2, u = coordinate across the gap. (The literal
"semicircle with its centre on the flat edge" needs r + g/2 <= p/2 and fails for r > 4.675 mm on the job's r grid; the two shapes differ by O(g^2/r).)

Solve: loaded Bloch sheet tensor at kt = beta (same building blocks as `sidequest_cell_search.bloch_loaded` = `tensor_sheet.bloch_sheet_tensor` + `loads.load_matrix`), but
with the job's Floquet cutoff F = 2n (`bloch_loaded` itself uses F = n). No solver change. Gap frame: a = along the gap (Y_par), u = across (Y_perp).

    python sub6_splitdisc_cell_scan.py pilot | grid | conv | rot | freq | collect   [--workers N]
"""
import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import hashlib
import json
import math
import platform
import subprocess
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
warnings.simplefilter("ignore")
import numpy as np

from src.basis_functions import MeshSettings, RooftopBasis
from src.composite_geometry import AxisRect
from src.constants import MU0
from src.dispersion.sheet_model import _y_air, _y_slab
from src.dispersion.tensor_sheet import _b_te, find_beta_tensor, rotation
from src.errors import LwmomError
from src.frequency import angular_frequency, free_space_wavenumber
from src.geometry import UnitCell
from src.green_function import LayeredGreenFunction
from src.loads import LumpedLoad, default_transverse_conducts, load_matrix
from src.materials import Metal, Substrate
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.periodic_bc import RectangularLattice
from src.polygon_geometry import Poly, PolyComposite, oriented_rect
from src.ris.diodes import DiodeModel, diode_load
from src.solver import solve

MM = 1e-3
P = 10.0  # mm
F0 = 5.8e9
OUT = ROOT / "results/sub6_splitdisc"
CK = OUT / "ckpt"
LABELS = ["MODEL_ONLY", "IDEAL_BIAS", "NOT_HFSS_VALIDATED"]

# ---------------------------------------------------------------- diode states (footprint family, gap g [mm], sheet along-gap width [mm])
FOOT_W = {"MA": 0.35, "SC79": 0.35}
GAP = {"MA": 0.33, "SC79": 0.65}
# name -> (family, kind, params). kind: "open" (no load, empty gap), "C" (series C only), "R" (series R only), "RLC"
STATES = {
    "OFF_OPEN_MA": ("MA", "open", {}),
    "OFF_MA": ("MA", "C", {"C": 28.16e-15}),
    "ON_MA": ("MA", "R", {"R": 2.44}),
    "OFF_OPEN": ("SC79", "open", {}),
    "OFF_SMP": ("SC79", "RLC", {"C": 0.19e-12, "L": 0.7e-9}),
    "OFF_SMP_0V": ("SC79", "RLC", {"C": 0.25e-12, "L": 0.7e-9}),  # ESTIMATE
    "ON_SMP": ("SC79", "RLC", {"R": 1.5, "L": 0.7e-9}),
    "OFF_BAR50": ("SC79", "RLC", {"C": 0.15e-12, "L": 0.6e-9}),
    "OFF_BAR50_HI": ("SC79", "RLC", {"C": 0.20e-12, "L": 0.6e-9}),
    "ON_BAR50": ("SC79", "RLC", {"R": 3.0, "L": 0.6e-9}),
    "ON_BAR50_MAX": ("SC79", "RLC", {"R": 4.5, "L": 0.6e-9}),
}
STATES_MA = ["OFF_OPEN_MA", "OFF_MA", "ON_MA"]
STATES_SC79_S = ["OFF_OPEN", "OFF_BAR50", "OFF_BAR50_HI", "ON_BAR50", "ON_BAR50_MAX", "OFF_SMP", "OFF_SMP_0V", "ON_SMP"]
STATES_BAR50 = ["OFF_BAR50", "OFF_BAR50_HI", "ON_BAR50", "ON_BAR50_MAX"]


def state_load(name, psi_deg=0.0, g_mm=None):
    """`LumpedLoad` of a state at the cell centre, conduction across the gap (psi + 90 deg); None for an open gap."""
    fam, kind, p = STATES[name]
    if kind == "open":
        return None
    g = (GAP[fam] if g_mm is None else g_mm) * MM
    d = "y" if psi_deg % 180.0 == 0.0 else float((psi_deg + 90.0) % 180.0)
    kw = dict(x=0.0, y=0.0, length=g, width=FOOT_W[fam] * MM, direction=d, transverse_conducts=default_transverse_conducts(d))
    if kind == "C":
        return LumpedLoad(C=p["C"], **kw)
    if kind == "R":
        return LumpedLoad(R=p["R"], **kw)
    m = DiodeModel(name, R_s=p.get("R", 0.0), L_s=p.get("L", 0.0), C_T=p.get("C"))
    return diode_load(m, 0.0, 0.0, g, FOOT_W[fam] * MM, d)


# ---------------------------------------------------------------- geometry
def half_disc(r, g, sign, psi_rad, n_arc=96):
    """Poly of the part of the disc of radius r with sign*u >= g/2 (u across the gap), rotated by psi (a, u -> x, y)."""
    t0 = math.asin(min(1.0, g / (2.0 * r)))
    th = np.linspace(t0, math.pi - t0, n_arc + 1)
    pts = [(r * math.cos(t), sign * r * math.sin(t)) for t in th]  # (a, u); closes along the chord u = sign*g/2
    c, s = math.cos(psi_rad), math.sin(psi_rad)
    return Poly(tuple((a * c - u * s, a * s + u * c) for a, u in pts))


def build_cell(variant, r1_mm, q, h_mm, g_mm, foot_mm, psi_deg=0.0, load=None, bridge=True, eps=2.2, tand=0.0009, p_mm=P):
    r1 = r1_mm * MM
    r2 = r1 if variant == "S" else q * r1
    psi = math.radians(psi_deg)
    metal = [half_disc(r1, g_mm * MM, +1, psi), half_disc(r2, g_mm * MM, -1, psi)]
    if bridge:  # footprint strip: conductor for the mesh, the load acts across it
        if psi_deg % 180.0 == 0.0:
            metal.append(AxisRect(0.0, 0.0, foot_mm * MM, g_mm * MM))
        else:
            metal.append(oriented_rect(0.0, 0.0, g_mm * MM, foot_mm * MM, psi + math.pi / 2))
    sub = Substrate("rt5880", eps, tand, h_mm * MM)
    return UnitCell(RectangularLattice(p_mm * MM, p_mm * MM), PolyComposite(metal=tuple(metal)), sub, Metal(), backing="pec", loads=(load,) if load is not None else ())


def cell_for_state(variant, r1_mm, q, h_mm, state, psi_deg=0.0, g_mm=None):
    fam = STATES[state][0]
    g = GAP[fam] if g_mm is None else g_mm
    ld = state_load(state, psi_deg, g)
    return build_cell(variant, r1_mm, q, h_mm, g, FOOT_W[fam], psi_deg, ld, bridge=ld is not None)


# ---------------------------------------------------------------- tensor
def bloch_loaded_f2n(unit_cell, f, n, kx0, ky0, floquet=None):
    """Loaded Bloch sheet tensor (2 x 2, S), F = 2n by default: `tensor_sheet.bloch_sheet_tensor` + `loads.load_matrix` (as `sidequest_cell_search.bloch_loaded`, whose cutoff is n)."""
    mesh = MeshSettings(n, n, 2 * n if floquet is None else floquet)
    green = LayeredGreenFunction(unit_cell.substrate, f, backing="pec")
    basis = RooftopBasis(unit_cell, mesh.cells_x, mesh.cells_y)
    z = assemble_impedance_matrix(basis, green, kx0, ky0, mesh.floquet_order)
    if unit_cell.loads:
        z = z + load_matrix(basis, list(unit_cell.loads), green.omega)
    v = np.stack([assemble_excitation_vector(basis, e, kx0, ky0) for e in np.eye(2)], axis=1)
    cur = solve(z, v)
    j = np.stack([basis.current_coefficients(cur[:, k], kx0, ky0) for k in range(2)], axis=1)
    e_tot = np.eye(2) + np.asarray(green.spectral_dyad(kx0, ky0), dtype=complex) @ j
    return j @ np.linalg.inv(e_tot)


def alpha_abs(y, phi, beta, sub, f):
    """Absorption alpha [Np/m] of the bound mode: Re F / (d Im F / d beta) at the root of F = Y_eff + Y_air + Y_slab (sheet resistance + dielectric loss in eps_c)."""
    om, k0, eps_c, h = angular_frequency(f), free_space_wavenumber(f), sub.eps_complex, sub.thickness
    r = rotation(phi)
    yl = r.T @ np.asarray(y, dtype=complex) @ r
    ykk, yks, yss = yl[0, 0], 0.5 * (yl[0, 1] + yl[1, 0]), yl[1, 1]
    coupled = abs(yks) > 1e-10 * (abs(ykk) + abs(yss))

    def big_f(b):
        ye = ykk - yks**2 / (yss + _b_te(b, om, eps_c, k0, h)) if coupled else ykk
        return ye + _y_air(b, om, k0) + _y_slab(b, om, eps_c, k0, h)

    db = 1e-4 * k0
    d = (big_f(beta + db).imag - big_f(beta - db).imag) / (2 * db)
    return float(abs(big_f(beta).real / d)) if d != 0 else float("nan")


def peak_rss_mb():
    """Peak working set of this process [MB] (Windows psapi; NaN elsewhere)."""
    try:
        import ctypes
        from ctypes import wintypes

        class PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD), ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                        ("a", ctypes.c_size_t), ("b", ctypes.c_size_t), ("c", ctypes.c_size_t), ("d", ctypes.c_size_t), ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]

        pm = PMC()
        pm.cb = ctypes.sizeof(PMC)
        k32, ps = ctypes.WinDLL("kernel32"), ctypes.WinDLL("psapi")
        k32.GetCurrentProcess.restype = wintypes.HANDLE
        ps.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD]
        ps.GetProcessMemoryInfo.restype = wintypes.BOOL
        ps.GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(pm), pm.cb)
        return pm.PeakWorkingSetSize / 2**20
    except Exception:  # noqa: BLE001
        return float("nan")


def solve_row(variant, r1, q, h, state, direction, psi, f, n, g_mm=None, beta0=1.07, tol=5e-4, maxit=6):
    """One row: tensor at kt = beta (self-consistent with the uniform lattice of THIS state) along the gap axis ('par', azimuth psi) or across it ('perp', psi + 90)."""
    c = cell_for_state(variant, r1, q, h, state, psi, g_mm)
    sub = c.substrate
    k0 = free_space_wavenumber(f)
    phi = math.radians(psi + (0.0 if direction == "par" else 90.0))
    beta, status, y, n_solve = beta0 * k0, "ITERATING", None, 0
    t0 = time.time()
    for it in range(1, maxit + 1):
        y = bloch_loaded_f2n(c, f, n, beta * math.cos(phi), beta * math.sin(phi))
        n_solve += 1
        try:
            b_new, _ = find_beta_tensor(y, phi, sub, f)
        except LwmomError:
            status = "BETA_NO_ROOT"
            break
        conv = abs(b_new - beta) < tol * k0
        beta_prev, beta = beta, b_new
        if conv:
            status = "BETA_CONVERGED"
            break
        if it == maxit:
            status = "BETA_NOT_CONVERGED"
    g_used = GAP[STATES[state][0]] if g_mm is None else g_mm
    row = {
        "variant": variant, "r1_mm": r1, "q": q, "h_mm": h, "g_mm": g_used, "foot_mm": FOOT_W[STATES[state][0]], "state": state, "direction": direction, "psi_deg": psi,
        "f_ghz": f / 1e9, "mesh": n, "floquet": 2 * n, "Y_mS": [[complex(v).real * 1e3, complex(v).imag * 1e3] for v in y.ravel()],  # row-major xx xy yx yy (lab frame)
        "beta_over_k0": beta / k0, "beta_status": status, "iterations": n_solve, "gap_mesh_cells": g_used / (P / n),
        "alpha_abs_np_per_m": alpha_abs(y, phi, beta, sub, f) if status != "BETA_NO_ROOT" else None,
        "seconds": time.time() - t0, "peak_rss_mb": peak_rss_mb(), "labels": LABELS,
    }
    return row


# ---------------------------------------------------------------- jobs / checkpoints
def key_of(t):
    return hashlib.sha1(json.dumps(t, sort_keys=True).encode()).hexdigest()[:16]


def run_unit(unit):
    """unit = dict(tag, variant, r1, q, h, psi, f, n, direction, states, [g_mm]); states run in sequence, beta warm-started; every row checkpointed."""
    CK.mkdir(parents=True, exist_ok=True)
    beta0, rows = unit.get("beta0", 1.07), []
    for st in unit["states"]:
        t = {k: unit[k] for k in ("tag", "variant", "r1", "q", "h", "psi", "f", "n", "direction")} | {"state": st, "g_mm": unit.get("g_mm")}
        fn = CK / f"{key_of(t)}.json"
        if fn.exists():
            row = json.loads(fn.read_text(encoding="utf-8"))
        else:
            row = solve_row(unit["variant"], unit["r1"], unit["q"], unit["h"], st, unit["direction"], unit["psi"], unit["f"], unit["n"], unit.get("g_mm"), beta0)
            row["tag"] = unit["tag"]
            tmp = fn.with_suffix(".tmp")
            tmp.write_text(json.dumps(row, default=float), encoding="utf-8", newline="\n")
            tmp.replace(fn)
        if row["beta_status"] != "BETA_NO_ROOT":
            beta0 = row["beta_over_k0"]
        rows.append(row)
    return rows


def units_grid():
    rs = [round(2.0 + 0.2 * i, 1) for i in range(14)]  # 2.0 .. 4.6 (p/2 - 0.3 = 4.7)
    out = []
    for h in (1.575, 3.175):
        for r in rs:
            for d in ("par", "perp"):
                base = dict(tag="grid", variant="S", r1=r, q=1.0, h=h, psi=0.0, f=F0, n=64, direction=d)
                out.append(base | dict(states=STATES_SC79_S))
                out.append(base | dict(states=STATES_MA))
                for q in (0.6, 0.7, 0.8, 0.9):
                    out.append(base | dict(variant="T", q=q, states=STATES_BAR50))
    return out


def _grid_beta0(variant, r, q, h, state, direction):
    """Converged beta/k0 of the same geometry/state/direction at mesh 64 (warm start of the refined-mesh runs); 1.07 if the grid row is missing."""
    t = {"tag": "grid", "variant": variant, "r1": r, "q": q, "h": h, "psi": 0.0, "f": F0, "n": 64, "direction": direction, "state": state, "g_mm": None}
    fn = CK / f"{key_of(t)}.json"
    if fn.exists():
        row = json.loads(fn.read_text(encoding="utf-8"))
        if row["beta_status"] != "BETA_NO_ROOT":
            return row["beta_over_k0"]
    return 1.07


SUB_R = (2.0, 3.4, 4.6)  # small / medium / largest r of variant S on the grid (rule fixed here, before the subsets run)
SUB_H = (1.575, 3.175)
CONV_STATES = ("OFF_MA", "ON_MA", "OFF_OPEN_MA", "OFF_BAR50", "ON_BAR50", "OFF_OPEN")  # MA at g 0.33 (hardest mesh case) and the primary BAR50 at g 0.65 (job correction 2026-11-07)
CONV_GROUPS = {"MA_g0.33": ("OFF_MA", "ON_MA", "OFF_OPEN_MA"), "BAR50_g0.65": ("OFF_BAR50", "ON_BAR50", "OFF_OPEN")}


def units_conv(n):
    return [dict(tag="conv", variant="S", r1=r, q=1.0, h=h, psi=0.0, f=F0, n=n, direction=d, states=[st], beta0=_grid_beta0("S", r, 1.0, h, st, d))
            for r in SUB_R for h in SUB_H for d in ("par", "perp") for st in CONV_STATES]


# ---- variant T: cross term (job correction 2: report, do not stop at H1_VIOLATED)
RECIP = OUT / "recip"
T_STATES = ("OFF_BAR50", "ON_BAR50")


def _offparts(row):
    """(sym, anti) = ((Y_xy + Y_yx)/2, (Y_xy - Y_yx)/2) of a row [S]."""
    y = _mat(row)
    return 0.5 * (y[0, 1] + y[1, 0]), 0.5 * (y[0, 1] - y[1, 0])


def pick_t_geom():
    """T grid geometry (r1, q, h) with the largest |Y_off| (max of sym / anti parts, both kt directions, OFF_BAR50 and ON_BAR50); None if no T rows yet."""
    best = None
    for r in load_rows():
        if r["tag"] != "grid" or r["variant"] != "T" or r["state"] not in T_STATES:
            continue
        sy, an = _offparts(r)
        v = max(abs(sy), abs(an))
        if best is None or v > best[0]:
            best = (v, r["r1_mm"], r["q"], r["h_mm"])
    return None if best is None else best[1:]


def units_tfollow(n):
    g = pick_t_geom()
    if g is None:
        return []
    r, q, h = g
    return [dict(tag="tfollow", variant="T", r1=r, q=q, h=h, psi=0.0, f=F0, n=n, direction=d, states=[st], beta0=_grid_beta0("T", r, q, h, st, d)) for d in ("par", "perp") for st in T_STATES]


def recip_task(a):
    """Y at +k and -k (same |kt| = grid beta of the state/direction) for one T geometry/state/direction: reciprocity needs Y(-k) = Y(k)^T."""
    variant, r, q, h, st, d = a
    RECIP.mkdir(parents=True, exist_ok=True)
    fn = RECIP / f"{variant}_{r}_{q}_{h}_{st}_{d}.json"
    if fn.exists():
        return json.loads(fn.read_text(encoding="utf-8"))
    k0 = free_space_wavenumber(F0)
    beta = _grid_beta0(variant, r, q, h, st, d) * k0
    c = cell_for_state(variant, r, q, h, st)
    phi = math.radians(0.0 if d == "par" else 90.0)
    out = {"variant": variant, "r1_mm": r, "q": q, "h_mm": h, "state": st, "direction": d, "mesh": 64, "beta_over_k0": beta / k0}
    for sign, key in ((1.0, "Y_plus_k_mS"), (-1.0, "Y_minus_k_mS")):
        y = bloch_loaded_f2n(c, F0, 64, sign * beta * math.cos(phi), sign * beta * math.sin(phi))
        out[key] = [[complex(v).real * 1e3, complex(v).imag * 1e3] for v in y.ravel()]
    fn.write_text(json.dumps(out), encoding="utf-8", newline="\n")
    return out


def units_recip():
    g = pick_t_geom()
    return [] if g is None else [("T", g[0], g[1], g[2], st, d) for st in T_STATES for d in ("par", "perp")]


# ---- variant T: dependence of the cross term on the angle between kt and the gap axis (job addition 2026-11-07)
ANGLE_DIR = OUT / "angle"
ANGLES = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0)


def angle_geoms():
    """The T geometry with the largest |Y_off| plus small (2.0 mm) and medium (3.4 mm) r1 at q = 0.6, both h (deduplicated)."""
    g = pick_t_geom()
    out = [] if g is None else [g]
    for r in (2.0, 3.4):
        for h in SUB_H:
            if (r, 0.6, h) not in out:
                out.append((r, 0.6, h))
    return out


def angle_task(a):
    """Full 2 x 2 tensor of a T cell (mesh 64, F = 2n) at kt = |beta| (the state's along-gap beta of the grid) in the direction `angle` from the gap axis."""
    r, q, h, st, ang = a
    ANGLE_DIR.mkdir(parents=True, exist_ok=True)
    fn = ANGLE_DIR / f"T_{r}_{q}_{h}_{st}_{ang}.json"
    if fn.exists():
        return json.loads(fn.read_text(encoding="utf-8"))
    k0 = free_space_wavenumber(F0)
    b_par, b_perp = _grid_beta0("T", r, q, h, st, "par"), _grid_beta0("T", r, q, h, st, "perp")
    beta = b_par * k0
    phi = math.radians(ang)
    y = bloch_loaded_f2n(cell_for_state("T", r, q, h, st), F0, 64, beta * math.cos(phi), beta * math.sin(phi))
    out = {"variant": "T", "r1_mm": r, "q": q, "h_mm": h, "state": st, "angle_deg": ang, "mesh": 64, "beta_used_over_k0": b_par, "beta_perp_grid_over_k0": b_perp,
           "Y_mS": [[complex(v).real * 1e3, complex(v).imag * 1e3] for v in y.ravel()], "labels": LABELS}
    fn.write_text(json.dumps(out), encoding="utf-8", newline="\n")
    return out


def units_angle():
    return [(r, q, h, st, ang) for (r, q, h) in angle_geoms() for st in T_STATES for ang in ANGLES]


def angle_analysis():
    """Per geometry / state / angle: mean (Y_xx+Y_yy)/2, Delta = Y_xx - Y_yy, symmetric off-diagonal S, antisymmetric part A; and the OFF -> ON ratio (|.| and phase) of each."""
    recs = {}
    for f in sorted(ANGLE_DIR.glob("*.json")) if ANGLE_DIR.exists() else []:
        o = json.loads(f.read_text(encoding="utf-8"))
        y = np.array([_cx(v) for v in o["Y_mS"]]).reshape(2, 2)
        recs[(o["r1_mm"], o["q"], o["h_mm"], o["state"], o["angle_deg"])] = {
            "mean": 0.5 * (y[0, 0] + y[1, 1]), "Delta": y[0, 0] - y[1, 1], "S": 0.5 * (y[0, 1] + y[1, 0]), "A": 0.5 * (y[0, 1] - y[1, 0]), "beta": o["beta_used_over_k0"]}
    out = []
    for (r, q, h) in sorted({k[:3] for k in recs}):
        row = {"r1_mm": r, "q": q, "h_mm": h, "angles": []}
        for ang in ANGLES:
            off, on = recs.get((r, q, h, "OFF_BAR50", ang)), recs.get((r, q, h, "ON_BAR50", ang))
            if off is None or on is None:
                continue
            e = {"angle_deg": ang}
            for comp in ("mean", "Delta", "S", "A"):
                a, b = off[comp], on[comp]
                e[comp] = {"OFF_mS": [a.real, a.imag], "ON_mS": [b.real, b.imag], "OFF_abs": abs(a), "ON_abs": abs(b), "ratio_abs_ON_OFF": abs(b) / abs(a) if abs(a) > 1e-12 else None,
                           "phase_change_deg": float(np.degrees(np.angle(b / a))) if abs(a) > 1e-12 and abs(b) > 1e-12 else None}
            e["A_over_Delta_OFF"] = abs(off["A"]) / abs(off["Delta"])
            e["A_over_Delta_ON"] = abs(on["A"]) / abs(on["Delta"])
            row["angles"].append(e)
        out.append(row)
    return out


def units_rot():
    return [dict(tag="rot", variant="S", r1=r, q=1.0, h=h, psi=psi, f=F0, n=64, direction=d, states=["OFF_MA", "ON_MA"], beta0=_grid_beta0("S", r, 1.0, h, "OFF_MA", d))
            for r in SUB_R for h in SUB_H for psi in (22.5, 45.0) for d in ("par", "perp")]


def units_freq():
    return [dict(tag="freq", variant="S", r1=r, q=1.0, h=h, psi=0.0, f=f, n=64, direction=d, states=["OFF_MA", "ON_MA", "OFF_BAR50", "ON_BAR50"], beta0=_grid_beta0("S", r, 1.0, h, "OFF_MA", d))
            for r in SUB_R for h in SUB_H for f in (5.6e9, 6.0e9) for d in ("par", "perp")]


def run_units(units, workers):
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(workers) as ex:
        futs = {ex.submit(run_unit, u): u for u in units}
        for fu in as_completed(futs):
            fu.result()
            done += 1
            print(f"[{time.time() - t0:8.0f}s] unit {done}/{len(units)}", flush=True)


def est_mem_gb(n, r):
    """Peak RAM estimate of one solve: about four n_unknowns x n_unknowns complex matrices (Z, load matrix, their sum, the solver copy) plus 0.6 GB (UNVERIFIED for 96 / 128: the
    rows record peak_rss_mb)."""
    nb = 2 * n * (n - 1) * min(1.0, math.pi * r * r / P**2)
    return 4 * nb * nb * 16 / 1e9 + 0.6


def run_subsets(workers=61, budget_gb=235.0):
    """All follow-up runs on the FULL pool (61 = the Windows ProcessPoolExecutor cap; no reserved cores). Only RAM limits concurrency: a task is started when its estimated
    peak memory fits in `budget_gb` together with the running ones (heaviest first). Every task runs in a fresh process (max_tasks_per_child=1), so memory is returned."""
    from concurrent.futures import FIRST_COMPLETED, wait

    tasks = []  # (est GB, function, argument)
    for n in (128, 96):
        tasks += [(est_mem_gb(n, u["r1"]), run_unit, u) for u in units_tfollow(n) + units_conv(n)]
    tasks += [(est_mem_gb(64, u["r1"]), run_unit, u) for u in units_rot() + units_freq()]
    tasks += [(1.5, recip_task, a) for a in units_recip()]
    tasks += [(1.5, angle_task, a) for a in units_angle()]
    tasks += [(est_mem_leak(a), leak_task, a) for a in units_leak()]
    tasks.sort(key=lambda t: -t[0])
    print(len(tasks), "tasks; est. total", round(sum(t[0] for t in tasks)), "GB-slots; heaviest", round(tasks[0][0], 1), "GB", flush=True)
    t0, running, done = time.time(), {}, 0
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        while tasks or running:
            used = sum(m for m, _ in running.values())
            for t in list(tasks):
                if len(running) >= workers:
                    break
                if used + t[0] <= budget_gb:
                    running[ex.submit(t[1], t[2])] = (t[0], t)
                    used += t[0]
                    tasks.remove(t)
            if not running:
                raise RuntimeError("a single task exceeds the RAM budget")
            fin, _ = wait(list(running), return_when=FIRST_COMPLETED)
            for fu in fin:
                fu.result()
                running.pop(fu)
                done += 1
                print(f"[{time.time() - t0:8.0f}s] done {done}, running {len(running)}, queued {len(tasks)}, est. RAM in use {sum(m for m, _ in running.values()):.0f} GB", flush=True)


# ---------------------------------------------------------------- re-prioritised queue (2026-11-07): h = 3.175 mm first
def _low_priority():
    """Worker initializer: BELOW_NORMAL priority class, so that another job on the same server wins the cores through the OS scheduler."""
    try:
        import ctypes

        k32 = ctypes.WinDLL("kernel32")
        k32.GetCurrentProcess.restype = ctypes.c_void_p
        k32.SetPriorityClass.argtypes = [ctypes.c_void_p, ctypes.c_uint32]
        k32.SetPriorityClass(k32.GetCurrentProcess(), 0x00004000)
    except Exception:  # noqa: BLE001
        pass


def avail_ram_gb():
    try:
        import ctypes

        class MS(ctypes.Structure):
            _fields_ = [("l", ctypes.c_uint32), ("load", ctypes.c_uint32), ("tp", ctypes.c_uint64), ("ap", ctypes.c_uint64), ("tpf", ctypes.c_uint64), ("apf", ctypes.c_uint64), ("tv", ctypes.c_uint64), ("av", ctypes.c_uint64), ("ae", ctypes.c_uint64)]

        m = MS()
        m.l = ctypes.sizeof(MS)
        ctypes.WinDLL("kernel32").GlobalMemoryStatusEx(ctypes.byref(m))
        return m.ap / 2**30
    except Exception:  # noqa: BLE001
        return 100.0


PRIO_R = (3.8, 4.0, 4.2, 4.4, 4.6)
PRIO_S_STATES = ("OFF_BAR50", "ON_BAR50", "OFF_OPEN", "OFF_BAR50_HI", "ON_BAR50_MAX")


def units_prio_s():
    """Item 1: variant S, h = 3.175 mm, r >= 3.8 mm, all BAR50 states and OFF_OPEN, both kt directions; one task per row (heaviest radius first)."""
    out = [dict(tag="grid", variant="S", r1=r, q=1.0, h=3.175, psi=0.0, f=F0, n=64, direction=d, states=[st], beta0=1.03) for r in reversed(PRIO_R) for st in PRIO_S_STATES for d in ("par", "perp")]
    return out


def units_prio_t():
    """Lower priority: variant T, h = 3.175 mm, r1 >= 3.8 mm, BAR50 states (as in the grid)."""
    return [dict(tag="grid", variant="T", r1=r, q=q, h=3.175, psi=0.0, f=F0, n=64, direction=d, states=[st], beta0=1.03)
            for r in reversed(PRIO_R) for q in (0.9, 0.8, 0.7, 0.6) for st in STATES_BAR50 for d in ("par", "perp")]


def run_prio(workers=61, frac_ram=0.7):
    """Item 1 then variant T at h = 3.175 mm, nothing else. Memory-aware (budget = frac_ram of the RAM that is free now), BELOW_NORMAL workers, fresh process per task."""
    from concurrent.futures import FIRST_COMPLETED, wait

    budget = frac_ram * avail_ram_gb()
    tasks = [(est_mem_gb(64, u["r1"]), u) for u in units_prio_s() + units_prio_t()]
    print(f"{len(tasks)} tasks ({len(units_prio_s())} item-1 rows first), RAM budget {budget:.0f} GB, {workers} workers", flush=True)
    t0, running, done = time.time(), {}, 0
    n1 = len(units_prio_s())
    with ProcessPoolExecutor(workers, max_tasks_per_child=1, initializer=_low_priority) as ex:
        while tasks or running:
            used = sum(m for m, _ in running.values())
            for t in list(tasks):
                if len(running) >= workers:
                    break
                if used + t[0] <= budget:
                    running[ex.submit(run_unit, t[1])] = (t[0], t[1])
                    used += t[0]
                    tasks.remove(t)
            if not running:
                raise RuntimeError("a single task exceeds the RAM budget")
            fin, _ = wait(list(running), return_when=FIRST_COMPLETED)
            for fu in fin:
                fu.result()
                running.pop(fu)
                done += 1
                print(f"[{time.time() - t0:8.0f}s] done {done}, running {len(running)}, queued {len(tasks)}", flush=True)
    print("PRIO_DONE", flush=True)


def provenance(extra):
    def sh(*a):
        try:
            return subprocess.run(a, capture_output=True, text=True, cwd=ROOT).stdout.strip()
        except Exception:  # noqa: BLE001
            return "NOT_AVAILABLE"

    import io
    from contextlib import redirect_stdout

    buf = io.StringIO()
    with redirect_stdout(buf):
        np.show_config()
    blas = [ln.strip() for ln in buf.getvalue().splitlines() if "name:" in ln or "version:" in ln][:4]
    return {
        "machine": platform.node(), "os": platform.platform(), "git_commit": sh("git", "rev-parse", "HEAD"), "git_branch": sh("git", "branch", "--show-current"),
        "dirty": bool(sh("git", "status", "--porcelain", "--untracked-files=no")), "python": sys.version, "numpy": np.__version__, "blas": blas,
        "threads": {k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")}, "labels": LABELS,
        "time_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), **extra,
    }


# ---------------------------------------------------------------- derived quantities (used by collect)
def gap_frame_tensor(rows_par, rows_perp, psi_deg):
    """(Y_par, Y_perp, Y_off_a, Y_off_b) in the gap frame from the two kt runs (S): par from the run along the gap axis, perp from the run across it."""
    r = rotation(math.radians(psi_deg))
    ya = r.T @ _mat(rows_par) @ r
    yb = r.T @ _mat(rows_perp) @ r
    return ya[0, 0], yb[1, 1], ya[0, 1], yb[1, 0]


def _mat(row):
    return np.array([complex(a, b) for a, b in row["Y_mS"]]).reshape(2, 2) * 1e-3


# ---------------------------------------------------------------- leakage of a hologram row (job addition 2026-11-07)
# Method of scripts/research/phase26_supercell_table.py (supercell of N cells with the axis rotating by 180 i / N deg along x, complex root of det(Y_s(beta) + diag(A, B)) by
# ANALYTIC CONTINUATION from real-beta evaluations, alpha_leak = alpha_row - alpha_uniform), ported to the split-disc cell (that script is hard-wired to the 10.5 GHz slotted cell).
# Differences, all stated: 9 real points / polynomial degrees 7 and 6 (not 7 / 6 / 5); the points are independent tasks (parallel); mesh 32 per cell (not 64: the N = 9 row has
# 9x the unknowns); Floquet orders (2 N m, 2 m) = the job's 2n rule; a small step (0.0012 k0) because the light-line branch point of Y_air sits only (beta/k0 - 1) k0 from the root.
LEAK_DIR = OUT / "leak"
LEAK_R = (2.0, 3.0, 3.6, 4.2, 4.6)
LEAK_STATES = ("OFF_BAR50", "ON_BAR50")
LEAK_M = 32
LEAK_STEP = 0.0012
LEAK_NPTS = 9
SIN25 = math.sin(math.radians(25.0))


def leak_n_cells(beta_over_k0):
    """Cells per modulation period for the 25 deg beam: lambda_mod = lambda0 / (beta/k0 - sin 25), rounded to a whole number of 10 mm cells (the supercell needs an integer)."""
    lam0 = 299792458.0 / F0
    return max(3, round(lam0 / (beta_over_k0 - SIN25) / (P * MM)))


def _shift(poly, dx):
    return Poly(tuple((x + dx, y) for x, y in poly.vertices))


def build_row(r_mm, h_mm, state, n_cells, psis=None):
    """Supercell UnitCell: n_cells split-disc cells (variant S, pitch 10 mm), cell i rotated by psi_i = 180 i / n_cells deg; one diode load per cell."""
    import dataclasses

    fam = STATES[state][0]
    g, w = GAP[fam], FOOT_W[fam]
    metal, loads = [], []
    for i in range(n_cells):
        psi = 180.0 * i / n_cells if psis is None else psis[i]
        cx = (-n_cells * P / 2 + (i + 0.5) * P) * MM
        c = build_cell("S", r_mm, 1.0, h_mm, g, w, psi, state_load(state, psi, g), bridge=True)
        metal += [_shift(q, cx) if isinstance(q, Poly) else AxisRect(q.cx + cx, q.cy, q.length, q.width) for q in c.geometry.metal]
        loads += [dataclasses.replace(ld, x=ld.x + cx) for ld in c.loads]
    sub = Substrate("rt5880", 2.2, 0.0009, h_mm * MM)
    return UnitCell(RectangularLattice(n_cells * P * MM, P * MM), PolyComposite(metal=tuple(metal)), sub, Metal(), backing="pec", loads=tuple(loads))


def f_det(uc, f, ncx, m, beta, orders):
    """(det(Y_s + diag(A, B)) complex, Y_s) at real beta along x for a supercell `uc` of ncx x 1 cells at mesh m per cell."""
    from src.dispersion.tensor_sheet import assemble_impedance_matrix_orders

    mesh = MeshSettings(ncx * m, m, max(orders))
    green = LayeredGreenFunction(uc.substrate, f, backing="pec")
    basis = RooftopBasis(uc, mesh.cells_x, mesh.cells_y)
    z = assemble_impedance_matrix_orders(basis, green, beta, 0.0, *orders)
    if uc.loads:
        z = z + load_matrix(basis, list(uc.loads), green.omega)
    v = np.stack([assemble_excitation_vector(basis, e, beta, 0.0) for e in np.eye(2)], axis=1)
    cur = solve(z, v)
    j = np.stack([basis.current_coefficients(cur[:, k], beta, 0.0) for k in range(2)], axis=1)
    y = j @ np.linalg.inv(np.eye(2) + np.asarray(green.spectral_dyad(beta, 0.0), dtype=complex) @ j)
    om, k0, sb = green.omega, free_space_wavenumber(f), uc.substrate
    a = _y_air(beta, om, k0) + _y_slab(beta, om, sb.eps_complex, k0, sb.thickness)
    b = _b_te(beta, om, sb.eps_complex, k0, sb.thickness)
    return complex(np.linalg.det(y + np.array([[a, 0.0], [0.0, b]]))), y


def leak_guess(r, h, st):
    """beta/k0 of the mean-sheet uniform lattice (grid rows, mesh 64): the centre of the evaluation window."""
    from src.dispersion import sheet_model as sm

    rows = {d: json.loads((CK / f"{key_of({'tag': 'grid', 'variant': 'S', 'r1': r, 'q': 1.0, 'h': h, 'psi': 0.0, 'f': F0, 'n': 64, 'direction': d, 'state': st, 'g_mm': None})}.json").read_text(encoding="utf-8")) for d in ("par", "perp")}
    yp, yq = _mat(rows["par"])[0, 0], _mat(rows["perp"])[1, 1]
    sub = Substrate("rt5880", 2.2, 0.0009, h * MM)
    return sm.find_beta_from_sheet_admittance(1j * 0.5 * (yp + yq).imag, sub, F0).beta_rad_per_m / free_space_wavenumber(F0), rows["par"]["beta_over_k0"], rows["perp"]["beta_over_k0"]


def leak_points(r, h, st):
    bg, bp, bq = leak_guess(r, h, st)
    n = leak_n_cells(bg)
    return bg, n, [bg + LEAK_STEP * k for k in range(-(LEAK_NPTS // 2), LEAK_NPTS // 2 + 1)]


def leak_task(a):
    """One real-beta evaluation: kind 'row' (N-cell rotating row) or 'uni' (the single cell at psi = 0 = the uniform lattice)."""
    kind, r, h, st, k, m = a
    LEAK_DIR.mkdir(parents=True, exist_ok=True)
    fn = LEAK_DIR / f"{kind}_{r}_{h}_{st}_{k}_{m}.json"
    if fn.exists():
        return json.loads(fn.read_text(encoding="utf-8"))
    bg, n, pts = leak_points(r, h, st)
    beta = pts[k] * free_space_wavenumber(F0)
    t0 = time.time()
    if kind == "row":
        uc, ncx = build_row(r, h, st, n), n
    else:
        uc, ncx = cell_for_state("S", r, 1.0, h, st), 1
    fv, y = f_det(uc, F0, ncx, m, beta, (2 * ncx * m, 2 * m))
    out = {"kind": kind, "r1_mm": r, "h_mm": h, "state": st, "k": k, "mesh_per_cell": m, "n_cells": n, "beta_over_k0": pts[k], "beta_guess_over_k0": bg, "f": [fv.real, fv.imag],
           "Y_mS": [[complex(v).real * 1e3, complex(v).imag * 1e3] for v in y.ravel()], "seconds": time.time() - t0, "peak_rss_mb": peak_rss_mb(), "labels": LABELS}
    fn.write_text(json.dumps(out), encoding="utf-8", newline="\n")
    return out


def units_leak(m=LEAK_M):
    out = []
    for h in SUB_H:
        for r in LEAK_R:
            for st in LEAK_STATES:
                for k in range(LEAK_NPTS):
                    out.append(("row", r, h, st, k, m))
                    out.append(("uni", r, h, st, k, m))
    return out


def est_mem_leak(a):
    kind, r, h, st, k, m = a
    if kind == "uni":
        return 0.6 + 4 * (2 * m * (m - 1) * min(1.0, math.pi * r * r / P**2)) ** 2 * 16 / 1e9
    n = 9
    nb = 2 * n * m * (m - 1) * min(1.0, math.pi * r * r / P**2)
    return 0.6 + 4 * nb * nb * 16 / 1e9


def _complex_root(betas, fv):
    """Complex root of f nearest the real-axis zero of Im f, by polynomial fits of degree 7 and 6 through the points (centred and scaled); returns (root over k0, spread over k0, real root over k0)."""
    x0, hh = betas[len(betas) // 2], betas[1] - betas[0]
    t = (np.asarray(betas) - x0) / hh
    res = []
    for d in (LEAK_NPTS - 2, LEAK_NPTS - 3):
        pol = np.polynomial.Polynomial.fit(t, fv, d)
        rts = pol.roots()
        res.append(rts[np.argmin(np.abs(rts))])
    pim = np.polynomial.Polynomial.fit(t, np.imag(fv), LEAK_NPTS - 3)
    rr = [x.real for x in pim.roots() if abs(x.imag) < 1e-6]
    real_root = x0 + hh * (min(rr, key=abs) if rr else float("nan"))
    root = x0 + hh * res[0]
    return root, abs(res[0] - res[1]) * hh, real_root


def leak_analysis():
    """Per (h, r, state): row and uniform complex roots, alpha_leak = alpha_row - alpha_uniform, beta, realised beam angle, validity flags; plus the first-order closed form."""
    from src.dispersion.leakage_model import leakage_alpha_first_order_tensor

    k0 = free_space_wavenumber(F0)
    out = []
    for h in SUB_H:
        for r in LEAK_R:
            for st in LEAK_STATES:
                rec = {"h_mm": h, "r1_mm": r, "state": st}
                try:
                    bg, n, pts = leak_points(r, h, st)
                except Exception as e:  # noqa: BLE001
                    rec["status"] = f"NO_GRID_ROW ({e.__class__.__name__})"
                    out.append(rec)
                    continue
                rec.update({"n_cells": n, "beta_guess_over_k0": bg, "supercell_period_mm": n * P})
                fr = {}
                for kind in ("row", "uni"):
                    pts_done = []
                    for k in range(LEAK_NPTS):
                        fn = LEAK_DIR / f"{kind}_{r}_{h}_{st}_{k}_{LEAK_M}.json"
                        if fn.exists():
                            o = json.loads(fn.read_text(encoding="utf-8"))
                            pts_done.append((o["beta_over_k0"] * k0, complex(*o["f"])))
                    if len(pts_done) == LEAK_NPTS:
                        b, fv = zip(*pts_done, strict=True)
                        root, spread, rreal = _complex_root(np.array(b), np.array(fv))
                        fr[kind] = {"beta_over_k0": root.real / k0, "alpha_np_per_m": float(-root.imag), "real_root_over_k0": rreal / k0, "spread_over_k0": spread / k0}
                rec["row"], rec["uniform"] = fr.get("row"), fr.get("uni")
                if len(fr) == 2:
                    b = fr["row"]["beta_over_k0"]
                    rec["alpha_leak_np_per_m"] = fr["row"]["alpha_np_per_m"] - fr["uni"]["alpha_np_per_m"]
                    rec["beta_over_k0"] = b
                    rec["beam_angle_deg"] = math.degrees(math.asin(min(1.0, b - 299792458.0 / F0 / (n * P * MM))))
                    off_c = abs(fr["row"]["real_root_over_k0"] - bg) / LEAK_STEP
                    flags = []
                    if not np.isfinite(fr["row"]["real_root_over_k0"]) or off_c > 2.5:
                        flags.append("ROOT_OFF_WINDOW")
                    if (b - 1.0) < 3 * LEAK_STEP:
                        flags.append("BRANCH_POINT_NEAR")
                    if abs(fr["row"]["alpha_np_per_m"]) / k0 > 3 * LEAK_STEP:
                        flags.append("ROOT_OUTSIDE_FIT_WINDOW")
                    if fr["row"]["spread_over_k0"] > 1e-4:
                        flags.append("POLY_DEGREE_UNSTABLE")
                    rec["flags"] = flags or ["OK"]
                    rec["status"] = "SUPERCELL_MOM_COMPUTED"
                else:
                    rec["status"] = "NOT_COMPLETE"
                # first-order closed form from the grid tensors (mesh 64): same beam angle as the supercell if available, else 25 deg
                try:
                    from src.dispersion import sheet_model as sm

                    rows = {d: json.loads((CK / f"{key_of({'tag': 'grid', 'variant': 'S', 'r1': r, 'q': 1.0, 'h': h, 'psi': 0.0, 'f': F0, 'n': 64, 'direction': d, 'state': st, 'g_mm': None})}.json").read_text(encoding="utf-8")) for d in ("par", "perp")}
                    yp, yq = _mat(rows["par"])[0, 0], _mat(rows["perp"])[1, 1]
                    sub = Substrate("rt5880", 2.2, 0.0009, h * MM)
                    beta = sm.find_beta_from_sheet_admittance(1j * 0.5 * (yp + yq).imag, sub, F0).beta_rad_per_m
                    rec["closed_form"] = {"alpha_first_order_np_per_m_25deg": leakage_alpha_first_order_tensor(beta, 1j * 0.5 * (yp + yq).imag, abs(0.5 * (yp - yq)), math.radians(25.0), sub, F0),
                                          "beta_mean_over_k0": beta / k0, "delta_half_mS": abs(0.5 * (yp - yq)) * 1e3, "method": "leakage_alpha_first_order_tensor, local tensors at mesh 64 (OFF/ON uniform-lattice tensors)"}
                except Exception as e:  # noqa: BLE001
                    rec["closed_form"] = f"NOT_AVAILABLE ({e.__class__.__name__})"
                out.append(rec)
    return out


# ---------------------------------------------------------------- collect: derived quantities and gates
DIODES = {  # family -> (OFF state, ON state, role)
    "BAR50": ("OFF_BAR50", "ON_BAR50", "primary (user's choice)"),
    "BAR50_WORST": ("OFF_BAR50_HI", "ON_BAR50_MAX", "sensitivity: OFF C 0.20 pF, ON R 4.5 ohm"),
    "MA": ("OFF_MA", "ON_MA", "reference only (user-measured, over budget)"),
    "SMP": ("OFF_SMP", "ON_SMP", "in-production fallback"),
    "SMP_0V": ("OFF_SMP_0V", "ON_SMP", "fallback, 0 V OFF capacitance ESTIMATE"),
}
SYM_TOL, ROT_TOL, CONV_BETA_TOL, CONV_DELTA_TOL = 0.02, 0.03, 0.005, 0.05


def load_rows():
    return [json.loads(f.read_text(encoding="utf-8")) for f in sorted(CK.glob("*.json"))]


def _cx(v):
    return complex(v[0], v[1])


def tensor_entry(rp, rq):
    """Gap-frame summary of one (geometry, state) from the 'par' and 'perp' rows."""
    yp, yq, oa, ob = gap_frame_tensor(rp, rq, rp["psi_deg"])
    dl = yp - yq
    return {
        "Y_par_mS": [yp.real * 1e3, yp.imag * 1e3], "Y_perp_mS": [yq.real * 1e3, yq.imag * 1e3], "Y_off_par_run_mS": [oa.real * 1e3, oa.imag * 1e3], "Y_off_perp_run_mS": [ob.real * 1e3, ob.imag * 1e3],
        "mean_mS": [0.5 * (yp + yq).real * 1e3, 0.5 * (yp + yq).imag * 1e3], "Delta_mS": [dl.real * 1e3, dl.imag * 1e3],
        "off_over_delta": max(abs(oa), abs(ob)) / abs(dl) if abs(dl) > 0 else None, "beta_par_over_k0": rp["beta_over_k0"], "beta_perp_over_k0": rq["beta_over_k0"],
        "alpha_par_np_per_m": rp["alpha_abs_np_per_m"], "alpha_perp_np_per_m": rq["alpha_abs_np_per_m"], "beta_status": [rp["beta_status"], rq["beta_status"]],
        "gap_mesh_cells": rp["gap_mesh_cells"], "gap_mesh_flag": "GAP_UNDER_1p5_MESH_CELLS" if rp["gap_mesh_cells"] < 1.5 else "OK", "iterations": [rp["iterations"], rq["iterations"]],
    }


def group_rows(rows):
    """{(tag, variant, r1, q, h, psi, f_ghz, mesh): {state: tensor_entry}} for every state that has both direction rows."""
    idx = {}
    for r in rows:
        idx.setdefault((r["tag"], r["variant"], r["r1_mm"], r["q"], r["h_mm"], r["psi_deg"], round(r["f_ghz"], 6), r["mesh"]), {}).setdefault(r["state"], {})[r["direction"]] = r
    out = {}
    for k, sts in idx.items():
        out[k] = {st: tensor_entry(d["par"], d["perp"]) for st, d in sts.items() if "par" in d and "perp" in d}
    return out


def _rho(ent, off, on):
    if off not in ent or on not in ent:
        return None
    d0, d1 = _cx(ent[off]["Delta_mS"]), _cx(ent[on]["Delta_mS"])
    return d1 / d0 if abs(d0) > 0 else None


def derived(ent):
    out = {}
    for fam, (off, on, role) in DIODES.items():
        rho = _rho(ent, off, on)
        if rho is None:
            continue
        a, b = ent[off], ent[on]
        dm = _cx(b["mean_mS"]) - _cx(a["mean_mS"])
        out[fam] = {
            "rho_re": rho.real, "rho_im": rho.imag, "rho_abs": abs(rho), "role": role, "delta_mean_mS": [dm.real, dm.imag],
            "delta_beta_par": b["beta_par_over_k0"] - a["beta_par_over_k0"], "delta_beta_perp": b["beta_perp_over_k0"] - a["beta_perp_over_k0"],
            "alpha_on_par": b["alpha_par_np_per_m"], "alpha_on_perp": b["alpha_perp_np_per_m"],
        }
    return out


def _rot_err(e_psi, e_0):
    """Gap-frame deviation ||(Y_par, Y_perp)(psi) - (Y_par, Y_perp)(0)|| / ||(Y_par, Y_perp)(0)|| (the rotated cell's tensor, rotated back, against the psi = 0 cell)."""
    a = np.array([_cx(e_psi["Y_par_mS"]), _cx(e_psi["Y_perp_mS"])])
    b = np.array([_cx(e_0["Y_par_mS"]), _cx(e_0["Y_perp_mS"])])
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def gates(grp):
    g = {}
    bad, worst, n_rows = [], {(v, d): 0.0 for v in ("S", "T") for d in ("along_gap", "across_gap")}, 0
    for (tag, var, r, q, h, psi, f, n), ent in grp.items():
        if tag != "grid":
            continue
        for st, e in ent.items():
            n_rows += 1
            for lab, key in (("along_gap", "Y_off_par_run_mS"), ("across_gap", "Y_off_perp_run_mS")):
                dl = abs(_cx(e["Delta_mS"]))
                v = abs(_cx(e[key])) / dl if dl > 0 else None
                if v is None:
                    continue
                worst[(var, lab)] = max(worst[(var, lab)], v)
                if v > SYM_TOL:
                    bad.append({"variant": var, "direction": lab, "r1_mm": r, "q": q, "h_mm": h, "state": st, "off_over_delta": v})
    by_var = {}
    for v in ("S", "T"):
        al = any(b["variant"] == v and b["direction"] == "along_gap" for b in bad)
        ac = any(b["variant"] == v and b["direction"] == "across_gap" for b in bad)
        by_var[v] = "H1_SUPPORTED" if not (al or ac) else ("H1_HOLDS_ACROSS_GAP / H1_VIOLATED_ALONG_GAP" if al and not ac else "H1_VIOLATED")
    g["SYMMETRY"] = {"by_variant": by_var, "state": "H1_SUPPORTED" if not bad else "H1_VIOLATED", "tolerance": SYM_TOL, "geometry_states_checked": n_rows,
                     "worst": {f"{v}_{d}": x for (v, d), x in worst.items()}, "violations": bad[:200], "n_violations": len(bad)}
    rot = []
    for (tag, var, r, q, h, psi, f, n), ent in grp.items():
        if tag != "rot":
            continue
        base = grp.get(("grid", var, r, q, h, 0.0, round(F0 / 1e9, 6), 64), {})
        rot += [{"r1_mm": r, "h_mm": h, "psi_deg": psi, "state": st, "err": _rot_err(e, base[st])} for st, e in ent.items() if st in base]
    mx = max((x["err"] for x in rot), default=None)
    g["ROTATION"] = {"state": "NOT_RUN" if mx is None else ("ROTATION_COVARIANT" if mx <= ROT_TOL else "ROTATION_DEVIATES"), "tolerance": ROT_TOL, "max_error": mx, "rows": rot}
    return g


def convergence(grp):
    out = []
    for (tag, var, r, q, h, psi, f, n), ent in grp.items():
        if tag != "conv" or n != 128:
            continue
        for st, e in ent.items():
            e96 = grp.get(("conv", var, r, q, h, psi, f, 96), {}).get(st)
            e64 = grp.get(("grid", var, r, q, h, 0.0, round(F0 / 1e9, 6), 64), {}).get(st)
            row = {"r1_mm": r, "h_mm": h, "state": st, "group": next(k for k, v in CONV_GROUPS.items() if st in v)}
            if e96:
                row["dbeta_par_96_128"] = abs(e["beta_par_over_k0"] / e96["beta_par_over_k0"] - 1)
                row["dbeta_perp_96_128"] = abs(e["beta_perp_over_k0"] / e96["beta_perp_over_k0"] - 1)
                row["dDelta_96_128"] = abs(_cx(e["Delta_mS"]) - _cx(e96["Delta_mS"])) / abs(_cx(e["Delta_mS"]))
                if e64:
                    row["dbeta_par_64_96"] = abs(e96["beta_par_over_k0"] / e64["beta_par_over_k0"] - 1)
                    row["dDelta_64_96"] = abs(_cx(e96["Delta_mS"]) - _cx(e64["Delta_mS"])) / abs(_cx(e96["Delta_mS"]))
            out.append(row)
    res = {"rows": out, "tolerances": {"beta": CONV_BETA_TOL, "Delta": CONV_DELTA_TOL}, "groups": {}}
    for name in CONV_GROUPS:
        ok = [x for x in out if x["group"] == name and "dDelta_96_128" in x]
        bad = [x for x in ok if max(x["dbeta_par_96_128"], x["dbeta_perp_96_128"]) > CONV_BETA_TOL or x["dDelta_96_128"] > CONV_DELTA_TOL]
        res["groups"][name] = {"state": "NOT_RUN" if not ok else ("MESH_CONVERGED" if not bad else "MESH_NOT_CONVERGED"), "n_rows": len(ok), "n_failing": len(bad),
                               "worst_dbeta": max((max(x["dbeta_par_96_128"], x["dbeta_perp_96_128"]) for x in ok), default=None), "worst_dDelta": max((x["dDelta_96_128"] for x in ok), default=None)}
    states = [v["state"] for v in res["groups"].values()]
    res["state"] = "NOT_RUN" if "NOT_RUN" in states else ("MESH_CONVERGED" if all(x == "MESH_CONVERGED" for x in states) else "MESH_NOT_CONVERGED")
    res["n_failing"] = sum(v["n_failing"] for v in res["groups"].values())
    return res


def t_crossterm(rows):
    """Variant T cross term: (a) Y_off along / across the gap, (b) symmetric and antisymmetric parts, (c) reciprocity at the worst geometry, (d) mesh trend, (e) OFF -> ON change."""
    parts, by = [], {}
    for r in rows:
        if r["tag"] == "grid" and r["variant"] == "T":
            by.setdefault((r["r1_mm"], r["q"], r["h_mm"], r["state"]), {})[r["direction"]] = r
    for (r1, q, h, st), d in sorted(by.items()):
        if "par" not in d or "perp" not in d:
            continue
        e = {"r1_mm": r1, "q": q, "h_mm": h, "state": st}
        dl = abs(_mat(d["par"])[0, 0] - _mat(d["perp"])[1, 1])
        for lab, key in (("along_gap", "par"), ("across_gap", "perp")):
            m = _mat(d[key])
            sy, an = 0.5 * (m[0, 1] + m[1, 0]), 0.5 * (m[0, 1] - m[1, 0])
            e[lab] = {"Y_xy_mS": [m[0, 1].real * 1e3, m[0, 1].imag * 1e3], "Y_yx_mS": [m[1, 0].real * 1e3, m[1, 0].imag * 1e3], "sym_mS": [sy.real * 1e3, sy.imag * 1e3],
                      "anti_mS": [an.real * 1e3, an.imag * 1e3], "max_over_Delta": max(abs(m[0, 1]), abs(m[1, 0])) / dl if dl > 0 else None,
                      "anti_over_Delta": abs(an) / dl if dl > 0 else None, "sym_over_Delta": abs(sy) / dl if dl > 0 else None}
        parts.append(e)
    mx = {k: max((x[k]["max_over_Delta"] or 0.0 for x in parts), default=None) for k in ("along_gap", "across_gap")}
    res = {"rows": parts, "max_Yoff_over_Delta": mx}
    if parts:
        a_ok, c_ok = mx["across_gap"] <= SYM_TOL, mx["along_gap"] <= SYM_TOL
        res["state"] = "H1_HOLDS" if a_ok and c_ok else ("H1_HOLDS_ACROSS_GAP / H1_VIOLATED_ALONG_GAP" if a_ok else "H1_VIOLATED_ACROSS_GAP" + ("" if c_ok else " AND ALONG_GAP"))
    else:
        res["state"] = "NOT_RUN"
    idx = {(x["r1_mm"], x["q"], x["h_mm"], x["state"]): x for x in parts}
    ch = []
    for (r1, q, h, st), x in idx.items():
        if st != "OFF_BAR50" or (r1, q, h, "ON_BAR50") not in idx:
            continue
        y = idx[(r1, q, h, "ON_BAR50")]
        row = {"r1_mm": r1, "q": q, "h_mm": h}
        for lab in ("along_gap", "across_gap"):
            for part in ("sym_mS", "anti_mS"):
                a, b = _cx(x[lab][part]), _cx(y[lab][part])
                k = f"{lab}_{part[:-3]}"
                row[k + "_OFF_abs_mS"], row[k + "_ON_abs_mS"] = abs(a), abs(b)
                row[k + "_ratio_ON_OFF"] = abs(b) / abs(a) if abs(a) > 1e-12 else None
                row[k + "_phase_change_deg"] = float(np.degrees(np.angle(b / a))) if abs(a) > 1e-12 and abs(b) > 1e-12 else None
        ch.append(row)
    res["off_to_on"] = ch
    rc = []
    for f in sorted(RECIP.glob("*.json")) if RECIP.exists() else []:
        o = json.loads(f.read_text(encoding="utf-8"))
        yp = np.array([_cx(v) for v in o["Y_plus_k_mS"]]).reshape(2, 2)
        ym = np.array([_cx(v) for v in o["Y_minus_k_mS"]]).reshape(2, 2)
        rc.append({k: o[k] for k in ("r1_mm", "q", "h_mm", "state", "direction")} | {"rel_err_Yminusk_vs_YplusT": float(np.linalg.norm(ym - yp.T) / np.linalg.norm(yp)),
                                                                                      "rel_err_Yminusk_vs_Yplus": float(np.linalg.norm(ym - yp) / np.linalg.norm(yp))})
    res["reciprocity"] = rc
    g, tr = pick_t_geom(), []
    for r in rows:
        if r["tag"] in ("tfollow", "grid") and r["variant"] == "T" and g and (r["r1_mm"], r["q"], r["h_mm"]) == g and r["state"] in T_STATES:
            sy, an = _offparts(r)
            tr.append({"state": r["state"], "direction": r["direction"], "mesh": r["mesh"], "sym_mS": [sy.real * 1e3, sy.imag * 1e3], "anti_mS": [an.real * 1e3, an.imag * 1e3]})
    res["mesh_trend"] = {"geometry_r1_q_h": g, "rows": sorted(tr, key=lambda x: (x["state"], x["direction"], x["mesh"]))}
    return res


def resonance(grp):
    """3-point extremum test (5.6 / 5.8 / 6.0 GHz) of Im Y_par and Im Y_perp on the frequency subset; the grid itself is single-frequency (NOT_EVALUATED there)."""
    flags = []
    for (tag, var, r, q, h, psi, f, n), ent in grp.items():
        if tag != "freq" or abs(f - 5.6) > 1e-6:
            continue
        for st in ent:
            es = [grp.get(("freq", var, r, q, h, psi, 5.6, 64), {}).get(st), grp.get(("grid", var, r, q, h, 0.0, 5.8, 64), {}).get(st), grp.get(("freq", var, r, q, h, psi, 6.0, 64), {}).get(st)]
            if any(e is None for e in es):
                continue
            for comp in ("Y_par_mS", "Y_perp_mS"):
                v = np.array([_cx(e[comp]).imag for e in es])
                if (v[1] - v[0]) * (v[2] - v[1]) < 0 or v[0] * v[2] < 0:
                    flags.append({"r1_mm": r, "h_mm": h, "state": st, "component": comp, "ImY_mS_5p6_5p8_6p0": v.tolist(), "flag": "NEAR_RESONANT"})
    return {"state": "NEAR_RESONANT_FOUND" if flags else "NO_EXTREMUM_IN_5.6-6.0_ON_SUBSET", "scope": "3 geometries x 2 h x 4 states, three frequencies; grid rows: NOT_EVALUATED (single frequency)", "flags": flags}


def collect():
    rows = load_rows()
    grp = group_rows(rows)
    geoms = []
    for (tag, var, r, q, h, psi, f, n), ent in sorted(grp.items(), key=lambda kv: str(kv[0])):
        geoms.append({"tag": tag, "variant": var, "r1_mm": r, "q": q, "h_mm": h, "psi_deg": psi, "f_ghz": f, "mesh": n, "states": ent, "derived": derived(ent)})
    g = gates(grp)
    g["CONVERGENCE"] = convergence(grp)
    g["RESONANCE"] = resonance(grp)
    g["T_CROSSTERM"] = t_crossterm(rows)
    g["LEAKAGE"] = {"state": "COMPUTED" if LEAK_DIR.exists() else "NOT_RUN", "method": "see header of the leakage block; supercell analytic continuation (phase26_supercell_table method) + first-order closed form", "rows": leak_analysis()}
    g["T_ANGLE"] = {"state": "NOT_RUN" if not ANGLE_DIR.exists() else "COMPUTED", "beta_rule": "|kt| = along-gap beta/k0 of the state (grid), angle from the gap axis", "geometries": angle_analysis()}
    out = {"labels": LABELS, "n_rows": len(rows), "geometries": geoms, "gates": g, "provenance": provenance({"job": "sub6-splitdisc-cell-scan", "n_rows": len(rows)})}
    (OUT / "cell_scan.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
    return out


def _pilot_job(a):
    d, st = a
    return solve_row("S", 4.0, 1.0, 3.175, st, d, 0.0, F0, 64)


def pilot():
    """Variant S, r = 4.0 mm, h = 3.175 mm, BAR50 OFF/ON, both kt directions, mesh 64 (4 rows in parallel); writes results/sub6_splitdisc/pilot.json."""
    jobs = [(d, st) for d in ("par", "perp") for st in ("OFF_BAR50", "ON_BAR50")]
    with ProcessPoolExecutor(len(jobs)) as ex:
        rows = list(ex.map(_pilot_job, jobs))
    for r in rows:
        print(r["direction"], r["state"], f"{r['seconds']:.0f}s iters={r['iterations']} -> {r['seconds'] / r['iterations']:.0f}s/solve beta/k0={r['beta_over_k0']:.4f} {r['beta_status']} rss={r['peak_rss_mb']:.0f}MB", flush=True)
    by = {(r["direction"], r["state"]): r for r in rows}
    for st in ("OFF_BAR50", "ON_BAR50"):
        yp, yq, oa, ob = gap_frame_tensor(by[("par", st)], by[("perp", st)], 0.0)
        dl = yp - yq
        print(st, f"Y_par={yp * 1e3:.4f} mS Y_perp={yq * 1e3:.4f} mS Delta={dl * 1e3:.4f} |Yoff|/|Delta|: {abs(oa) / abs(dl):.2e} {abs(ob) / abs(dl):.2e}", flush=True)
    (OUT / "pilot.json").write_text(json.dumps({"rows": rows, "provenance": provenance({"job": "pilot"})}, indent=1, default=float), encoding="utf-8", newline="\n")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "pilot"
    w = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 1
    OUT.mkdir(parents=True, exist_ok=True)
    if cmd == "grid":
        us = units_grid()
        print(len(us), "units", sum(len(u["states"]) for u in us), "rows", flush=True)
        run_units(us, w)
    elif cmd == "prio":
        run_prio()
    elif cmd == "subsets":
        run_subsets()
    elif cmd == "collect":
        o = collect()
        print({k: (v["state"], v.get("by_variant")) for k, v in o["gates"].items()}, o["n_rows"], "rows")
    elif cmd == "pilot":
        pilot()
