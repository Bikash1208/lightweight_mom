"""Phase 26 S: substrate study -- reachable conversion on thicker / higher-eps boards (model only).

    python scripts/research/phase26_substrate_study.py prelim <out.json>
    python scripts/research/phase26_substrate_study.py table <name> <out.json> --slot L,W [--cells 3,4,5,6,8] [--mesh 24]
    python scripts/research/phase26_substrate_study.py reach <out.json>

`prelim`: per candidate board the TE1 cutoff (grounded slab: c / (4 h sqrt(eps - 1)); excluded if <= 10.54 GHz),
bare TM0 beta, the pitch/diameter/slot family re-sized for the board and its smallest features, and the
pin-feed launch fraction (`pin_launch`, bare slab; pin radius 0.305 mm; pin top flush with the surface and
1.22 mm above it = the Phase 26 A reading).
Family sizing rule (the only one tried -- stated, not optimised): pitch P = 2 pi / (3 K_max) with K_max =
(1.04 beta_bare + k0 sin 30 deg), rounded to 0.1 mm, so the shortest modulation period of the hexagon is three
cells; circle D = 0.9375 P (6.0 / 6.4 of the 5880 family); the four slots keep the 5880 family's L/D and W/D.
`table`: the SOLVER_SUPERCELL table of ONE slot for the board (as `phase26_supercell_table`), `reach`:
best-case conversion over the hexagon from all tables of all boards (as `phase26_table_reach`).
Labels: SOLVER_SUPERCELL / MODEL_ONLY / NOT_HFSS_VALIDATED / GAMMA_0_COLLINEAR_ONLY.
"""

import glob
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.basis_functions import MeshSettings, RooftopBasis  # noqa: E402
from src.dispersion import tensor_sheet as ts  # noqa: E402
from src.dispersion.analytic_tm0 import tm0_beta  # noqa: E402
from src.dispersion.pin_launch import Launcher, PinLaunch  # noqa: E402
from src.dispersion.sheet_model import _y_air, _y_slab  # noqa: E402
from src.dispersion.slotted_cells import CellRow, SlottedCircularPatch  # noqa: E402
from src.frequency import angular_frequency, free_space_wavenumber  # noqa: E402
from src.geometry import UnitCell  # noqa: E402
from src.hologram.leaky_forward import polygon_ray_distance  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402

F = 10.54e9
K0 = free_space_wavenumber(F)
C0 = 299792458.0
OUT = ROOT / "results/phase26/substrate"
# name: (eps_r, tan_delta, h_mm). tan_delta: manufacturer values at 10 GHz (RT5880 0.0009, RO3006 0.0020,
# RO3010 0.0022, RT6010 0.0023) -- datasheet nominals, not measured.
BOARDS = {
    "rt5880_1p575": (2.2, 0.0009, 1.575),
    "rt5880_3p175": (2.2, 0.0009, 3.175),
    "ro3006_1p28": (6.15, 0.0020, 1.28),
    "ro3006_1p52": (6.15, 0.0020, 1.52),
    "ro3010_1p27": (10.2, 0.0022, 1.27),
    "rt6010_1p90": (10.2, 0.0023, 1.90),
}
# slot (L/D, W/D) of the 5880 family at D = 6.0: 4.8x2.0, 4.8x1.6, 4.6x1.6, 4.4x1.2
SLOT_REL = {"deep": (0.8, 2.0 / 6.0), "mid": (0.8, 1.6 / 6.0), "mid2": (4.6 / 6.0, 1.6 / 6.0), "weak": (4.4 / 6.0, 1.2 / 6.0)}


def te1_cutoff_ghz(eps, h_m):
    return C0 / (4 * h_m * np.sqrt(eps - 1)) / 1e9


def target_beta(beta_bare):
    """Design loaded beta of the uniform (slot-0, deep) cell: bare + 0.012 k0, as the 5880 1.575 mm control
    (1.0305 vs bare 1.0185), never below 1.03 k0."""
    return max(1.03 * K0, beta_bare + 0.012 * K0)


def family(name):
    eps, tand, h = BOARDS[name]
    beta = tm0_beta(F, eps, h * 1e-3)
    sizing = OUT / "sizing.json"
    if sizing.exists() and name in json.loads(sizing.read_text()):
        s = json.loads(sizing.read_text())[name]
        pitch, d = s["pitch_mm"], s["diameter_mm"]
    else:
        k_max = 1.04 * beta + K0 * 0.5
        pitch = round(2 * np.pi / (3 * k_max) * 1e3, 1)
        d = round(0.9375 * pitch * 20) / 20
    slots = {}
    for k, (lr, wr) in SLOT_REL.items():
        L, W = round(lr * d * 20) / 20, round(wr * d * 20) / 20
        chord_half = np.sqrt((d / 2) ** 2 - (W / 2) ** 2)
        slots[k] = {"L_mm": L, "W_mm": W, "bridge_mm": float(chord_half - L / 2)}
    return eps, tand, h, beta, pitch, d, slots


def size(out):
    """Pitch P = 2 pi / (3 (beta_target + k0/2)); diameter by a scan of D/P in 0.025 steps: the largest D
    (smallest gap >= 0.35 mm) whose uniform deep-slot cell has beta <= beta_target + 0.06 k0. A failed Bloch solve (no bound root) is recorded, never skipped silently."""
    res = {}
    for name, (eps, tand, h) in BOARDS.items():
        bare = tm0_beta(F, eps, h * 1e-3)
        bt = target_beta(bare)
        pitch = round(2 * np.pi / (3 * (bt + 0.5 * K0)) * 1e3, 1)
        sub = Substrate(name, eps, tand, h * 1e-3)
        rows = []
        for ratio in np.arange(0.25, 0.951, 0.025):
            d = round(ratio * pitch * 40) / 40
            if pitch - d < 0.35 - 1e-9:
                continue
            L, W = round(0.8 * d * 20) / 20, round(d / 3.0 * 20) / 20
            geo = SlottedCircularPatch(d * 1e-3, L * 1e-3, W * 1e-3, 0.0)
            cell = UnitCell(RectangularLattice(pitch * 1e-3, pitch * 1e-3), geo, sub, Metal(), backing="pec")
            try:
                b, _ = ts.bloch_beta(cell, F, MeshSettings(24, 24, 24), 0.0, n_grid=5, refine=1)
                rows.append({"D_mm": d, "beta_over_k0": b / K0})
            except Exception as exc:  # noqa: BLE001
                rows.append({"D_mm": d, "error": type(exc).__name__})
        ok = [r for r in rows if "beta_over_k0" in r]
        pick = None
        if ok:
            # the 5880 1.575 control sits just below the resonance edge (D 6.0 on 6.4 mm; 6.2 has no bound
            # root): pick the LARGEST rooted D with beta <= target + 0.06 k0
            cand = [r for r in ok if bare - 1e-3 * K0 <= r["beta_over_k0"] * K0 <= bt + 0.06 * K0]  # a root below the bare slab's beta is spurious
            pick = max(cand, key=lambda r: r["D_mm"]) if cand else None
        res[name] = {"pitch_mm": pitch, "diameter_mm": pick["D_mm"] if pick else None, "target_beta_over_k0": bt / K0,
                     "bare_beta_over_k0": bare / K0, "scan": rows}
        print(name, pitch, pick, flush=True)
        Path(out).write_text(json.dumps(res, indent=1))


def prelim(out):
    res = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"], "frequency_ghz": F / 1e9, "boards": {}}
    for name, (eps, tand, h) in BOARDS.items():
        e, t, hh, beta, pitch, d, slots = family(name)
        rec = {"eps_r": e, "tan_delta": t, "h_mm": hh, "te1_cutoff_ghz": float(te1_cutoff_ghz(e, hh * 1e-3))}
        rec["te1_propagates_at_10p54"] = bool(rec["te1_cutoff_ghz"] <= F / 1e9)
        rec["beta_bare_over_k0"] = beta / K0
        rec["pitch_mm"], rec["diameter_mm"], rec["gap_mm"] = pitch, d, pitch - d
        rec["slots"] = slots
        rec["smallest_feature_mm"] = float(min([pitch - d] + [s["W_mm"] for s in slots.values()] + [s["bridge_mm"] for s in slots.values()]))
        rec["launch"] = {}
        for tag, top in (("pin_top_flush", 0.0), ("pin_top_1p22mm_above", 1.22e-3)):
            pl = PinLaunch(F, e, hh * 1e-3, Launcher(top, 0.305e-3, 0.0), tan_delta=t, n_seg=24)
            s = pl.split()
            rec["launch"][tag] = {k: float(v) for k, v in s.items()}
        res["boards"][name] = rec
        print(name, json.dumps({k: rec[k] for k in ("te1_cutoff_ghz", "te1_propagates_at_10p54", "beta_bare_over_k0", "pitch_mm", "diameter_mm", "smallest_feature_mm")}),
              {k: round(v["launch_fraction"], 3) for k, v in rec["launch"].items()}, flush=True)
    Path(out).write_text(json.dumps(res, indent=1))


class Row:
    def __init__(self, geos, sub, pitch, mesh_n):
        n = len(geos)
        self.sub, self.pitch = sub, pitch
        self.uc = UnitCell(RectangularLattice(n * pitch, pitch), CellRow(tuple(geos), pitch), sub, Metal(), backing="pec")
        self.mesh = MeshSettings(n * mesh_n, mesh_n, max(n * mesh_n, mesh_n))
        self.basis = RooftopBasis(self.uc, self.mesh.cells_x, self.mesh.cells_y)
        self.orders = (n * mesh_n, mesh_n)
        self.omega = angular_frequency(F)

    def _y_bg(self, beta):
        s = self.sub
        a = _y_air(beta, self.omega, K0) + _y_slab(beta, self.omega, s.eps_complex, K0, s.thickness)
        b = ts._b_te(beta, self.omega, s.eps_complex, K0, s.thickness)
        return a, b

    def y(self, beta):
        return ts.bloch_sheet_tensor(self.uc, F, self.mesh, beta, 0.0, basis=self.basis, orders=self.orders)

    def f(self, beta):
        a, b = self._y_bg(beta)
        return complex(np.linalg.det(self.y(beta) + np.array([[a, 0.0], [0.0, b]])))

    def f_imag(self, beta):
        a, b = self._y_bg(beta)
        return float(np.linalg.det(self.y(beta).imag + np.array([[a.imag, 0.0], [0.0, b.imag]])))

    def real_root(self, guess, lo, hi):
        """Secant from `guess`; falls back to a scan of `[lo, hi]` for the sign change nearest the guess."""
        b0, b1 = guess, guess * 1.004
        f0, f1 = self.f_imag(b0), self.f_imag(b1)
        for _ in range(10):
            if f1 == f0:
                break
            b2 = b1 - f1 * (b1 - b0) / (f1 - f0)
            b0, f0, b1 = b1, f1, b2
            f1 = self.f_imag(b1)
            if abs(b1 - b0) < 1e-4 * K0:
                break
        if not (lo < b1 < hi):
            raise RuntimeError(f"secant left [{lo / K0:.3f}, {hi / K0:.3f}] k0: {b1 / K0:.3f}")
        return float(b1)

    def complex_root(self, beta_r, step_k0=0.004, m=3, deg=6):
        h = step_k0 * K0
        xs = beta_r + h * np.arange(-m, m + 1)
        vals = np.array([self.f(x) for x in xs])
        p6 = np.polynomial.Polynomial.fit((xs - beta_r) / h, vals, deg)
        p5 = np.polynomial.Polynomial.fit((xs - beta_r) / h, vals, deg - 1)
        r6 = p6.roots()
        r6 = r6[np.argmin(np.abs(r6))]
        r5 = p5.roots()
        r5 = r5[np.argmin(np.abs(r5))]
        return beta_r + h * complex(r6), abs(r6 - r5) * h


def table(name, out, slot, cells, mesh_n):
    eps, tand, h, beta_bare, pitch, d, slots = family(name)
    sub = Substrate(name, eps, tand, h * 1e-3)
    p = pitch * 1e-3
    L, W = slot
    res = {"labels": ["SOLVER_SUPERCELL", "MODEL_ONLY", "NOT_HFSS_VALIDATED", "GAMMA_0_COLLINEAR_ONLY"], "board": name,
           "eps_r": eps, "h_mm": h, "pitch_mm": pitch, "diameter_mm": d, "slot_L_W_mm": [L, W], "mesh_per_cell": mesh_n, "points": []}
    lam0 = C0 / F
    lo, hi = 1.004 * K0, 0.995 * np.sqrt(eps) * K0
    for n in cells:
        t0 = time.time()
        geos = [SlottedCircularPatch(d * 1e-3, L * 1e-3, W * 1e-3, np.radians(180.0 * i / n)) for i in range(n)]
        row = Row(geos, sub, p, mesh_n)
        uni = Row([SlottedCircularPatch(d * 1e-3, L * 1e-3, W * 1e-3, 0.0)] * n, sub, p, mesh_n)
        rec = {"n_cells": n, "period_mm": n * pitch, "K_over_k0": lam0 / (n * p) / 1.0 * (2 * np.pi / lam0) / K0 * 1.0}
        rec["K_over_k0"] = (2 * np.pi / (n * p)) / K0
        try:
            # guess: the uniform row's own root, then the rotating row from it
            bu = uni.real_root(beta_bare * 1.04, lo, hi)
            br = row.real_root(bu * 1.01, lo, hi)
            bc, spread = row.complex_root(br)
            bcu, _ = uni.complex_root(bu)
            rec.update({"beta_real_root_over_k0": br / K0, "beta_over_k0": bc.real / K0, "alpha_np_per_m": float(-bc.imag),
                        "beta_uniform_over_k0": bcu.real / K0, "alpha_uniform_np_per_m": float(-bcu.imag),
                        "alpha_leak_np_per_m": float(-bc.imag + bcu.imag), "root_spread_over_k0": spread / K0})
        except Exception as exc:  # reported, never filled in
            rec["error"] = f"{type(exc).__name__}: {exc}"
        rec["seconds"] = time.time() - t0
        print(rec, flush=True)
        res["points"].append(rec)
        Path(out).write_text(json.dumps(res, indent=1))


def reach(out):
    R = 0.150
    v = np.array([[R * np.cos(np.radians(60 * i)), R * np.sin(np.radians(60 * i))] for i in range(6)])
    psi = np.radians(np.arange(0, 360, 1.0))
    red = polygon_ray_distance(v, (0.0, 0.0), psi)
    pre = json.loads((OUT / "prelim.json").read_text())
    res = {"labels": ["SOLVER_SUPERCELL", "MODEL_ONLY", "GAMMA_0_COLLINEAR_ONLY"], "boards": {}}
    for name, b in pre["boards"].items():
        pts = []
        for f in sorted(glob.glob(str(OUT / f"table_{name}_*.json"))):
            pts += [p for p in json.loads(Path(f).read_text())["points"] if "alpha_leak_np_per_m" in p]
        n_all = len(pts)
        # QC: positive alpha, polynomial roots (deg 6 vs 5) agree to 1e-3 k0, complex root within 0.03 k0 of the
        # real-axis root (a root that jumps is the row's own resonance, not the continued bound mode)
        pts = [p for p in pts if p["alpha_leak_np_per_m"] > 0 and p["root_spread_over_k0"] < 1e-3
               and abs(p["beta_over_k0"] - p["beta_real_root_over_k0"]) < 0.03]
        if not pts:
            res["boards"][name] = {"state": "NO_TABLE_POINTS", "te1_propagates": b["te1_propagates_at_10p54"]}
            continue
        ks = sorted({round(p["K_over_k0"], 4) for p in pts})
        am = [max(p["alpha_leak_np_per_m"] for p in pts if round(p["K_over_k0"], 4) == k) for k in ks]
        beta = float(np.mean([p["beta_over_k0"] for p in pts]))
        kk = np.sqrt(beta**2 + 0.25 - beta * np.cos(psi - np.radians(225)))
        covered = (kk.min() >= ks[0] - 0.05) and (kk.max() <= ks[-1] + 0.05)
        a = np.interp(kk, ks, am)
        eta = 1 - np.exp(-2 * a * red)
        lf = b["launch"]["pin_top_flush"]["launch_fraction"]
        res["boards"][name] = {
            "n_points": len(pts), "n_points_excluded_by_qc": n_all - len(pts), "K_grid_over_k0": ks, "alpha_max_np_per_m": am, "beta_mean_over_k0": beta,
            "K_range_needed_over_k0": [float(kk.min()), float(kk.max())], "K_range_covered_by_table": bool(covered),
            "alpha_available_mean": float(a.mean()), "alpha_needed_85pct_mean": float((-np.log(0.15) / (2 * red)).mean()),
            "conversion_best_case_mean": float(eta.mean()), "conversion_range": [float(eta.min()), float(eta.max())],
            "gain_estimate_dbi": float(28.5 + 10 * np.log10(lf * eta.mean())),
            "mesh_stable_features_ge_0p35mm": bool(b["smallest_feature_mm"] >= 0.35),
            "launch_fraction_flush": lf,
            "launch_fraction_1p22mm": b["launch"]["pin_top_1p22mm_above"]["launch_fraction"],
            "te1_propagates": b["te1_propagates_at_10p54"], "smallest_feature_mm": b["smallest_feature_mm"],
        }
        print(name, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in res["boards"][name].items() if k in ("conversion_best_case_mean", "launch_fraction_flush", "alpha_available_mean", "K_range_covered_by_table")}, flush=True)
    Path(out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    a = sys.argv[1:]
    OUT.mkdir(parents=True, exist_ok=True)
    if a[0] == "size":
        size(a[1])
    elif a[0] == "prelim":
        prelim(a[1])
    elif a[0] == "table":
        slot = tuple(float(v) for v in a[a.index("--slot") + 1].split(","))
        cells = [int(v) for v in a[a.index("--cells") + 1].split(",")] if "--cells" in a else [3, 4, 5, 6, 8]
        mesh = int(a[a.index("--mesh") + 1]) if "--mesh" in a else 24
        table(a[1], a[2], slot, cells, mesh)
    elif a[0] == "reach":
        reach(a[1])
