"""Phase 28 (model side): symmetric ONE-diode scalar cell, 9.0 mm period, 10.0 GHz, for RT/duroid 5880 3.175 mm (first choice) and RO4003C 1.524 mm (backup). MODEL_ONLY, IDEALISED, NOT_HFSS_VALIDATED.

    python scripts/research/phase28_scalar_ring_cell.py screen  <out.jsonl>   [--workers N]   # variants x diode cases x 2 states, mesh 45
    python scripts/research/phase28_scalar_ring_cell.py confirm <screen.jsonl> <out.jsonl> [--workers N]   # best per (substrate, diode case): meshes 60/90, +-20 % tolerance, band 9.0 - 11.0 GHz
    python scripts/research/phase28_scalar_ring_cell.py report  <screen.jsonl> <confirm.jsonl> <out.json>

Cell: a closed square metal LOOP (outer side So, strip width wr, D4-symmetric), optionally a centre square patch (side sc); the ONE diode sits in the middle of the top side (current along x, region g x wr, `diode_load`).
The loop is geometrically closed (as in the earlier necks) and the diode load is the series element: OFF = (nearly) open, ON = a few ohm. Dimensions are multiples of 0.6 mm (self-similar staircase on meshes 45/60/90).
Diode cases (from the start, part of the resonance design): MADP-000907-14020 with L_s = 0 and 0.2 nH (ON 5.2 ohm, OFF 0.025 pF; the L_s values are ASSUMED); packaged SC-79 SMP1340-type: R = 2 ohm (read in [78]), L = 0.7 nH
(ASSUMED bound), C_T = 0.19 pF (PROXY: SMP1345-079LF value, the SMP1340 C_T is NOT_AVAILABLE_FROM_CURRENT_DATA).
Leakage: (a) first-order closed form from the two scalar sheets (x-x element), (b) EXACT in the modulation depth within the local sheet model: `floquet_complex_beta` of the periodic two-level sheet
(6 cells per period, duty 0.5, the 32 deg radiating harmonic), using the complex x-x sheet admittances (conductance included). Required: alpha 1.5056 Np/m, DeltaB 4.3 mS, mean b_bar 0.93 mS (RT5880 3.175).
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
import warnings  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from dataclasses import dataclass  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

F0 = 10.0e9
PITCH = 9.0e-3
SUBSTRATES = {"rt5880_3p175": (2.2, 0.0009, 3.175e-3), "ro4003c_1p524": (3.38, 0.0027, 1.524e-3)}
ALPHA_REQ, DB_REQ_MS = 1.5056, 4.3
VARIANTS = {  # So, wr, sc [mm]
    "R1_So6.6_w0.6": (6.6, 0.6, 0.0),
    "R2_So7.8_w0.6": (7.8, 0.6, 0.0),
    "R3_So5.4_w0.6": (5.4, 0.6, 0.0),
    "R4_So6.6_w1.2": (6.6, 1.2, 0.0),
    "R5_So7.8_w1.2": (7.8, 1.2, 0.0),
    "R6_So6.6_w0.6_c1.8": (6.6, 0.6, 1.8),
}
VARIANTS_BY_SUB = {"rt5880_3p175": list(VARIANTS), "ro4003c_1p524": ["R1_So6.6_w0.6", "R2_So7.8_w0.6", "R4_So6.6_w1.2", "R6_So6.6_w0.6_c1.8"]}
CASES = ["MADP_Ls0", "MADP_Ls0.2nH", "SC79_R2_L0.7nH_C0.19pF"]
GAP = 0.6e-3


@dataclass(frozen=True)
class LoopGeometry:
    outer: float
    width: float
    centre: float
    kind = "loop_patch"
    polarity = "metal"
    rotation = 0.0

    @property
    def feature_area(self):
        return self.outer**2

    def corners(self):
        h = self.outer / 2
        return np.array([[h, h], [-h, h], [-h, -h], [h, -h]])

    def in_feature(self, x, y):
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        h = self.outer / 2 + 1e-12
        box = (np.abs(x) <= h) & (np.abs(y) <= h)
        inner = (np.abs(x) < self.outer / 2 - self.width - 1e-12) & (np.abs(y) < self.outer / 2 - self.width - 1e-12)
        ring = box & ~inner
        patch = (np.abs(x) <= self.centre / 2 + 1e-12) & (np.abs(y) <= self.centre / 2 + 1e-12) if self.centre > 0 else np.zeros(x.shape, bool)
        return ring | patch


def diode(case, on, scale=None):
    from src.ris.diodes import DiodeModel

    sr, sl, sc = (scale or {}).get("R", 1.0), (scale or {}).get("L", 1.0), (scale or {}).get("C", 1.0)
    if case.startswith("MADP"):
        ls = (0.2e-9 if "0.2" in case else 0.0) * sl
        return DiodeModel("MADP ON", R_s=5.2 * sr, L_s=ls) if on else DiodeModel("MADP OFF", L_s=ls, C_T=0.025e-12 * sc)
    return DiodeModel("SC79 ON", R_s=2.0 * sr, L_s=0.7e-9 * sl) if on else DiodeModel("SC79 OFF", L_s=0.7e-9 * sl, C_T=0.19e-12 * sc)


def build(variant, substrate, case, on, scale=None):
    from src.geometry import UnitCell
    from src.materials import Metal, Substrate
    from src.periodic_bc import RectangularLattice
    from src.ris.diodes import diode_load

    so, wr, sc = (v * 1e-3 for v in VARIANTS[variant])
    er, td, h = SUBSTRATES[substrate]
    sub = Substrate(substrate, er, td, h)
    load = diode_load(diode(case, on, scale), 0.0, so / 2 - wr / 2, GAP, wr, direction="x", transverse_conducts=False)
    return UnitCell(RectangularLattice(PITCH, PITCH), LoopGeometry(so, wr, sc), sub, Metal(), backing="pec", loads=(load,)), sub


def solve(job):
    variant, substrate, case, on, mesh, f_hz, scale = job
    warnings.simplefilter("ignore")
    from src.hologram.coded_cell_eval import sheet_tensor

    uc, _ = build(variant, substrate, case, on, scale)
    y = sheet_tensor(uc, f_hz, mesh)
    return {"variant": variant, "substrate": substrate, "case": case, "on": on, "mesh": mesh, "f_hz": f_hz, "scale": scale,
            "y_s": [[[float(v.real), float(v.imag)] for v in row] for row in y]}


def key(r):
    return json.dumps([r["variant"], r["substrate"], r["case"], r["on"], r["mesh"], r["f_hz"], r["scale"]], sort_keys=True)


def run(jobs, out, workers):
    out = Path(out)
    done = set()
    if out.exists():
        done = {key(json.loads(line)) for line in out.read_text(encoding="utf-8").splitlines()}
    todo = [j for j in jobs if key(dict(zip(("variant", "substrate", "case", "on", "mesh", "f_hz", "scale"), j, strict=True))) not in done]
    print(len(todo), "solves", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex, out.open("a", encoding="utf-8") as fh:
        for r in ex.map(solve, todo):
            fh.write(json.dumps(r) + "\n")
            fh.flush()


def rows_of(path):
    p = Path(path)
    return [json.loads(line) for line in p.read_text(encoding="utf-8").splitlines()] if p.exists() else []


def pair_metrics(r_on, r_off):
    from src.dispersion.leakage_model import floquet_complex_beta, leakage_alpha_first_order, state_surface_reactances, two_level_fourier
    from src.dispersion.sheet_model import find_beta_from_sheet_admittance
    from src.frequency import free_space_wavenumber
    from src.hologram import coded_feasibility as cf
    from src.materials import Substrate

    f = r_on["f_hz"]
    er, td, h = SUBSTRATES[r_on["substrate"]]
    sub = Substrate("s", er, td, h)
    cf.set_frequency(f)
    cf.SUB = sub
    k0 = free_space_wavenumber(f)
    on = np.array([[complex(*c) for c in row] for row in r_on["y_s"]])
    off = np.array([[complex(*c) for c in row] for row in r_off["y_s"]])
    b_on, b_off = on[0, 0].imag, off[0, 0].imag
    db = abs(b_on - b_off)
    out = {"b_on_ms": b_on * 1e3, "b_off_ms": b_off * 1e3, "b_bar_ms": 0.5 * (b_on + b_off) * 1e3, "d_b_ms": db * 1e3, "g_on_ms": on[0, 0].real * 1e3, "g_off_ms": off[0, 0].real * 1e3,
           "aniso_on_ms": (on[0, 0].imag - on[1, 1].imag) * 1e3, "aniso_off_ms": (off[0, 0].imag - off[1, 1].imag) * 1e3, "xy_on_ms": abs(on[0, 1]) * 1e3, "xy_off_ms": abs(off[0, 1]) * 1e3}
    ybar = 0.5 * (on[0, 0] + off[0, 0])
    try:
        beta = find_beta_from_sheet_admittance(1j * ybar.imag, sub, f).beta_rad_per_m
        out["beta_over_k0"] = beta / k0
        xs = state_surface_reactances([1j * b_on, 1j * b_off], beta, sub, f)
        out["dx_switch_ohm"] = float(abs(xs[0] - xs[1]))
        g1 = cf.first_harmonic_gain(cf.uniform_levels(1))
        ks = min(leakage_alpha_first_order(beta, 1j * ybar.imag, 1.0, math.radians(30.0), math.radians(z), sub, f) for z in (0.0, 90.0))
        out["alpha_first_order"] = float(ks * (db * g1 / 4.0) ** 2)
        out["modulation_depth"] = float((db / 2) / abs(1j * ybar.imag + cf._slab_tm(beta)))
        yc = two_level_fourier(off[0, 0], on[0, 0], 6, duty=0.5)
        b, _ = floquet_complex_beta(yc, 6 * PITCH, 0.0, sub, f, beta, n_max=6)
        out["alpha_exact"] = float(-b.imag)
        out["beta_real_over_k0_exact"] = float(b.real / k0)
        out["beam_angle_deg_exact"] = float(math.degrees(math.asin(min(abs(b.real - 2 * math.pi / (6 * PITCH)) / k0, 1.0))))
    except Exception as exc:  # noqa: BLE001
        out["note"] = f"{type(exc).__name__}: {exc}"
    return out


def table(rows):
    g = {}
    for r in rows:
        g.setdefault((r["variant"], r["substrate"], r["case"], r["mesh"], r["f_hz"], json.dumps(r["scale"], sort_keys=True)), {})[r["on"]] = r
    return {k: pair_metrics(v[True], v[False]) for k, v in g.items() if True in v and False in v}


def score(m):
    if not m.get("alpha_exact") or m["alpha_exact"] <= 0 or not m.get("d_b_ms"):
        return 99.0
    return abs(math.log(m["alpha_exact"] / ALPHA_REQ)) + 0.5 * abs(math.log(m["d_b_ms"] / DB_REQ_MS)) + (0.0 if 0.95 < m.get("beta_over_k0", 0) < 1.25 else 1.0)


def best(screen):
    t = table(screen)
    b = {}
    for (variant, sub, case, mesh, f, sc), m in t.items():
        if mesh != 45 or f != F0:
            continue
        s = score(m)
        if (sub, case) not in b or s < b[(sub, case)][1]:
            b[(sub, case)] = (variant, s)
    return b


def main(argv):
    cmd, w = argv[0], (int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 40)
    if cmd == "screen":
        jobs = [(v, s, c, on, 45, F0, None) for s in SUBSTRATES for v in VARIANTS_BY_SUB[s] for c in CASES for on in (True, False)]
        run(jobs, argv[1], w)
    elif cmd == "confirm":
        jobs = []
        for (sub, case), (variant, _) in best(rows_of(argv[1])).items():
            if sub != "rt5880_3p175":
                continue
            for on in (True, False):
                for mesh in (60, 90):
                    jobs.append((variant, sub, case, on, mesh, F0, None))
                for sc in ({"R": 0.8, "L": 0.8, "C": 0.8}, {"R": 1.2, "L": 1.2, "C": 1.2}, {"C": 0.8}, {"C": 1.2}):
                    jobs.append((variant, sub, case, on, 60, F0, sc))
                for f in (9.0, 9.25, 9.5, 9.75, 10.25, 10.5, 10.75, 11.0):
                    jobs.append((variant, sub, case, on, 60, f * 1e9, None))
        run(jobs, argv[2], w)
    elif cmd == "report":
        screen, conf = rows_of(argv[1]), rows_of(argv[2])
        out = {"labels": ["MODEL_ONLY", "IDEALISED_CELLS", "NOT_HFSS_VALIDATED"], "required": {"alpha_np_per_m": ALPHA_REQ, "d_b_ms": DB_REQ_MS, "b_bar_ms": 0.93},
               "screen": {"|".join(map(str, k)): v for k, v in table(screen).items()}, "confirm": {"|".join(map(str, k)): v for k, v in table(conf).items()},
               "best": {"|".join(k): v[0] for k, v in best(screen).items()}}
        Path(argv[3]).write_text(json.dumps(out, indent=1), encoding="utf-8")
        print("wrote", argv[3])


if __name__ == "__main__":
    main(sys.argv[1:])
