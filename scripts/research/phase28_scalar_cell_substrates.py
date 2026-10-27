"""Phase 28 (model side): required switching contrast of a 1-bit SCALAR cell for the sector-fed design, and what a simple one-diode split-patch cell reaches on three substrates.

    python scripts/research/phase28_scalar_cell_substrates.py screen  <out.jsonl> [--workers N]    # 4 variants x 3 substrates x 3 diode models x 2 states, mesh 45
    python scripts/research/phase28_scalar_cell_substrates.py confirm <screen.jsonl> <out.jsonl> [--workers N]   # best variant per substrate, meshes 60 and 90
    python scripts/research/phase28_scalar_cell_substrates.py report  <screen.jsonl> <confirm.jsonl> <out.json>

Cell (a CHOICE for this screen, in the spirit of [78], not its geometry): period 9.0 mm, a square metal patch of side S split by a gap g across x (|x| < g/2) except a bridging strip of width w at the centre
that carries the ONE diode (lumped load, current along x, the TM0 field direction). State 1 = diode ON, state 0 = OFF. The sheet tensor is the normal-incidence one (`coded_cell_eval`); the x-x element is the
scalar sheet seen by a surface wave travelling along x. Dimensions are multiples of 0.6 mm (self-similar staircase on meshes 45, 60, 90).
Diode models: MADP-000907-14020 (repo: ON 5.2 ohm, OFF 0.025 pF, L_s = 0); PKG: packaged SMP1340-type: ON R = 2 ohm (read in [78]), OFF C_T = 0.19 pF (ASSUMED: the SMP1345-079LF value, the SMP1340 C_T
is NOT in the repo or in [78]'s text), package inductance L_s = 0.3 or 0.7 nH in both states (ASSUMED bounds). Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET.
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

F_HZ = 10.0e9
PITCH = 9.0e-3
SUBSTRATES = {
    "ro3006_6p15_1p27": (6.15, 0.0022, 1.27e-3),
    "ro4003c_3p38_1p524": (3.38, 0.0027, 1.524e-3),
    "rt5880_2p2_3p175": (2.2, 0.0009, 3.175e-3),
}
VARIANTS = {  # side S, gap g, strap width w [mm]
    "V1_S8.4_g0.6_w0.6": (8.4, 0.6, 0.6),
    "V2_S8.4_g1.2_w0.6": (8.4, 1.2, 0.6),
    "V3_S7.2_g0.6_w0.6": (7.2, 0.6, 0.6),
    "V4_S8.4_g0.6_w1.2": (8.4, 0.6, 1.2),
}
MODELS = ["MADP_Ls0", "PKG_Ls0.3", "PKG_Ls0.7"]
ALPHA_REQUIRED = 1.5056  # Np/m: sector-fed design, converted fraction 0.504 (results/phase28/f10/feeds/sector_feed.json)


@dataclass(frozen=True)
class SplitPatchGeometry:
    side: float
    gap: float
    strap: float
    kind = "split_patch"
    polarity = "metal"
    rotation = 0.0

    @property
    def feature_area(self):
        return self.side**2

    def corners(self):
        h = self.side / 2
        return np.array([[h, h], [-h, h], [-h, -h], [h, -h]])

    def in_feature(self, x, y):
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        h = self.side / 2 + 1e-12
        inside = (np.abs(x) <= h) & (np.abs(y) <= h)
        cut = (np.abs(x) < self.gap / 2 - 1e-12) & (np.abs(y) > self.strap / 2 + 1e-12)
        return inside & ~cut


def diode(model, on):
    from src.ris.diodes import DiodeModel

    if model == "MADP_Ls0":
        return DiodeModel("MADP ON", R_s=5.2) if on else DiodeModel("MADP OFF", C_T=0.025e-12)
    ls = float(model.split("Ls")[1]) * 1e-9
    return DiodeModel("PKG ON", R_s=2.0, L_s=ls) if on else DiodeModel("PKG OFF", L_s=ls, C_T=0.19e-12)


def unit_cell(variant, substrate, model, on):
    from src.geometry import UnitCell
    from src.materials import Metal, Substrate
    from src.periodic_bc import RectangularLattice
    from src.ris.diodes import diode_load

    s, g, w = (v * 1e-3 for v in VARIANTS[variant])
    er, td, h = SUBSTRATES[substrate]
    sub = Substrate(substrate, er, td, h)
    load = diode_load(diode(model, on), 0.0, 0.0, g, w, direction="x", transverse_conducts=False)
    return UnitCell(RectangularLattice(PITCH, PITCH), SplitPatchGeometry(s, g, w), sub, Metal(), backing="pec", loads=(load,)), sub


def solve(job):
    variant, substrate, model, on, mesh = job
    warnings.simplefilter("ignore")
    from src.hologram.coded_cell_eval import sheet_tensor

    uc, sub = unit_cell(variant, substrate, model, on)
    y = sheet_tensor(uc, F_HZ, mesh)
    return {
        "variant": variant,
        "substrate": substrate,
        "model": model,
        "on": on,
        "mesh": mesh,
        "y_s": [[[float(v.real), float(v.imag)] for v in row] for row in y],
    }


def _run(jobs, out, workers):
    out = Path(out)
    done = set()
    if out.exists():
        for line in out.read_text(encoding="utf-8").splitlines():
            r = json.loads(line)
            done.add((r["variant"], r["substrate"], r["model"], r["on"], r["mesh"]))
    todo = [j for j in jobs if j not in done]
    print(len(todo), "solves", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex, out.open("a", encoding="utf-8") as fh:
        for r in ex.map(solve, todo):
            fh.write(json.dumps(r) + "\n")
            fh.flush()


def required_contrast():
    """Per substrate: bare-slab beta/k0, the target mean loading, and the 1-bit scalar contrast the design needs for ALPHA_REQUIRED (30 deg beam)."""
    from scipy.optimize import brentq

    from src.dispersion.sheet_model import find_beta_from_sheet_admittance
    from src.frequency import free_space_wavenumber
    from src.hologram import coded_feasibility as cf
    from src.materials import Substrate

    cf.set_frequency(F_HZ)
    k0 = free_space_wavenumber(F_HZ)
    out = {}
    for name, (er, td, h) in SUBSTRATES.items():
        sub = Substrate(name, er, td, h)
        cf.SUB = sub
        bare = find_beta_from_sheet_admittance(1e-9j, sub, F_HZ).beta_rad_per_m / k0
        target = max(1.053, bare + 0.015)
        row = {"eps_r": er, "thickness_mm": h * 1e3, "beta_over_k0_bare_slab": bare, "beta_over_k0_target": target}
        try:
            f = lambda b: find_beta_from_sheet_admittance(1j * b * 1e-3, sub, F_HZ).beta_rad_per_m / k0 - target  # noqa: E731
            row["b_bar_for_target_ms"] = brentq(f, 0.01, 40.0, xtol=1e-9)
            c = cf.contrast_from_alpha(ALPHA_REQUIRED, target, 30.0, 1)["scalar"]
            row.update({"required_dx_switch_ohm": c["dx_switch_ohm"], "required_d_b_full_swing_ms": c["d_b_full_swing_ms"]})
        except Exception as exc:  # noqa: BLE001
            row["note"] = f"{type(exc).__name__}: {exc}"
        out[name] = row
    return out


def metrics(r_on, r_off, substrate):
    """Scalar (x-x) switching metrics of the two states."""
    from src.dispersion.leakage_model import leakage_alpha_first_order, state_surface_reactances
    from src.dispersion.sheet_model import find_beta_from_sheet_admittance
    from src.frequency import free_space_wavenumber
    from src.hologram import coded_feasibility as cf
    from src.materials import Substrate

    er, td, h = SUBSTRATES[substrate]
    sub = Substrate(substrate, er, td, h)
    cf.SUB = sub
    k0 = free_space_wavenumber(F_HZ)
    yo = np.array([[complex(*c) for c in row] for row in r_on["y_s"]])
    yf = np.array([[complex(*c) for c in row] for row in r_off["y_s"]])
    ybar = 0.5 * (yo[0, 0] + yf[0, 0])
    out = {"b_on_ms": yo[0, 0].imag * 1e3, "b_off_ms": yf[0, 0].imag * 1e3, "g_on_ms": yo[0, 0].real * 1e3, "g_off_ms": yf[0, 0].real * 1e3,
           "aniso_on_ms": (yo[0, 0].imag - yo[1, 1].imag) * 1e3, "aniso_off_ms": (yf[0, 0].imag - yf[1, 1].imag) * 1e3}
    db = abs(yo[0, 0].imag - yf[0, 0].imag)
    out["d_b_ms"] = db * 1e3
    try:
        beta = find_beta_from_sheet_admittance(1j * ybar.imag, sub, F_HZ).beta_rad_per_m
        xs = state_surface_reactances([1j * yo[0, 0].imag, 1j * yf[0, 0].imag], beta, sub, F_HZ)
        th = math.radians(30.0)
        g1 = cf.first_harmonic_gain(cf.uniform_levels(1))
        c1 = db * g1 / 4.0
        k_s = min(leakage_alpha_first_order(beta, 1j * ybar.imag, 1.0, th, math.radians(z), sub, F_HZ) for z in (0.0, 90.0))
        depth = (db / 2) / abs(1j * ybar.imag + cf._slab_tm(beta))
        out.update({"modulation_depth": float(depth), "first_order_valid": bool(depth <= 0.5 and 1.0 < beta / k0 < 1.2), "beta_over_k0": beta / k0, "dx_switch_ohm": float(abs(xs[0] - xs[1])), "alpha_np_per_m": float(k_s * c1**2), "loss_ratio_g_on_over_db": float(yo[0, 0].real / max(db, 1e-12))})
    except Exception as exc:  # noqa: BLE001
        out.update({"beta_over_k0": None, "note": f"NO_BOUND_MODE: {type(exc).__name__}"})
    return out


def load_rows(path):
    return [json.loads(line) for line in Path(path).read_text(encoding="utf-8").splitlines()]


def table(rows):
    res = {}
    for r in rows:
        res.setdefault((r["variant"], r["substrate"], r["model"], r["mesh"]), {})[r["on"]] = r
    return {k: metrics(v[True], v[False], k[1]) for k, v in res.items() if True in v and False in v}


def best_per_substrate(screen_rows):
    """Per substrate, the MADP variant whose scalar leakage is closest (log scale) to the required one, first-order-valid cases preferred (penalty 1 otherwise)."""
    t = table(screen_rows)
    best = {}
    for (variant, sub, model, mesh), m in t.items():
        if model != "MADP_Ls0" or not m.get("alpha_np_per_m") or m["alpha_np_per_m"] <= 0:
            continue
        score = abs(math.log(m["alpha_np_per_m"] / ALPHA_REQUIRED)) + (0.0 if m.get("first_order_valid") else 1.0)
        if sub not in best or score < best[sub][1]:
            best[sub] = (variant, score)
    return {s: v[0] for s, v in best.items()}


def main(argv):
    cmd = argv[0]
    w = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 40
    if cmd == "screen":
        jobs = [(v, s, m, on, 45) for v in VARIANTS for s in SUBSTRATES for m in MODELS for on in (True, False)]
        _run(jobs, argv[1], w)
    elif cmd == "confirm":
        best = best_per_substrate(load_rows(argv[1]))
        jobs = [(v, s, m, on, mesh) for s, v in best.items() for m in ("MADP_Ls0", "PKG_Ls0.3") for on in (True, False) for mesh in (60, 90)]
        _run(jobs, argv[2], w)
    elif cmd == "report":
        screen, confirm = load_rows(argv[1]), load_rows(argv[2]) if Path(argv[2]).exists() else []
        out = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "PROVISIONAL_DATASHEET"], "alpha_required_np_per_m": ALPHA_REQUIRED, "required": required_contrast(),
               "screen": {"|".join(map(str, k)): v for k, v in table(screen).items()}, "confirm": {"|".join(map(str, k)): v for k, v in table(confirm).items()},
               "best_variant_per_substrate": best_per_substrate(screen)}
        Path(argv[3]).write_text(json.dumps(out, indent=1), encoding="utf-8")
        print("wrote", argv[3])


if __name__ == "__main__":
    main(sys.argv[1:])
