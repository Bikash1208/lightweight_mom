"""Phase 20 R3b: CP-2 re-optimisation at 10.5 GHz with separate axis / diagonal stub parameters. OUTSIDE the solver core (Nelder-Mead). Checkpointed: every evaluation is appended to a
JSON-lines cache (fsync) and replayed on restart, so an interrupted run resumes exactly where it stopped.
Labels: MODEL_ONLY with the user-approved HFSS_CALIBRATED delta of the BASELINE cell transferred (HFSS_CALIBRATED_TRANSFERRED), PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.

    python scripts/research/phase20_r3b_opt.py eval <x...>                 # one evaluation (prints the breakdown and timing)
    python scripts/research/phase20_r3b_opt.py run <start_name> [maxfev]   # Nelder-Mead from a named start
    python scripts/research/phase20_r3b_opt.py many <maxfev> <workers> <start> ...   # several starts in threads sharing one process pool (server)

x = circumradius, axis line, axis width, diagonal line, diagonal width [mm]; period 15 mm, MADP-000907-14020, RT5880 1.575 mm.

Objective J = sum_c w_c * v_c(+ small step / deviation terms) over the conditions below (v = normalised violation of the R3 targets: ladder deviation <= 15, step error <= 20, min |r_RR| >= -1.5 dB,
spread <= 1 dB, opposite-handed <= -15 dB), plus a lower-Q term (phase slope of the four r_RR across 10.25-10.75 GHz) and band-edge terms. Conditions (user decision 2026-10-23):
  nominal 10.5 GHz at (theta, phi) = (0,0) (0,0)w3, (30,0)w2, (30,45)w1, (30,90)w1;  corners at 10.5 GHz at (0,0) and (30,45), w1 each: C_T +20 % (datasheet max), C_T -10 % (ASSUMED; the datasheet gives no minimum),
  L_s +0.1 nH (assumed; MADP L_s NOT_AVAILABLE, nominal 0);  band edges 10.25 / 10.75 GHz nominal at (0,0), (30,0).
Corrections added to the mesh-60 / F = 2n MoM: baseline mesh correction (c_corr - c_raw of the F = 2n library at the same point) and, at phi = 0 only, the HFSS delta table of the baseline cell.
"""

import contextlib
import json
import math
import os
import sys
import threading
import time
import warnings
from dataclasses import replace
from multiprocessing import Pool
from pathlib import Path

for _k in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from phase20_r4_metrics import point_metrics  # noqa: E402

from src.basis_functions import MeshSettings  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.ris.cells import Cp2Cell, DiodeKit, madp_kit  # noqa: E402
from src.ris.hfss_calibration import DeltaTable  # noqa: E402
from src.ris.polarization import circular_reflection, to_structure_frame  # noqa: E402
from src.ris.state_library import StateLibrary  # noqa: E402
from src.sparameters import compute_s_parameters_oblique_over_loads  # noqa: E402

PERIOD = 15e-3
MESH = 60
PL = ReferencePlanes(0.0, None, 1.575e-3)
OUT = ROOT / "results" / "phase20_r3b"
CACHE = OUT / "opt_cache.jsonl"
BOUNDS = [
    (3.0, 5.5),
    (0.5, 1.9),
    (0.6, 2.2),
    (0.5, 6.0),
    (0.6, 3.0),
]  # R, axis line (edge-limited), axis width, diag line, diag width [mm]
BASELINE = [4.601469, 1.748, 1.374926, 1.748, 1.374926]
NAMES = ["circumradius", "line_axis", "width_axis", "line_diag", "width_diag"]
TARGETS = {"dev": 15.0, "step": 20.0, "amp_min": -1.5, "spread": 1.0, "opp": -15.0}

# conditions: (name, f_ghz, theta, phi, corner, weight)
NOMINAL = [
    ("n00", 10.5, 0.0, 0.0, 3.0),
    ("n30_0", 10.5, 30.0, 0.0, 2.0),
    ("n30_45", 10.5, 30.0, 45.0, 1.0),
    ("n30_90", 10.5, 30.0, 90.0, 1.0),
]
CORNER_POINTS = [(0.0, 0.0), (30.0, 45.0)]
EDGES = [(10.25, 0.0, 0.0), (10.25, 30.0, 0.0), (10.75, 0.0, 0.0), (10.75, 30.0, 0.0)]


def kits():
    nom = madp_kit()
    return {
        "nominal": nom,
        "C_T+20%": DiodeKit(nom.name, nom.on, nom.off.scaled(C_T=1.2)),
        "C_T-10%(assumed)": DiodeKit(nom.name, nom.on, nom.off.scaled(C_T=0.9)),
        "L_s+0.1nH(assumed)": madp_kit(0.1e-9),
    }


def make_cell(x):
    r, la, wa, ld, wd = [float(v) * 1e-3 for v in x]
    c = Cp2Cell(period=PERIOD, circumradius=r, stub_width=wa, line=la, line_diag=ld, stub_width_diag=wd)
    room_axis = PERIOD / 2 - 0.25e-3 - 1e-6 - (c.stub_extent() - c.line)
    # the tip corners of a rotated stub rectangle (half width w/2) must stay inside the cell: r_end cos45 + (w/2) sin45 <= P/2 - 0.25 mm
    room_diag = ((PERIOD / 2 - 0.25e-3 - 1e-6) - 0.5 * wd * math.sin(math.pi / 4)) / math.cos(math.pi / 4) - (
        c.apothem + 2 * c.pad + c.diode_span
    )
    return replace(c, line=max(0.25e-3, min(c.line, room_axis)), line_diag=max(0.25e-3, min(ld, room_diag)))


def solve_task(task):
    x, f, th, ph, kit_names = task
    cell = make_cell(x)
    ks = kits()
    sets = [cell.unit_cell(st, ks[k]).loads for k in kit_names for st in cell.states()]
    base = cell.unit_cell(cell.states()[0], ks["nominal"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = compute_s_parameters_oblique_over_loads(
            base,
            f * 1e9,
            MeshSettings(MESH, MESH, 2 * MESH),
            PL,
            math.radians(th),
            math.radians(ph),
            sets,
            warn_grating=False,
        )
    return np.array([r.s11 for r in res])  # (len(kit_names) * 4, 2, 2)


def violation(m, edge=False):
    if edge:  # band edges: relaxed targets, no opposite-handed term
        return max(0.0, m["step"] - 25.0) / 25.0 + max(0.0, -2.0 - m["amp_min"]) / 2.0
    return (
        max(0.0, m["dev"] - 15.0) / 15.0
        + max(0.0, m["step"] - 20.0) / 20.0
        + max(0.0, -1.5 - m["amp_min"]) / 1.5
        + max(0.0, m["spread"] - 1.0)
        + max(0.0, m["opp"] + 15.0) / 5.0
        + m["step"] / 200.0
        + m["dev"] / 600.0
    )


class Evaluator:
    def __init__(self, delta_path=None):
        self.lib = StateLibrary.load(ROOT / "results/phase20_r4/cp2_library_v1.npz")
        self.delta = DeltaTable.load(delta_path or os.environ.get("R3B_DELTA") or OUT / "hfss_delta_table.npz")

    def correction(self, f, th, ph):
        """Baseline mesh correction (c_corr - c_raw) and, at phi = 0, the HFSS delta; both (4, 2, 2) in the structure frame."""
        c = np.zeros((4, 2, 2), dtype=complex)
        for s in range(4):
            cc = self.lib.structure_circular(f, th, ph, s + 1, key="c_corr")[0]
            cr = self.lib.structure_circular(f, th, ph, s + 1, key="c_raw")[0]
            c[s] = cc - cr
            if ph == 0.0:
                c[s] += self.delta.at(f, th, s + 1)[0]
        return c

    def circ(self, s11, f, th, ph, k):
        """Corrected structure-frame circular responses of kit block k (4 states)."""
        c = np.array([to_structure_frame(circular_reflection(s11[4 * k + i]), math.radians(ph)) for i in range(4)])
        return c + self.correction(f, th, ph)

    def evaluate(self, x, pool=None):
        tasks, keys = [], []
        for name, f, th, ph, _ in NOMINAL:
            corner = name in ("n00", "n30_45")
            kn = ["nominal"] + (["C_T+20%", "C_T-10%(assumed)", "L_s+0.1nH(assumed)"] if corner else [])
            tasks.append((tuple(x), f, th, ph, kn))
            keys.append((name, f, th, ph, kn))
        for f, th, ph in EDGES:
            tasks.append((tuple(x), f, th, ph, ["nominal"]))
            keys.append((f"e{f:g}_{th:g}", f, th, ph, ["nominal"]))
        res = pool.map(solve_task, tasks) if pool is not None else [solve_task(t) for t in tasks]
        J, parts, metrics = 0.0, {}, {}
        for (name, f, th, ph, kn), s11 in zip(keys, res, strict=True):
            for k, kname in enumerate(kn):
                c = self.circ(s11, f, th, ph, k)
                m = point_metrics(c)
                tag = f"{name}|{kname}"
                metrics[tag] = m
                if name.startswith("e"):
                    w, v = 0.5, violation(m, edge=True)
                elif kname == "nominal":
                    w, v = next(w for n, *_, w in NOMINAL if n == name), violation(m)
                else:
                    w, v = 1.0, violation(m)
                J += w * v
                parts[tag] = w * v
                if kname == "nominal" and name.startswith("e"):
                    metrics[tag]["rr_phase"] = np.degrees(np.angle(c[:, 0, 0])).tolist()
                if name == "n00" and kname == "nominal":
                    metrics[tag]["rr_phase"] = np.degrees(np.angle(c[:, 0, 0])).tolist()
        # lower-Q term: phase slope of the four r_RR between 10.25 and 10.75 GHz at theta = 0 [deg/GHz], normalised by 150
        p_lo = np.array(metrics["e10.25_0|nominal"]["rr_phase"])
        p_hi = np.array(metrics["e10.75_0|nominal"]["rr_phase"])
        slope = np.abs((p_hi - p_lo + 180.0) % 360.0 - 180.0) / 0.5
        q = float(slope.mean()) / 150.0
        J += q
        parts["Q_slope"] = q
        metrics["slope_deg_per_ghz"] = slope.tolist()
        return float(J), parts, metrics


def cached():
    out = {}
    if CACHE.exists():
        for line in CACHE.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                out[tuple(r["x"])] = r
    return out


_APPEND_LOCK = threading.Lock()


def append(rec):
    with _APPEND_LOCK, open(CACHE, "a") as fh:
        fh.write(json.dumps(rec) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def f_obj(x, ev, pool, cache, log):
    x = np.clip(np.asarray(x, dtype=float), [b[0] for b in BOUNDS], [b[1] for b in BOUNDS])
    key = tuple(round(float(v), 4) for v in x)
    pen = 0.0
    if key in cache:
        rec = cache[key]
    else:
        t0 = time.monotonic()
        try:
            J, parts, metrics = ev.evaluate(key, pool)
        except ValueError as e:  # geometry outside the cell / invalid: large penalty (never a crash)
            J, parts, metrics = 1e3, {"invalid": str(e)[:120]}, {}
        rec = {"x": list(key), "J": J, "parts": parts, "metrics": metrics, "seconds": time.monotonic() - t0}
        cache[key] = rec
        append(rec)
    log.append(rec["J"])
    return rec["J"] + pen


STARTS = {
    "baseline": BASELINE,
    "diag_long": [4.6, 1.748, 1.375, 3.0, 1.375],
    "diag_longer_wide": [4.6, 1.748, 1.375, 4.5, 2.0],
    "diag_short": [4.6, 1.748, 1.375, 1.0, 1.375],
    "diag_wide": [4.6, 1.748, 1.375, 1.748, 2.2],
    "small_r": [4.0, 1.748, 1.2, 2.5, 1.2],
}
STEPS = [0.3, 0.2, 0.2, 0.8, 0.4]


def trajectory(start, maxfev=120):
    """Replay a start from the cache only (Nelder-Mead is deterministic): the ordered list of distinct x evaluated by that start. Raises KeyError if the run has not covered it yet."""
    from scipy.optimize import minimize

    cache = cached()
    seen = []

    def f(x, *_):
        x = np.clip(np.asarray(x, dtype=float), [b[0] for b in BOUNDS], [b[1] for b in BOUNDS])
        key = tuple(round(float(v), 4) for v in x)
        if key not in cache:
            raise KeyError(key)  # the run had not evaluated this point yet: stop the replay here
        if key not in seen:
            seen.append(key)
        return cache[key]["J"]

    x0 = np.array(STARTS[start], dtype=float)
    simplex = np.vstack(
        [x0] + [x0 + np.eye(5)[i] * STEPS[i] * (1 if x0[i] + STEPS[i] <= BOUNDS[i][1] else -1) for i in range(5)]
    )
    with contextlib.suppress(KeyError):  # run still in progress: the replay stops at the first uncached point
        minimize(
            f,
            x0,
            method="Nelder-Mead",
            options={"maxfev": maxfev, "xatol": 0.01, "fatol": 1e-3, "initial_simplex": simplex},
        )
    return seen


def run_one(start, maxfev, ev, pool, cache):
    from scipy.optimize import minimize

    x0 = np.array(STARTS[start], dtype=float)
    simplex = np.vstack(
        [x0] + [x0 + np.eye(5)[i] * STEPS[i] * (1 if x0[i] + STEPS[i] <= BOUNDS[i][1] else -1) for i in range(5)]
    )
    log = []
    res = minimize(
        f_obj,
        x0,
        args=(ev, pool, cache, log),
        method="Nelder-Mead",
        options={"maxfev": maxfev, "xatol": 0.01, "fatol": 1e-3, "initial_simplex": simplex},
    )
    best = min((r for r in cache.values()), key=lambda r: r["J"])
    (OUT / f"opt_{start}.json").write_text(
        json.dumps(
            {
                "start": start,
                "x_final": list(map(float, res.x)),
                "J_final": float(res.fun),
                "nfev": int(res.nfev),
                "best_overall": best["x"],
                "J_best_overall": best["J"],
                "names": NAMES,
            },
            indent=1,
        )
    )
    print(start, "done: J", res.fun, "x", res.x, "nfev", res.nfev, flush=True)
    return res


def run(start, maxfev=120):
    OUT.mkdir(parents=True, exist_ok=True)
    ev = Evaluator()
    cache = cached()
    with Pool(int(os.environ.get("R3B_POOL", "8"))) as pool:
        run_one(start, maxfev, ev, pool, cache)


def run_many(starts, maxfev=120, workers=40):
    """Resume several starts from the cache in threads that share ONE process pool (`workers` processes, each single-threaded BLAS): the 8 condition solves of every evaluation of every start run concurrently."""
    from concurrent.futures import ThreadPoolExecutor

    OUT.mkdir(parents=True, exist_ok=True)
    ev = Evaluator()
    cache = cached()
    with Pool(workers) as pool, ThreadPoolExecutor(len(starts)) as ex:
        futs = [ex.submit(run_one, s, maxfev, ev, pool, cache) for s in starts]
        for f in futs:
            f.result()


if __name__ == "__main__":
    if sys.argv[1] == "eval":
        x = [float(v) for v in sys.argv[2:7]] if len(sys.argv) > 2 else BASELINE
        ev = Evaluator()
        with Pool(8) as pool:
            t0 = time.monotonic()
            J, parts, metrics = ev.evaluate(tuple(x), pool)
        print("J", J, "seconds", time.monotonic() - t0)
        for k, v in parts.items():
            print(k, round(v, 3))
        for k, m in metrics.items():
            if isinstance(m, dict):
                print(k, {a: round(b, 2) for a, b in m.items() if a != "rr_phase"})
        print("slope", metrics["slope_deg_per_ghz"])
    elif sys.argv[1] == "many":  # many <maxfev> <workers> <start> [<start> ...]
        run_many(sys.argv[4:], int(sys.argv[2]), int(sys.argv[3]))
    else:
        run(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 120)
