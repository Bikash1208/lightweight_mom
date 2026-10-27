"""Phase 28: loop cell R4 (outer 6.6 mm, strip 1.2 mm, 9.0 mm period, RT/duroid 5880 3.175 mm) with the MEASURED MA4AGP907 diode (`src/ris/diodes_measured.py`, MEASURED_VENDOR_SPARAM).
ON at 5 mA and 10 mA, OFF at 0 V, Z(f) taken from the extracted table (frequency-dependent LumpedLoad); the diode footprint is the cell's 0.6 x 1.2 mm region. MODEL_ONLY.

    python scripts/research/phase28_loop_measured.py solve <out.jsonl> [--workers N]
    python scripts/research/phase28_loop_measured.py report <rows.jsonl> <out.json>

States: ON_5mA, ON_10mA, OFF_0V (full measured Z(f)); ON_5mA_Ronly, ON_10mA_Ronly (constant real part at 10 GHz) and OFF_0V_Conly (ideal capacitor at the 10 GHz value): the HFSS-representable forms (R only / C only).
Scaled states (tolerance): R x 0.8 / 1.2 on an ON state, C x 0.8 / 1.2 on the OFF state.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts" / "research"))
import phase28_scalar_ring_cell as rc  # noqa: E402

VARIANT, SUBSTRATE = "R4_So6.6_w1.2", "rt5880_3p175"
FREQS = [9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0]


def load_for(state, scale):
    from src.loads import LumpedLoad
    from src.ris import diodes_measured as dm

    so, wr = 6.6e-3, 1.2e-3
    x, y, ln = 0.0, so / 2 - wr / 2, rc.GAP
    if state.endswith("Ronly"):
        r = dm.z_of_f(state.replace("_Ronly", ""))(10e9).real * (scale or {}).get("R", 1.0)
        return LumpedLoad(x=x, y=y, length=ln, width=wr, direction="x", R=float(r), transverse_conducts=False)
    if state.endswith("Conly"):
        zi = dm.z_of_f(state.replace("_Conly", ""))(10e9).imag / (scale or {}).get("C", 1.0)
        c = -1.0 / (2 * np.pi * 10e9 * zi)
        return LumpedLoad(x=x, y=y, length=ln, width=wr, direction="x", C=float(c), transverse_conducts=False)
    return dm.load(state, x, y, ln, wr, "x", scale)


def solve(job):
    state, mesh, f_hz, scale = job
    warnings.simplefilter("ignore")
    from dataclasses import replace

    from src.hologram.coded_cell_eval import sheet_tensor

    uc, _ = rc.build(VARIANT, SUBSTRATE, "MADP_Ls0", True)
    uc = replace(uc, loads=(load_for(state, scale),))
    y = sheet_tensor(uc, f_hz, mesh)
    return {"state": state, "mesh": mesh, "f_hz": f_hz, "scale": scale, "y_s": [[[float(v.real), float(v.imag)] for v in row] for row in y]}


def jobs():
    out = []
    for s in ("ON_5mA", "ON_10mA", "OFF_0V"):
        out += [(s, 60, f * 1e9, None) for f in FREQS] + [(s, 90, 10e9, None)]
    for s, sc in (("ON_5mA", "R"), ("ON_10mA", "R"), ("OFF_0V", "C")):
        out += [(s, 60, 10e9, {sc: 0.8}), (s, 60, 10e9, {sc: 1.2})]
    out += [("ON_5mA_Ronly", 60, 10e9, None), ("ON_10mA_Ronly", 60, 10e9, None), ("OFF_0V_Conly", 60, 10e9, None)]
    return out


def jobs_extra():
    """HFSS-representable (R only / C only) states over the package band, for the 10 mA package default, and mesh 90 at 10 GHz."""
    out = []
    for f in (9.5, 9.75, 10.25, 10.5):
        out += [("ON_10mA_Ronly", 60, f * 1e9, None), ("OFF_0V_Conly", 60, f * 1e9, None)]
    out += [("ON_10mA_Ronly", 90, 10e9, None), ("OFF_0V_Conly", 90, 10e9, None)]
    return out


def key(r):
    return json.dumps([r["state"], r["mesh"], r["f_hz"], r["scale"]], sort_keys=True)


def rows_of(p):
    return [json.loads(l) for l in Path(p).read_text(encoding="utf-8").splitlines()] if Path(p).exists() else []


def pair(rows, on, off, mesh, f, on_scale=None, off_scale=None):
    idx = {key(r): r for r in rows}
    a = idx.get(key({"state": on, "mesh": mesh, "f_hz": f, "scale": on_scale}))
    b = idx.get(key({"state": off, "mesh": mesh, "f_hz": f, "scale": off_scale}))
    if a is None or b is None:
        return None
    base = {"substrate": SUBSTRATE, "f_hz": f}
    return rc.pair_metrics({**base, "y_s": a["y_s"]}, {**base, "y_s": b["y_s"]})


def report(rows):
    out = {"labels": ["MODEL_ONLY", "MEASURED_VENDOR_SPARAM", "NOT_HFSS_VALIDATED"], "cases": {}}
    old = json.loads((ROOT / "results/phase28/f10/ring_cell/report.json").read_text(encoding="utf-8"))["confirm"]["R4_So6.6_w1.2|rt5880_3p175|MADP_Ls0|90|10000000000.0|null"]
    out["old_model_5p2ohm_25fF_mesh90"] = old
    for on, tag in (("ON_5mA", "5mA"), ("ON_10mA", "10mA")):
        c = {"mesh90_10ghz": pair(rows, on, "OFF_0V", 90, 10e9), "mesh60_10ghz": pair(rows, on, "OFF_0V", 60, 10e9),
             "band": {f"{f:g}": pair(rows, on, "OFF_0V", 60, f * 1e9) for f in FREQS},
             "tolerance": {"R0.8": pair(rows, on, "OFF_0V", 60, 10e9, {"R": 0.8}), "R1.2": pair(rows, on, "OFF_0V", 60, 10e9, {"R": 1.2}),
                           "C0.8": pair(rows, on, "OFF_0V", 60, 10e9, None, {"C": 0.8}), "C1.2": pair(rows, on, "OFF_0V", 60, 10e9, None, {"C": 1.2}),
                           "R0.8_C0.8": pair(rows, on, "OFF_0V", 60, 10e9, {"R": 0.8}, {"C": 0.8}), "R1.2_C1.2": pair(rows, on, "OFF_0V", 60, 10e9, {"R": 1.2}, {"C": 1.2})},
             "R_and_C_only": pair(rows, on + "_Ronly", "OFF_0V_Conly", 60, 10e9)}
        out["cases"][tag] = c
    return out


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] in ("solve", "solve-extra"):
        done = {key(r) for r in rows_of(a[1])}
        todo = [j for j in (jobs() if a[0] == "solve" else jobs_extra()) if key({"state": j[0], "mesh": j[1], "f_hz": j[2], "scale": j[3]}) not in done]
        print(len(todo), "solves", flush=True)
        with ProcessPoolExecutor(max_workers=int(a[a.index("--workers") + 1]) if "--workers" in a else 40) as ex, open(a[1], "a") as fh:
            for r in ex.map(solve, todo):
                fh.write(json.dumps(r) + "\n"); fh.flush()
    else:
        Path(a[2]).write_text(json.dumps(report(rows_of(a[1])), indent=1, default=float), encoding="utf-8")
        print("wrote", a[2])
