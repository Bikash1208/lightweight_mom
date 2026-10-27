"""Phase 28 (cell redesign, <= 2 diodes): mesh-64 follow-up of the best two-diode topology of the screen, G-family variants, thicker boards and the package-inductance sweep, plus the
8-diode G on the thicker boards as a reference. python scripts/research/phase28_two_diode_followup.py <out.jsonl> [--workers N]
Substrates are LISTED options (RT/duroid 5880 1.575 mm = the board of every earlier Phase; RT/duroid 5880 3.175 mm; RO4003C 1.524 mm), none is chosen. MODEL_ONLY."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import json, sys, warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts" / "research"))
BEST = "MMGAMMBG"
SUBS = {"rt5880_1p575": (2.2, 0.0009, 1.575e-3), "rt5880_3p175": (2.2, 0.0009, 3.175e-3), "ro4003c_1p524": (3.38, 0.0027, 1.524e-3)}

def jobs():
    import phase28_two_diode_screen as s
    G = s.BASES["G"]
    out = []
    for sub in SUBS:
        for a in (0, 1):
            for b in (0, 1):
                out.append({"name": f"T1_{BEST}_{sub}", "kind": "two", "roles": BEST, "params": G, "sub": sub, "a": a, "b": b, "ls": None})
    var = {"V1_hub1p0": {"hub_half": 1.0}, "V2_diag3p0": {"diag_outer": 3.0}, "V3_hub1p0_diag3p0": {"hub_half": 1.0, "diag_outer": 3.0}, "V4_neckw0p6": {"neck_width": 0.6}}
    for n, o in var.items():
        for a in (0, 1):
            for b in (0, 1):
                out.append({"name": f"{n}_{BEST}", "kind": "two", "roles": BEST, "params": {**G, **o}, "sub": "rt5880_1p575", "a": a, "b": b, "ls": None})
    for ls in (0.1, 0.2, 0.3):
        for a in (0, 1):
            for b in (0, 1):
                out.append({"name": f"LS_{ls}nH_{BEST}", "kind": "two", "roles": BEST, "params": G, "sub": "rt5880_1p575", "a": a, "b": b, "ls": ls})
    for sub in ("rt5880_3p175", "ro4003c_1p524"):
        for st in (0, 1):
            out.append({"name": f"REF_G8diode_{sub}", "kind": "g8", "state": st, "sub": sub})
    return out

def run(job):
    warnings.simplefilter("ignore")
    import phase28_two_diode_screen as s
    from src.hologram.coded_cell_eval import principal_metrics, sheet_tensor
    from src.materials import Substrate
    er, td, h = SUBS[job["sub"]]
    sub = Substrate(job["sub"], er, td, h)
    if job["kind"] == "two":
        uc = s.unit_cell(job["params"], job["roles"], job["a"], job["b"], l_s_nh=job["ls"], substrate=sub)
    else:
        from phase28_cell_characterise import make_cell, make_kit
        from dataclasses import replace
        spec = json.loads((ROOT / "results/phase28/f10/cells/spec_G_6p4_bowtie.json").read_text())
        cell = make_cell(spec)
        uc = replace(cell.unit_cell(cell.states()[job["state"]], make_kit("madp")), substrate=sub)
    y = sheet_tensor(uc, 10.0e9, 64)
    m = principal_metrics(y, sub, 10.0e9)
    return {**{k: v for k, v in job.items() if k != "params"}, "params_mm": job.get("params"), **{k: m.get(k) for k in ("b_bar_ms", "delta_ms", "axis_deg", "beta_over_k0", "dx_aniso_ohm", "dx_aniso_at_design_loading_ohm", "alpha_beam_np_per_m")}}

if __name__ == "__main__":
    a = sys.argv[1:]
    js = jobs()
    print(len(js), "solves", flush=True)
    with ProcessPoolExecutor(max_workers=int(a[a.index("--workers") + 1]) if "--workers" in a else 40) as ex, open(a[0], "a") as fh:
        for r in ex.map(run, js):
            fh.write(json.dumps(r) + "\n"); fh.flush()
