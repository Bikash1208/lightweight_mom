"""Run the sweep defined in a config file and write the unit-cell database.

Run from the project root:
    python examples/run_sweep.py [examples/sweep_config.yaml] [results/sweep_example]

Re-running the same command resumes from results/.../checkpoint.jsonl (finished results are
reused, so an interrupted run continues where it stopped); pass --restart to discard it.
Outputs: checkpoint.jsonl, database.json (complete), database.csv (flat table).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import units
from src.config import load_config
from src.database import Database
from src.runner import run_sweep
from src.sweep import plan_sweep

ROOT = Path(__file__).resolve().parent.parent


def main(argv):
    restart = "--restart" in argv
    args = [a for a in argv if not a.startswith("--")]
    config_path = Path(args[0]) if args else ROOT / "examples" / "sweep_config.yaml"
    out_dir = Path(args[1]) if len(args) > 1 else ROOT / "results" / "sweep_example"

    cfg = load_config(config_path)
    plan = plan_sweep(cfg)
    f = plan.frequencies
    print(f"Sweep: {len(plan.cells)} cell(s) ({len(plan.valid_cells)} valid) x {f.size} frequencies "
          f"({units.to_ghz(f[0]):.6g}..{units.to_ghz(f[-1]):.6g} GHz) = {plan.n_solves} solves")

    total = plan.n_solves
    def progress(e):
        if e["done"] == total or e["done"] % 10 == 0:
            print(f"  {e['done']}/{total} results", end="\r", flush=True)

    result = run_sweep(cfg, out_dir, resume=not restart, progress=progress)
    db = result.database
    print(f"\ncomputed {result.computed}, reused from checkpoint {result.reused}, "
          f"dropped unreadable checkpoint lines {result.dropped}")
    print(f"wrote {result.json_path}\n      {result.csv_path}")

    last = len(db.frequencies) - 1
    print(f"\nS21 at {units.to_ghz(db.frequencies[last]):.6g} GHz (magnitude, phase deg; exp(+jwt); "
          f"planes at the substrate faces)")
    print("  length[mm] width[mm] |  xx: mag  phase |  yx: mag | RR: mag  phase | LR: mag | "
          "dFloq   dMesh  ok")
    for c in db.cells:
        if c.status != "ok":
            print(f"  cell {c.index}: SKIPPED ({c.reason})")
            continue
        r = db.record(c.index, last)
        ref = c.validation["refinement"]
        flag = "yes" if (r.diagnostics["checks_passed"] and (ref is None or ref["passed"])) else "NO"
        errs = f"{ref['floquet_error']:.4f} {ref['mesh_error']:.4f}" if ref else "  -      -  "
        print(f"  {units.to_mm(c.parameters['length_m']):9.3g} {units.to_mm(c.parameters['width_m']):8.3g} | "
              f"{r.magnitude('S21_xx'):.4f} {r.phase_deg('S21_xx'):8.2f} | {r.magnitude('S21_yx'):.4f} | "
              f"{r.magnitude('S21_RR'):.4f} {r.phase_deg('S21_RR'):8.2f} | {r.magnitude('S21_LR'):.4f} | {errs}  {flag}")
    failed = db.failed_validation()
    print(f"\ncells failing validation: {len(failed)}")

    back = Database.from_json(result.json_path)      # the file round-trips exactly
    assert back == db


if __name__ == "__main__":
    main(sys.argv[1:])
