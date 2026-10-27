"""Phase 28 (HFSS-first path, item 4): screening of variants of the cell G and its package-inductance sensitivity, 10.0 GHz.

    python scripts/research/phase28_g_refine.py <out.jsonl> <jobs.json> [--workers N]

`jobs.json`: list of {"name": ..., "params": {overrides of G's RationalSectorCell fields in mm}, "mesh": 64, "ls_nh": null|x, "states": [0, 1]}.
States 0 and 1 stand for all four (C4 pinwheel: 0/2 and 1/3 are equal in the G library). Every finished evaluation is one line of `out.jsonl` (resumable). Screening is at mesh 64
(G's own 64 -> 96 change at 10.0 GHz is 1.3 % in dX); only the variant that is kept gets a mesh-96 confirmation. Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET.
Variant dimensions are multiples of 0.2 mm so that the rooftop staircase is self-similar on both mesh 64 (0.1 mm) and mesh 96 (0.0667 mm)."""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

SPEC = json.loads((ROOT / "results/phase28/f10/cells/spec_G_6p4_bowtie.json").read_text(encoding="utf-8"))
F_HZ = 10.0e9


def run(job):
    warnings.simplefilter("ignore")
    from phase28_cell_characterise import evaluate

    spec = json.loads(json.dumps(SPEC))
    for k, v in (job.get("params") or {}).items():
        spec["params"][k] = v * 1e-3
    scale = {"madp_Ls_nH": job["ls_nh"]} if job.get("ls_nh") is not None else None
    t0 = time.time()
    r = evaluate((spec, job["state"], job["mesh"], F_HZ, scale, 1))
    r.update({"name": job["name"], "params_mm": job.get("params") or {}, "ls_nh": job.get("ls_nh"), "wall_s": time.time() - t0})
    return r


def main(out, jobs_file, workers):
    out = Path(out)
    done = set()
    if out.exists():
        for line in out.read_text(encoding="utf-8").splitlines():
            r = json.loads(line)
            done.add((r["name"], r["state"], r["mesh"]))
    jobs = []
    for v in json.loads(Path(jobs_file).read_text(encoding="utf-8")):
        for s in v.get("states", [0, 1]):
            j = {**v, "state": s, "mesh": v.get("mesh", 64)}
            if (j["name"], s, j["mesh"]) not in done:
                jobs.append(j)
    print(len(jobs), "jobs", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex, out.open("a", encoding="utf-8") as fh:
        for r in ex.map(run, jobs):
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            print("done", r["name"], r["state"], r["mesh"], f"{r['wall_s']:.0f}s", flush=True)


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], int(a[a.index("--workers") + 1]) if "--workers" in a else 8)
