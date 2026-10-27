"""Phase 28 (cell redesign, <= 2 diodes): topology screen in the rational-sector family at mesh 32 (relative screening, NOT converged values).

    python scripts/research/phase28_two_diode_screen.py <out.jsonl> [--workers N]

The eight neck sites of `RationalSectorCell` get a role each: two sites hold the diodes A and B (MADP-000907-14020, L_s = 0), every other site is a static metal bridge ('M', ideal ON)
or a static gap ('G', ideal OFF). The four states are the on/off combinations (A, B) in {0,1}^2. For every (base geometry, roles, state) the normal-incidence tensor sheet is solved
at mesh 32 and reduced to (b_bar, delta, axis, alpha_beam at the cell's OWN loading). Resumable (one JSON line per solve). Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET."""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import itertools  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import warnings  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

F_HZ = 10.0e9
MESH = 32
G_PARAMS = {"hub_half": 1.2, "neck_gap": 0.4, "neck_width": 0.4, "gap_t": 0.2, "arm_outer": 3.0, "diag_outer": 2.8, "clip": 3.0}  # mm
BASES = {
    "G": G_PARAMS,
    "H": {**G_PARAMS, "hub_half": 1.0, "diag_outer": 3.0},
}


def configs():
    for i, j in itertools.combinations(range(8), 2):
        others = [k for k in range(8) if k not in (i, j)]
        for fill in itertools.product("MG", repeat=6):
            roles = ["?"] * 8
            roles[i], roles[j] = "A", "B"
            for k, r in zip(others, fill, strict=True):
                roles[k] = r
            yield "".join(roles)


def unit_cell(params_mm, roles, a, b, period_mm=6.4, l_s_nh=None, substrate=None):
    from src.hologram.coded_cell_candidates import RationalSectorCell
    from src.ris.cells import ideal_kit, madp_kit
    from src.ris.diodes import diode_load
    from src.geometry import UnitCell
    from src.materials import Metal
    from src.periodic_bc import RectangularLattice
    from src.ris.cells import RT5880

    cell = RationalSectorCell(period=period_mm * 1e-3, mode="bowtie", **{k: v * 1e-3 for k, v in params_mm.items()})
    mk = madp_kit(None if l_s_nh is None else l_s_nh * 1e-9)
    ik = ideal_kit()
    loads = []
    for r, (x, y, length, width, d) in zip(roles, cell.diode_sites(), strict=True):
        model = {"A": mk.on if a else mk.off, "B": mk.on if b else mk.off, "M": ik.on, "G": ik.off}[r]
        loads.append(diode_load(model, x, y, length, width, direction=d, transverse_conducts=False))
    return UnitCell(RectangularLattice(cell.period, cell.period), cell.geometry(), substrate or RT5880, Metal(), backing="pec", loads=tuple(loads))


def solve(job):
    base, roles, a, b = job
    warnings.simplefilter("ignore")
    from src.hologram.coded_cell_eval import principal_metrics, sheet_tensor

    uc = unit_cell(BASES[base], roles, a, b)
    y = sheet_tensor(uc, F_HZ, MESH)
    m = principal_metrics(y, uc.substrate, F_HZ)
    keep = ("b_bar_ms", "delta_ms", "axis_deg", "beta_over_k0", "dx_aniso_ohm", "dx_aniso_at_design_loading_ohm", "alpha_beam_np_per_m", "g_max_ms")
    return {"base": base, "roles": roles, "a": a, "b": b, **{k: m.get(k) for k in keep}}


def main(out, workers):
    out = Path(out)
    done = set()
    if out.exists():
        for line in out.read_text(encoding="utf-8").splitlines():
            r = json.loads(line)
            done.add((r["base"], r["roles"], r["a"], r["b"]))
    jobs = [(base, roles, a, b) for base in BASES for roles in configs() for a in (0, 1) for b in (0, 1) if (base, roles, a, b) not in done]
    print(len(jobs), "solves to run", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex, out.open("a", encoding="utf-8") as fh:
        for n, r in enumerate(ex.map(solve, jobs, chunksize=8)):
            fh.write(json.dumps(r) + "\n")
            if n % 200 == 0:
                fh.flush()
                print(n, flush=True)


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], int(a[a.index("--workers") + 1]) if "--workers" in a else 32)
