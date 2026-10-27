"""Phase 28 item 3: full characterisation of one diode-switched cell candidate at 10.0 GHz.

    python scripts/research/phase28_cell_characterise.py <spec.json> <out.json> [--workers 8]

`spec.json`: {"name": ..., "period_mm": 6.4, "family": "sector_hub", "params": {...SectorHubCell fields in metres...}, "kit": "madp"|"smp"|"ideal",
"meshes": [32, 48, 64], "i_forward_ma": 10}. Everything is evaluated with `coded_cell_eval` (normal-incidence tensor sheet, lumped diode loads):

* mesh convergence of every state (b_bar, delta, dX_aniso at the cell's own mean loading and at the design loading, alpha_beam) over `meshes` (the last two meshes are
  >= 48 on a 6.4 mm cell: 48, 64; 60, 80 on an 8.0 mm cell), relative change between the last two meshes;
* the four states on the converged mesh: balance (spread of b_bar and delta over the states) and the rotation check;
* diode tolerance: +-20 % on C_T, L_s, R_s of the kit (all three together in the worst-case directions and each alone), states 0 and 1;
* contrast bandwidth: the four states over 9.0 .. 11.0 GHz in 0.25 GHz steps;
* hardware: diodes per cell, forward diodes per state, forward current and DC power per cell (same formula as the item-2 report).
Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET (diode values), VIA_OMITTED_APPROXIMATION; MADP package inductance NOT_AVAILABLE (L_s = 0 assumed).
"""

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

F0 = 10.0e9
TARGET_OHM, STRETCH_OHM = 103.0, 132.0


def make_cell(spec):
    from src.hologram.coded_cell_candidates import ManhattanHubCell, RationalSectorCell, SectorHubCell

    families = {
        "sector_hub": SectorHubCell,
        "rational_sector_hub": RationalSectorCell,
        "manhattan_hub": ManhattanHubCell,
    }
    return families[spec.get("family", "sector_hub")](period=spec["period_mm"] * 1e-3, **spec["params"])


def make_kit(name, scale=None):
    """Diode kit; `scale` = {param: factor} applied to ON/OFF model parameters (tolerance study, ASSUMED_SENSITIVITY_NOT_DEVICE_DATA)."""
    from src.ris.cells import DiodeKit, ideal_kit, madp_kit, smp1345_kit

    if (
        scale and "madp_Ls_nH" in scale
    ):  # ASSUMED package inductance (the datasheet gives none): a sensitivity, not device data
        return madp_kit(l_s=scale["madp_Ls_nH"] * 1e-9)
    kit = {"madp": madp_kit, "smp": smp1345_kit, "ideal": ideal_kit}[name]()
    if not scale:
        return kit
    on, off = kit.on, kit.off
    f_on = {k: v for k, v in scale.items() if k in ("R_s", "L_s")}
    f_off = {k: v for k, v in scale.items() if k in ("L_s", "C_T")}
    if f_on:
        on = on.scaled(**f_on)
    if f_off and off.C_T is not None:
        off = off.scaled(**f_off)
    elif "L_s" in f_off:
        off = off.scaled(L_s=f_off["L_s"])
    return DiodeKit(f"{kit.name} x{scale}", on, off)


def evaluate(job):
    spec, state_idx, mesh, f_hz, scale, ffac = job
    warnings.simplefilter("ignore")
    from src.hologram.coded_cell_eval import principal_metrics, sheet_tensor

    cell = make_cell(spec)
    kit = make_kit(spec["kit"], scale)
    uc = cell.unit_cell(cell.states()[state_idx], kit)
    t0 = time.time()
    y = sheet_tensor(uc, f_hz, mesh, ffac)
    m = principal_metrics(y, uc.substrate, f_hz)
    m["seconds"] = time.time() - t0
    return {"state": state_idx, "mesh": mesh, "f_ghz": f_hz / 1e9, "scale": scale, "floquet_factor": ffac, **m}


CONV_TOL = 0.05  # relative change of a metric between the last two meshes accepted as converged (an ASSUMPTION, stated in the report)
CONV_TOL_ALPHA = 0.10  # alpha_beam ~ delta**2 times a smooth function of the loading, so a 5 % tolerance on delta is a 10 % tolerance on alpha_beam
METRICS = ("b_bar_ms", "delta_ms", "dx_aniso_ohm", "dx_aniso_at_design_loading_ohm", "alpha_beam_np_per_m")


def _get(rows, **kw):
    kw.setdefault("floquet_factor", 1)
    out = [r for r in rows if all(r.get(k) == v for k, v in kw.items())]
    if len(out) != 1:
        raise KeyError((kw, len(out)))
    return out[0]


def derive(doc):
    """Everything the report prints that is not a raw evaluation: convergence changes and verdicts, state balance, target verdicts, tolerance changes, bandwidth,
    hardware counts. (Computed here so that the report can be generated from the JSON alone.)"""
    spec, rows = doc["spec"], doc["rows"]
    meshes = spec["meshes"]
    m_last, m_prev = meshes[-1], meshes[-2]
    m_main = spec.get(
        "mesh_main", m_last
    )  # the mesh of the reported states, tolerance and bandwidth (can be below the finest convergence mesh)
    d = {
        "mesh_last": m_last,
        "mesh_main": m_main,
        "mesh_prev": m_prev,
        "convergence_tolerance": CONV_TOL,
        "convergence": {},
        "states": {},
        "target": {},
        "tolerance": [],
        "bandwidth": {},
        "hardware": {},
    }
    for met in METRICS:
        ch = []
        for s in range(4):
            a, b = (
                _get(rows, state=s, mesh=m_last, f_ghz=10.0, scale=None)[met],
                _get(rows, state=s, mesh=m_prev, f_ghz=10.0, scale=None)[met],
            )
            ch.append(abs(a - b) / max(abs(a), 1e-12))
        tol = CONV_TOL_ALPHA if met == "alpha_beam_np_per_m" else CONV_TOL
        d["convergence"][met] = {"max_relative_change": max(ch), "tolerance": tol, "converged": max(ch) <= tol}
    d["convergence_all"] = all(v["converged"] for v in d["convergence"].values())
    d["floquet"] = {"mesh": m_prev, "factor": 2, "max_relative_change": {}}
    for met in METRICS:
        ch = []
        for s in range(4):
            a, b = (
                _get(rows, state=s, mesh=m_prev, f_ghz=10.0, scale=None, floquet_factor=2)[met],
                _get(rows, state=s, mesh=m_prev, f_ghz=10.0, scale=None)[met],
            )
            ch.append(abs(a - b) / max(abs(a), 1e-12))
        d["floquet"]["max_relative_change"][met] = max(ch)
    d["floquet"]["converged"] = all(
        v <= (CONV_TOL_ALPHA if k == "alpha_beam_np_per_m" else CONV_TOL)
        for k, v in d["floquet"]["max_relative_change"].items()
    )
    last = [_get(rows, state=s, mesh=m_main, f_ghz=10.0, scale=None) for s in range(4)]
    for met in METRICS + ("beta_over_k0",):
        vals = [r[met] for r in last]
        d["states"][met] = {
            "min": min(vals),
            "max": max(vals),
            "spread_relative": (max(vals) - min(vals)) / max(abs(v) for v in vals),
        }
    t, st = doc["targets_ohm"]["target"], doc["targets_ohm"]["stretch"]
    for key in ("dx_aniso_ohm", "dx_aniso_at_design_loading_ohm"):
        mn = min(r[key] for r in last)
        d["target"][key] = {"min_over_states": mn, "meets_target": mn >= t, "meets_stretch": mn >= st}
    nominal = {s: _get(rows, state=s, mesh=m_main, f_ghz=10.0, scale=None) for s in (0, 1)}
    d["ls_sensitivity"] = []
    for r in rows:
        if r["scale"] is None or r["f_ghz"] != 10.0:
            continue
        nom = nominal[r["state"]]
        if "madp_Ls_nH" in r["scale"]:
            d["ls_sensitivity"].append(
                {
                    "madp_Ls_nH": r["scale"]["madp_Ls_nH"],
                    "state": r["state"],
                    **{
                        met: r[met]
                        for met in ("delta_ms", "dx_aniso_ohm", "dx_aniso_at_design_loading_ohm", "alpha_beam_np_per_m")
                    },
                    **{
                        "change_" + met: (r[met] - nom[met]) / abs(nom[met])
                        for met in ("delta_ms", "dx_aniso_ohm", "dx_aniso_at_design_loading_ohm", "alpha_beam_np_per_m")
                    },
                }
            )
            continue
        d["tolerance"].append(
            {
                "scale": r["scale"],
                "state": r["state"],
                **{
                    met: r[met]
                    for met in ("delta_ms", "dx_aniso_ohm", "dx_aniso_at_design_loading_ohm", "alpha_beam_np_per_m")
                },
                **{
                    "change_" + met: (r[met] - nom[met]) / abs(nom[met])
                    for met in ("delta_ms", "dx_aniso_ohm", "dx_aniso_at_design_loading_ohm", "alpha_beam_np_per_m")
                },
            }
        )
    d["tolerance_worst_change"] = {
        met: max(abs(x["change_" + met]) for x in d["tolerance"])
        for met in ("delta_ms", "dx_aniso_ohm", "dx_aniso_at_design_loading_ohm", "alpha_beam_np_per_m")
    }
    freqs = sorted(
        {r["f_ghz"] for r in rows if r["scale"] is None and r["mesh"] == m_main and r["floquet_factor"] == 1}
    )
    bw = []
    for f_ in freqs:
        rr = [_get(rows, state=s, mesh=m_main, f_ghz=f_, scale=None) for s in range(4)]
        bw.append(
            {
                "f_ghz": f_,
                "min_dx_aniso_ohm": min(r["dx_aniso_ohm"] for r in rr),
                "min_dx_design_ohm": min(r["dx_aniso_at_design_loading_ohm"] for r in rr),
                "min_alpha_np_per_m": min(r["alpha_beam_np_per_m"] for r in rr),
            }
        )
    d["bandwidth"]["rows"] = bw
    nom_dx = next(x for x in bw if x["f_ghz"] == 10.0)["min_dx_aniso_ohm"]
    for tag, key, thr in (
        ("above_target_design_loading", "min_dx_design_ohm", t),
        ("within_70pct_of_nominal_own_loading", "min_dx_aniso_ohm", 0.7 * nom_dx),
    ):
        ok = [x["f_ghz"] for x in bw]
        i0 = ok.index(10.0)
        lo = hi = i0
        good = [x[key] >= thr for x in bw]
        if good[i0]:
            while lo > 0 and good[lo - 1]:
                lo -= 1
            while hi < len(bw) - 1 and good[hi + 1]:
                hi += 1
            d["bandwidth"][tag] = {
                "threshold_ohm": thr,
                "f_low_ghz": bw[lo]["f_ghz"],
                "f_high_ghz": bw[hi]["f_ghz"],
                "width_ghz": bw[hi]["f_ghz"] - bw[lo]["f_ghz"],
            }
        else:
            d["bandwidth"][tag] = {"threshold_ohm": thr, "f_low_ghz": None, "f_high_ghz": None, "width_ghz": 0.0}
    n = doc["n_cells_150mm_hexagon"]
    fwd = max(doc["forward_diodes_per_state"])
    hw = {
        "diodes_per_cell": doc["n_diodes_per_cell"],
        "forward_diodes_per_state": fwd,
        "n_cells": n,
        "diodes_total": n * doc["n_diodes_per_cell"],
        "forward_diodes_total": n * fwd,
    }
    hw["dc_power_w_at_i_forward"] = hw["forward_diodes_total"] * doc["i_forward_ma"] * 1e-3 * doc["supply_rail_v"]
    hw["dc_power_w_at_1ma"] = hw["forward_diodes_total"] * 1e-3 * doc["supply_rail_v"]
    d["hardware"] = hw
    return d


def main(spec_path, out_path, workers):
    spec = json.loads(Path(spec_path).read_text(encoding="utf-8"))
    meshes = spec["meshes"]
    mesh_main = spec.get("mesh_main", meshes[-1])
    jobs = []
    for mesh in meshes:  # convergence: all four states
        for s in range(4):
            jobs.append((spec, s, mesh, F0, None, 1))
    for s in range(4):  # Floquet cutoff 2 n at the second-finest mesh
        jobs.append((spec, s, meshes[-2], F0, None, 2))
    tol_cases = [
        None,
        {"C_T": 0.8},
        {"C_T": 1.2},
        {"L_s": 0.8},
        {"L_s": 1.2},
        {"R_s": 0.8},
        {"R_s": 1.2},
        {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2},
        {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8},
    ]
    if spec["kit"] == "madp":
        tol_cases += [{"madp_Ls_nH": 0.3}, {"madp_Ls_nH": 0.7}]
    for sc in tol_cases[1:]:
        for s in (0, 1):
            jobs.append((spec, s, mesh_main, F0, sc, 1))
    freqs = [
        9.0e9 + 0.25e9 * i for i in range(9) if abs(9.0e9 + 0.25e9 * i - F0) > 1.0
    ]  # 10.0 GHz itself is the convergence row
    for f in freqs:
        for s in range(4):
            jobs.append((spec, s, mesh_main, f, None, 1))
    t0 = time.time()
    # resumable: every finished evaluation is appended to <out>.rows.jsonl; a rerun reuses them
    cache_path = Path(str(out_path) + ".rows.jsonl")
    done = {}
    if cache_path.exists():
        for line in cache_path.read_text(encoding="utf-8").splitlines():
            rec = json.loads(line)
            done[rec["key"]] = rec["row"]

    def key_of(job):
        return json.dumps([job[1], job[2], job[3], job[4], job[5]], sort_keys=True)

    todo = [j for j in jobs if key_of(j) not in done]
    print(f"{len(jobs)} evaluations, {len(done)} cached, {len(todo)} to run", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool, open(cache_path, "a", encoding="utf-8") as cache:
        for job, row in zip(todo, pool.map(evaluate, todo), strict=True):
            cache.write(json.dumps({"key": key_of(job), "row": row}) + chr(10))
            cache.flush()
            done[key_of(job)] = row
    rows = [done[key_of(j)] for j in jobs]
    Path(str(out_path) + ".raw.json").write_text(json.dumps({"spec": spec, "rows": rows}, indent=1), encoding="utf-8")
    kit_name = spec["kit"]
    cell = make_cell(spec)
    from src.hologram import coded_feasibility as cf

    n_cells = int(cf.lattice(cf.hexagon_m(0.150), spec["period_mm"] * 1e-3)[0].size)
    doc = {
        "label": "PHASE28_ITEM3_CELL_CHARACTERISATION",
        "labels": [
            "MODEL_ONLY",
            "NOT_HFSS_VALIDATED",
            "PROVISIONAL_DATASHEET",
            "VIA_OMITTED_APPROXIMATION",
            "IDEALISED_NORMAL_INCIDENCE_TENSOR_SHEET",
        ],
        "spec": spec,
        "targets_ohm": {"target": TARGET_OHM, "stretch": STRETCH_OHM},
        "frequency_ghz": F0 / 1e9,
        "kit": kit_name,
        "n_diodes_per_cell": cell.n_diodes,
        "forward_diodes_per_state": [sum(v == "ON" for v in st) for st in cell.states()],
        "rows": rows,
        "seconds": time.time() - t0,
        "n_cells_150mm_hexagon": n_cells,
        "supply_rail_v": 3.3,
        "i_forward_ma": spec.get("i_forward_ma", 10.0),
    }
    doc["derived"] = derive(doc)
    Path(out_path).write_text(json.dumps(doc, indent=1), encoding="utf-8")
    print("wrote", out_path, f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], int(a[a.index("--workers") + 1]) if "--workers" in a else 8)
