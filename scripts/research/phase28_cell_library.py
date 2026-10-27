"""Phase 28 (HFSS-first path): state library of the cell G (6.4 mm rational-sector hub, 'bowtie' mode, MADP-000907-14020, L_s = 0) at 9.5 .. 10.5 GHz.

    python scripts/research/phase28_cell_library.py <out_dir> [--workers N]

Per state (0 .. 3 = slot axis 0, 45, 90, 135 deg) and frequency (9.5, 9.75, 10.0, 10.25, 10.5 GHz) the normal-incidence tensor sheet is solved on mesh 64 and mesh 96
(Floquet cutoff = mesh, the item-3 setting); both the 2 x 2 reflection matrix (x/y) and the sheet admittance tensor are stored. At 10.0 GHz an oblique solve
(theta 30 deg, phi 225 deg, TE/TM) is stored for the single-cell HFSS check. Resumable: every finished job is one line of `rows.jsonl`.
Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET (diode values), package inductance L_s = 0 nH (NOT_AVAILABLE_FROM_CURRENT_DATA, assumed zero)."""

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
FREQS_GHZ = [9.5, 9.75, 10.0, 10.25, 10.5]
MESHES = [64, 96]
OBLIQUE = {"theta_deg": 30.0, "phi_deg": 225.0, "f_ghz": 10.0, "mesh": 64}
CONV_TOL, CONV_TOL_ALPHA = 0.05, 0.10  # as the item-3 characterisation
METRICS = ("b_bar_ms", "delta_ms", "dx_aniso_ohm", "dx_aniso_at_design_loading_ohm", "alpha_beam_np_per_m")


def _c(z):
    return [float(z.real), float(z.imag)]


def job_normal(job):
    state, mesh, f_ghz = job
    warnings.simplefilter("ignore")
    import numpy as np
    from phase28_cell_characterise import make_cell, make_kit

    from src.basis_functions import MeshSettings
    from src.hologram.coded_cell_eval import principal_metrics
    from src.dispersion.tensor_sheet import reflection_matrix, sheet_admittance_tensor_from_gamma

    cell = make_cell(SPEC)
    uc = cell.unit_cell(cell.states()[state], make_kit(SPEC["kit"], None))
    t0 = time.time()
    g = reflection_matrix(uc, f_ghz * 1e9, MeshSettings(mesh, mesh, mesh))
    y = sheet_admittance_tensor_from_gamma(g, uc.substrate, f_ghz * 1e9)
    m = principal_metrics(y, uc.substrate, f_ghz * 1e9)
    return {
        "kind": "normal",
        "state": state,
        "mesh": mesh,
        "f_ghz": f_ghz,
        "gamma_xy": [[_c(g[i, j]) for j in range(2)] for i in range(2)],
        "y_s": [[_c(np.asarray(y)[i, j]) for j in range(2)] for i in range(2)],
        "seconds": time.time() - t0,
        **m,
    }


def job_oblique(job):
    state = job[1]
    warnings.simplefilter("ignore")
    import math

    import numpy as np
    from phase28_cell_characterise import make_cell, make_kit

    from src.basis_functions import MeshSettings
    from src.sparameters import ReferencePlanes, compute_s_parameters_oblique

    cell = make_cell(SPEC)
    uc = cell.unit_cell(cell.states()[state], make_kit(SPEC["kit"], None))
    mesh = OBLIQUE["mesh"]
    planes = ReferencePlanes(top_air_thickness=0.0, substrate_thickness=uc.substrate.thickness, bottom_air_thickness=None)
    t0 = time.time()
    sp = compute_s_parameters_oblique(
        uc,
        OBLIQUE["f_ghz"] * 1e9,
        MeshSettings(mesh, mesh, mesh),
        planes,
        math.radians(OBLIQUE["theta_deg"]),
        math.radians(OBLIQUE["phi_deg"]),
        warn_grating=False,
        warn_surface_wave=False,
    )
    s11 = np.asarray(sp.s11, dtype=complex)
    return {
        "kind": "oblique",
        "state": state,
        "mesh": mesh,
        "f_ghz": OBLIQUE["f_ghz"],
        "theta_deg": OBLIQUE["theta_deg"],
        "phi_deg": OBLIQUE["phi_deg"],
        "basis": "TE/TM (compute_s_parameters_oblique)",
        "s11_te_tm": [[_c(s11[i, j]) for j in range(2)] for i in range(2)],
        "seconds": time.time() - t0,
    }


def run(job):
    return job_oblique(job) if job[0] == "obl" else job_normal(job)


def _key(r):
    return ("obl", r["state"]) if r["kind"] == "oblique" else (r["state"], r["mesh"], r["f_ghz"])


def derive(rows):
    """Convergence mesh 64 -> 96 at every frequency (max over states), the criterion of the item-3 characterisation."""
    out = {}
    for f in FREQS_GHZ:
        ch = {}
        for met in METRICS:
            worst = 0.0
            for s in range(4):
                a = next(r for r in rows if r["kind"] == "normal" and (r["state"], r["mesh"], r["f_ghz"]) == (s, 96, f))[met]
                b = next(r for r in rows if r["kind"] == "normal" and (r["state"], r["mesh"], r["f_ghz"]) == (s, 64, f))[met]
                worst = max(worst, abs(a - b) / max(abs(a), 1e-12))
            tol = CONV_TOL_ALPHA if met == "alpha_beam_np_per_m" else CONV_TOL
            ch[met] = {"max_relative_change": worst, "tolerance": tol, "converged": worst <= tol}
        out[f"{f:g}"] = {"metrics": ch, "converged": all(v["converged"] for v in ch.values())}
    return out


def main(out_dir, workers):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows_path = out / "rows.jsonl"
    rows = [json.loads(line) for line in rows_path.read_text(encoding="utf-8").splitlines()] if rows_path.exists() else []
    done = {_key(r) for r in rows}
    jobs = [(s, m, f) for m in MESHES[::-1] for f in FREQS_GHZ for s in range(4)] + [("obl", s) for s in range(4)]
    todo = [j for j in jobs if (j if j[0] == "obl" else j) not in done]
    print(f"{len(todo)} jobs to run, {len(done)} done", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex, rows_path.open("a", encoding="utf-8") as fh:
        for r in ex.map(run, todo):
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            rows.append(r)
            print("done", _key(r), f"{r['seconds']:.0f}s", flush=True)
    doc = {
        "label": "PHASE28_G_STATE_LIBRARY",
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "PROVISIONAL_DATASHEET", "VIA_OMITTED_APPROXIMATION"],
        "package_inductance_nh": 0.0,
        "package_inductance_note": "NOT_AVAILABLE_FROM_CURRENT_DATA; L_s = 0 nH assumed (no datasheet value)",
        "spec": SPEC,
        "frequencies_ghz": FREQS_GHZ,
        "meshes": MESHES,
        "floquet_cutoff": "equal to the mesh (factor 1); the factor-2 check at mesh 64 and 10.0 GHz is in cells/cell_G_6p4_bowtie.json (derived.floquet)",
        "convergence_mesh_64_to_96": derive(rows),
        "rows": sorted(rows, key=lambda r: (r["kind"], r["state"], r["mesh"], r["f_ghz"])),
    }
    (out / "g_library.json").write_text(json.dumps(doc, indent=1), encoding="utf-8")
    print("wrote", out / "g_library.json")


if __name__ == "__main__":
    args = sys.argv[1:]
    w = int(args[args.index("--workers") + 1]) if "--workers" in args else 16
    main(args[0], w)
