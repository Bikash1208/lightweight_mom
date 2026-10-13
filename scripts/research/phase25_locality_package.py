"""Phase 25 follow-up 3: the HFSS neighbour-coupling package `hfss_exports/cell_pair_locality_check/`.

    python scripts/research/phase25_locality_package.py [pair name ...]   # default: all pairs

For the CURRENT pair (6.4 mm cell, empty | 5.85 mm square, the Phase 21-23 boards) and the RECOMMENDED pair
(8.0 mm cell, empty | 6.85 mm square, `results/phase25/gap_pair_study.json`): single square cell, 1x2 (square |
empty along x) and 2x2 (checkerboard) supercells as STEP, plus the MoM reference CSV (plane wave, normal
incidence, (0,0)-order S11_xx at the patch plane and at the 5 mm library plane, 10.33 / 10.54 / 10.75 GHz).
Every MoM number is a full supercell solve; the forward (isolated-cell) prediction of the same quantity is in
the same row, so the HFSS run decides directly whether the neighbour coupling is as large as the MoM says.
"""

import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase25_block_locality as bl  # noqa: E402

from src.export.array_export import write_hfss_step  # noqa: E402
from src.frequency import free_space_wavenumber  # noqa: E402

THICK_MM = 1.575
FREQS = (10.33e9, 10.54e9, 10.75e9)
PAIRS = {
    # name: (period_mm, HFSS square side mm, MoM mesh cells per period, mesh cells across the square)
    "current_P6.4_L5.85": (6.4, 5.85, 48, 44),  # the mesh represents 5.867 mm
    "candidate_separation_best_P8.0_L6.85": (8.0, 6.85, 28, 24),  # the mesh represents 6.857 mm; 172 deg, poor locality
    "candidate_locality_best_P6.4_L5.33": (6.4, 5.33, 48, 40),  # the 1.07 mm-gap pair; ~65 deg, best locality
}
CASES = {
    "1x2_square_empty": ((1, 2), np.array([[1, 0], [1, 0]])),  # (n_y, n_x) STEP grid; MoM: rows identical
    "2x2_checkerboard": ((2, 2), np.array([[1, 0], [0, 1]])),
}


def step_grid(shape, lev, side_mm):
    ny, nx = shape
    return [
        [(side_mm * 1e-3, side_mm * 1e-3) if lev[i % lev.shape[0], j % lev.shape[1]] else (0.0, 0.0) for j in range(nx)]
        for i in range(ny)
    ]


def main(out_dir, only=None):
    import json

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    manifest = {}
    if only and (out / "mom_reference.csv").exists():  # add pairs to an existing package
        with open(out / "mom_reference.csv", newline="") as fh:
            rows = [r for r in csv.DictReader(fh) if r["pair"] not in only]
        manifest = json.loads((out / "manifest.json").read_text())
    for pair, (period, side, mesh, cells) in PAIRS.items():
        if only and pair not in only:
            continue
        pd = out / pair
        pd.mkdir(exist_ok=True)
        p_m = period * 1e-3
        files = {}
        single = [[(side * 1e-3, side * 1e-3)]]
        write_hfss_step(single, p_m, p_m, THICK_MM * 1e-3, pd / "single_square_cell.step")
        files["single_square_cell"] = "single_square_cell.step"
        for case, (shape, lev) in CASES.items():
            write_hfss_step(step_grid(shape, lev, side), p_m, p_m, THICK_MM * 1e-3, pd / f"{case}.step")
            files[case] = f"{case}.step"
        manifest[pair] = {"period_mm": period, "square_side_mm": side, "files": files}
        for f in FREQS:
            k0 = free_space_wavenumber(f)
            to_top = np.exp(-2j * k0 * 5e-3)
            res = {
                name: bl.run(mesh, cells, {name: lev}, 2, f, tag=f"{pair} {f / 1e9:.2f} ", period_mm=period)
                for name, (_, lev) in CASES.items()
            }
            any_res = next(iter(res.values()))
            for tag, g in (
                ("isolated_square", any_res["isolated_square"]),
                ("isolated_empty", any_res["isolated_empty"]),
            ):
                gc = g[0] * np.exp(1j * np.radians(g[1]))
                rows.append(
                    {
                        "pair": pair,
                        "case": tag,
                        "frequency_GHz": f / 1e9,
                        "period_mm": period,
                        "hfss_square_side_mm": side,
                        "mom_square_side_mm": any_res["square_side_mm"],
                        "mom_mesh_cells_per_period": mesh,
                        "S11_xx_mag": abs(gc),
                        "S11_xx_phase_deg_ref0mm": np.degrees(np.angle(gc)),
                        "S11_xx_phase_deg_ref5mm": np.degrees(np.angle(gc * to_top)),
                        "isolated_model_mag": "",
                        "isolated_model_phase_deg_ref0mm": "",
                    }
                )
            for name, r in res.items():
                o = r["cases"][name]["orders"]["0,0"]
                gc = o["mom"][0] * np.exp(1j * np.radians(o["mom"][1]))
                rows.append(
                    {
                        "pair": pair,
                        "case": name,
                        "frequency_GHz": f / 1e9,
                        "period_mm": period,
                        "hfss_square_side_mm": side,
                        "mom_square_side_mm": r["square_side_mm"],
                        "mom_mesh_cells_per_period": mesh,
                        "S11_xx_mag": abs(gc),
                        "S11_xx_phase_deg_ref0mm": np.degrees(np.angle(gc)),
                        "S11_xx_phase_deg_ref5mm": np.degrees(np.angle(gc * to_top)),
                        "isolated_model_mag": o["forward"][0],
                        "isolated_model_phase_deg_ref0mm": o["forward"][1],
                    }
                )
    with open(out / "mom_reference.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1))
    print(f"wrote {len(rows)} reference rows to {out / 'mom_reference.csv'}")


if __name__ == "__main__":
    main(str(ROOT / "hfss_exports" / "cell_pair_locality_check"), only=sys.argv[1:] or None)
