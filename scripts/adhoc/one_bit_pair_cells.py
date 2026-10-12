"""Single-cell HFSS check package for the two 10.54 GHz 1-bit cell pairs (STATUS.md open step 2).

Pairs (both RT/duroid 5880, 1.575 mm, PEC ground, normal incidence):
  coded4   -- the 4-port coded boards' pair (hfss_exports/four_port_1bit_gap*/): 6.4 mm cell,
              state 0 = empty (bare grounded substrate), state 1 = 5.56 mm square patch with an
              annular slot R = 1.922 mm, r = 1.179 mm.
  walk     -- the Phase 15 Part E walk-through pair (docs/history/PHASE15_POLARIZATION_VORTEX.md):
              8 mm cell, square patches L = 5 mm and L = 7 mm.

Writes, per state, a single-cell STEP (same convention as the board packages: zero-thickness
sheets, src.export.array_export.write_hfss_step on a 1x1 grid) and one MoM reference CSV with S11
(xx, yy) at the 5 mm library plane AND the patch plane, for two mesh densities, plus a pair
separation CSV. No solver change: every number comes from src.sparameters.s_parameters_from_config
(the empty state from the analytic stack, exactly as examples/make_paper_board_26ghz.py does).

Usage (project root):  python scripts/adhoc/one_bit_pair_cells.py <out_dir> [f_min_GHz f_max_GHz step_GHz]
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from src.config import config_from_dict  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.green_function import LayeredGreenFunction  # noqa: E402
from src.hologram.two_state_cells import LIBRARY_REFERENCE_OFFSET_M, patch_plane_gamma  # noqa: E402
from src.materials import Substrate  # noqa: E402
from src.phase import principal_phase_deg, wrapped_difference_deg  # noqa: E402
from src.sparameters import s_parameters_from_config, scattering_from_fields  # noqa: E402

EPS_R, TAN_D, THICK_MM = 2.2, 0.0009, 1.575
DESIGN_F_HZ = 10.54e9

# (case id, pair, state, period_mm, geometry dict in mm or None for the empty cell, meshes)
CASES = (
    ("coded4_state0_empty", "coded4", 0, 6.4, None, ()),
    (
        "coded4_state1_slotted",
        "coded4",
        1,
        6.4,
        {"type": "patch_annular_slot", "a": 5.56, "b": 5.56, "R": 1.922, "r": 1.179},
        (32, 48),
    ),
    ("walk_state_L5mm", "walk", 0, 8.0, {"type": "rectangle", "length": 5.0, "width": 5.0}, (32, 64)),
    ("walk_state_L7mm", "walk", 1, 8.0, {"type": "rectangle", "length": 7.0, "width": 7.0}, (32, 64)),
)


def _base_config(period_mm, geometry_mm, mesh_cells, f_hz):
    return {
        "period_x": period_mm * 1e-3,
        "period_y": period_mm * 1e-3,
        "substrate": {"name": "user", "eps_r": EPS_R, "tan_delta": TAN_D, "thickness": THICK_MM * 1e-3, "backing": "pec"},
        "geometry": {k: (v * 1e-3 if k != "type" else v) for k, v in geometry_mm.items()},
        "excitation": {"frequency": f_hz},
        "mesh": {"cells_x": mesh_cells, "cells_y": mesh_cells, "floquet_order": mesh_cells},
        "reference_planes": {"top_air_thickness": LIBRARY_REFERENCE_OFFSET_M},
    }


def s11_top(case, mesh_cells, f_hz):
    """(S11_xx, S11_yy) referenced 5 mm above the patch, as the auto-built libraries store it."""
    _, _, _, period_mm, geom, _ = case
    if geom is None:
        sub = Substrate("user", EPS_R, TAN_D, THICK_MM * 1e-3)
        green = LayeredGreenFunction(sub, f_hz, backing="pec")
        planes = ReferencePlanes(
            top_air_thickness=LIBRARY_REFERENCE_OFFSET_M, bottom_air_thickness=None, substrate_thickness=THICK_MM * 1e-3
        )
        s = scattering_from_fields(green, planes, np.zeros((1, 2, 2), dtype=complex)).s11
    else:
        s = s_parameters_from_config(config_from_dict(_base_config(period_mm, geom, mesh_cells, f_hz))).s11
    return complex(s[0, 0]), complex(s[1, 1])


def frequencies(f_min_ghz=10.0, f_max_ghz=11.1, step_ghz=0.05):
    n = int(round((f_max_ghz - f_min_ghz) / step_ghz)) + 1
    f = [round(f_min_ghz + i * step_ghz, 6) * 1e9 for i in range(n)]
    if not any(abs(x - DESIGN_F_HZ) < 1.0 for x in f):
        f.append(DESIGN_F_HZ)
    return sorted(f)


def solve_all(freqs_hz, log=print):
    rows = []
    for case in CASES:
        meshes = case[5] or (0,)  # 0 = analytic (no mesh)
        for m in meshes:
            t0 = time.time()
            for f in freqs_hz:
                gxx, gyy = s11_top(case, m, f)
                pxx, pyy = patch_plane_gamma(gxx, f, LIBRARY_REFERENCE_OFFSET_M), patch_plane_gamma(gyy, f, LIBRARY_REFERENCE_OFFSET_M)
                rows.append(
                    {
                        "case": case[0],
                        "pair": case[1],
                        "state": case[2],
                        "period_mm": case[3],
                        "mesh_cells": m if m else "analytic",
                        "frequency_GHz": round(f / 1e9, 6),
                        "S11_xx_mag": abs(gxx),
                        "S11_xx_phase_deg_ref5mm": float(principal_phase_deg(gxx)),
                        "S11_xx_phase_deg_ref0mm": float(principal_phase_deg(pxx)),
                        "S11_yy_mag": abs(gyy),
                        "S11_yy_phase_deg_ref5mm": float(principal_phase_deg(gyy)),
                        "S11_yy_phase_deg_ref0mm": float(principal_phase_deg(pyy)),
                    }
                )
            log(f"{case[0]} mesh={m or 'analytic'}: {len(freqs_hz)} freqs in {time.time() - t0:.0f} s")
    return rows


def pair_separation(rows):
    """Per pair, per mesh pairing (finest-with-finest, coarsest-with-coarsest), per frequency: the
    patch-plane xx phase separation (wrapped to [0, 180]) and the |S11| ratio in dB."""
    out = []
    for pair in ("coded4", "walk"):
        sel = [r for r in rows if r["pair"] == pair]
        meshes = sorted({r["mesh_cells"] for r in sel if r["mesh_cells"] != "analytic"})
        for level, m in (("coarse", meshes[0]), ("fine", meshes[-1])):
            s0 = {r["frequency_GHz"]: r for r in sel if r["state"] == 0 and r["mesh_cells"] in (m, "analytic")}
            s1 = {r["frequency_GHz"]: r for r in sel if r["state"] == 1 and r["mesh_cells"] in (m, "analytic")}
            for f in sorted(s0):
                a, b = s0[f], s1[f]
                sep = abs(float(wrapped_difference_deg(a["S11_xx_phase_deg_ref0mm"], b["S11_xx_phase_deg_ref0mm"])))
                out.append(
                    {
                        "pair": pair,
                        "mesh_level": level,
                        "mesh_cells": m,
                        "frequency_GHz": f,
                        "state0_phase_deg_ref0mm": a["S11_xx_phase_deg_ref0mm"],
                        "state1_phase_deg_ref0mm": b["S11_xx_phase_deg_ref0mm"],
                        "separation_deg": sep,
                        "mag_ratio_dB": float(20 * np.log10(b["S11_xx_mag"] / a["S11_xx_mag"])),
                    }
                )
    return out


def write_csv(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})


def write_steps(out_dir):
    from src.export.array_export import write_hfss_step

    for case_id, _, _, period_mm, geom, _ in CASES:
        if geom is None:
            cell = (0.0, 0.0)
        elif geom["type"] == "rectangle":
            cell = (geom["length"] * 1e-3, geom["width"] * 1e-3)
        else:
            cell = {"kind": "annular_slot", "a_m": geom["a"] * 1e-3, "b_m": geom["b"] * 1e-3, "R_m": geom["R"] * 1e-3, "r_m": geom["r"] * 1e-3}
        write_hfss_step([[cell]], period_mm * 1e-3, period_mm * 1e-3, THICK_MM * 1e-3, out_dir / f"{case_id}.step")


def main(argv):
    out = Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    fargs = [float(x) for x in argv[2:5]] if len(argv) >= 5 else []
    freqs = frequencies(*fargs)
    write_steps(out)
    rows = solve_all(freqs)
    write_csv(rows, out / "mom_reference_one_bit_cells.csv")
    write_csv(pair_separation(rows), out / "mom_pair_separation.csv")


if __name__ == "__main__":
    main(sys.argv)
