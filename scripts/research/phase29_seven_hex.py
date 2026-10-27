"""Phase 29: seven-hexagon cluster, tiles alone + pair excitation + tilt/size limits (MODEL_ONLY). Writes results/phase29/seven_hex.json.
Usage: python scripts/research/phase29_seven_hex.py [--sizes 150 75] [--tilts 45 30 20 15 10]"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

from src.hologram import seven_hex_cluster as c
from src.hologram import tensor_cp_design as t
from src.hologram.tensor_leaky_cp import SupercellTable

OUT = Path(__file__).resolve().parents[2] / "results/phase29/seven_hex.json"


def hpbw_deg(g_db, theta, phi):
    """Equivalent circular -3 dB beamwidth from the solid angle of the -3 dB region about the peak (sin-theta weighted)."""
    tg, _ = np.meshgrid(theta, phi, indexing="ij")
    w = np.sin(np.radians(np.maximum(tg, 0.25))) * np.radians(theta[1] - theta[0]) * np.radians(phi[1] - phi[0])
    om = np.sum(w[g_db >= g_db.max() - 3.0])
    return float(np.degrees(2 * np.sqrt(om / np.pi)))


def run_config(table, green, r_mm, tilt, theta, phi, log):
    cen = c.tile_centres_mm(r_mm)
    own = c.owned_cells_mm(cen, r_mm)
    tiles, fields = [], []
    for j, (x, y, az) in enumerate(cen):
        beam = (0.0, 0.0) if j == 0 else (float(tilt), az)
        tile = c.design_tile(table, (x, y), r_mm, *beam, cells=own[j])
        tiles.append(tile)
        fields.append(c.tile_field(tile, green, theta, phi))
    alone = []
    for j, f in enumerate(fields):
        g = c.rhcp_gain_dbi(*f)
        lb = c.lobes(g, theta, phi, within_db=3.0)
        mt = t.pol_metrics(f[0], f[1], theta, phi, tiles[j]["beam"])
        alone.append(
            {
                "tile": j,
                "n_cells": int(tiles[j]["x_mm"].size),
                "eta_conv": tiles[j]["eta_conv"],
                "beam_design_theta_phi": list(tiles[j]["beam"]),
                "peak_gain_dbi": lb[0][0],
                "peak_theta_phi_deg": [lb[0][1], lb[0][2]],
                "hpbw_deg": hpbw_deg(g, theta, phi),
                "rhcp_sll_db": mt["lobes_rel_db_theta_phi_deg_from_peak"][0][0]
                if mt["lobes_rel_db_theta_phi_deg_from_peak"]
                else None,
                "ar_peak_db": mt["ar_peak_db"],
                "ar_worst_minus3db_db": mt["ar_worst_in_minus3db_beam_db"],
            }
        )
    pairs = {}
    names = [(0, k, f"C+O{k}") for k in range(1, 7)] + [(k, k % 6 + 1, f"O{k}+O{k % 6 + 1}") for k in range(1, 7)]
    for a, b, name in names:
        rows = c.pair_scan(fields[a], fields[b], theta, phi, tiles[a]["beam"], tiles[b]["beam"])
        sep = c.angle_between(tiles[a]["beam"], tiles[b]["beam"])
        pairs[name] = {
            "beam_separation_deg": sep,
            "centre_distance_mm": float(np.hypot(cen[a][0] - cen[b][0], cen[a][1] - cen[b][1])),
            "states": sorted({r["state"] for r in rows}),
            "n_phases_single_in_between": sum(r["state"] == "SINGLE_BEAM_IN_BETWEEN" for r in rows),
            "best_peak_gain_dbi": max(r["peak_gain_dbi"] for r in rows),
            "min_n_lobes_3db": min(r["n_lobes_3db"] for r in rows),
            "sweep": rows,
        }
    log(f"R {r_mm} tilt {tilt}: " + ", ".join(f"{k}:{v['n_phases_single_in_between']}/12" for k, v in pairs.items()))
    return {
        "tile_radius_mm": r_mm,
        "outer_tilt_deg": tilt,
        "centres_mm": [list(map(float, q[:2])) for q in cen],
        "tiles_alone": alone,
        "pairs": pairs,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=float, nargs="+", default=[150.0, 75.0])
    ap.add_argument("--tilts", type=float, nargs="+", default=[45.0, 30.0, 20.0, 15.0, 10.0])
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()
    table = SupercellTable.from_files(t.TABLE_GLOB)
    green = c.green_function()
    theta, phi = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 2.0)
    t0 = time.time()
    res = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "SOLVER_SUPERCELL", "INTER_TILE_COUPLING_NOT_MODELLED"],
        "frequency_ghz": t.F / 1e9,
        "grid": {"theta_step": 1.0, "phi_step": 2.0},
        "gain_anchor": "single-hexagon HFSS run (OFFSET_DB); tilts/sizes other than 30 deg / 150 mm are model extrapolations",
        "configs": [],
    }
    for r in a.sizes:
        for tilt in a.tilts:
            res["configs"].append(
                run_config(
                    table, green, r, tilt, theta, phi, lambda s: print(f"[{time.time() - t0:5.0f}s] {s}", flush=True)
                )
            )
    res["seconds"] = time.time() - t0
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
