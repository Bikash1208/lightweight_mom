"""Phase 26: model-only aperture-size study of the single-hexagon tensor-CP design (no export).

    python scripts/research/phase26_scaling_study.py <out.json>

Hexagon side (= circumradius) 150, 200, 250, 300 mm, same cell family (SOLVER_SUPERCELL table, nearest grid, mesh-stable
slots), same pin/board/feed at the centre, beam theta 30 / phi 225, 10.54 GHz. For each size the design march is
re-run (eta_design 0.7/0.9/0.99, taper 0.5/0.8; best gain with RHCP SLL <= -13 dB), the far field of the cell currents
is summed, and the figures are anchored to the current HFSS run (150 mm: 23.03 dBi, launch 0.62):

  gain_anchored = D_ap(R) + 10 log(0.62 eta(R)) + offset,   offset = 23.03 - [D_ap(150) + 10 log(0.62 eta(150))]
  main-beam share = k_main * 0.62 eta(R) within8(R),  k_main from the 150 mm HFSS value 0.271
  non-beam power  = 1 - main; grazing+back and LHCP shares = the HFSS split of the non-beam power (0.731, 0.425) -- an
  ASSUMPTION (pin-direct and unconverted edge power given the same polarisation/angle mix) until the E-field line separates them.

The model's idealised aperture and the calibration offset are carried to the larger apertures unchanged (not validated).
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase26_tensor_cp_hexagon as cp  # noqa: E402
from src.hologram.tensor_leaky_cp import SupercellTable, cell_design, march_rays  # noqa: E402

LAUNCH = 0.62
HFSS_GAIN = 23.03
HFSS_MAIN, HFSS_GB, HFSS_LHCP = 0.271, 0.264 + 0.269, 0.310


def run(side_mm, table, green, theta, phi):
    cp.R_HEX = side_mm * 1e-3
    poly = cp.hexagon()
    x, y = cp.lattice_cells(poly)
    best = None
    for eta_d in (0.7, 0.9, 0.99):
        for taper in (0.5, 0.8):
            m = march_rays(table, poly, (0.0, 0.0), cp.K0, cp.THETA0, cp.PHI0, eta_d, taper_c=taper, n_psi=720)
            d = cell_design(m, table, x, y, (0.0, 0.0), cp.K0, cp.THETA0, cp.PHI0, hand=+1)
            et, ep = cp.pattern(d, x, y, green, theta, phi)
            mt = cp.metrics(et, ep, theta, phi)
            rec = {"eta_design": eta_d, "taper_c": taper, "eta_conv": m["eta_conv"],
                   "eta_by_ray": [float(m["eta_conv_by_ray"].min()), float(m["eta_conv_by_ray"].max())],
                   "d_ap_dbi": mt["directivity_ap_rhcp_dbi"], "sll_db": mt["rhcp_sll_db"], "xpd_db": mt["xpd_at_peak_db"],
                   "ar_peak_db": mt["ar_at_peak_db"], "ar_beam_max_db": mt["ar_max_over_minus3db_beam_db"],
                   "within8": mt["within_8deg_fraction"], "within8_rhcp": mt["within_8deg_rhcp_fraction"],
                   "peak_theta_phi": [mt["peak_rhcp_theta_deg"], mt["peak_rhcp_phi_deg"]], "n_cells": int(x.size)}
            g = rec["d_ap_dbi"] + 10 * np.log10(LAUNCH * rec["eta_conv"])
            rec["gain_unanchored_dbi"] = float(g)
            print(f"  side {side_mm}: eta_d {eta_d} taper {taper}: eta {rec['eta_conv']:.3f} D_ap {rec['d_ap_dbi']:.2f} G {g:.2f} SLL {rec['sll_db']:.1f}", flush=True)
            key = (rec["sll_db"] <= -13.0, g)
            if best is None or key > best[0]:
                best = (key, rec)
    return best[1]


def main(out):
    table = SupercellTable.from_files(cp.TABLE_GLOB)
    green = cp.LayeredGreenFunction(cp.SUB, cp.F, backing="pec")
    theta, phi = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 2.0)
    rows = {}
    for side in (150, 200, 250, 300):
        rows[side] = run(side, table, green, theta, phi)
    off = HFSS_GAIN - rows[150]["gain_unanchored_dbi"]
    k_main = HFSS_MAIN / (LAUNCH * rows[150]["eta_conv"] * rows[150]["within8"])
    f_gb = HFSS_GB / (1 - HFSS_MAIN)
    f_lh = HFSS_LHCP / (1 - HFSS_MAIN)
    area = lambda s: 1.5 * np.sqrt(3.0) * s * s  # noqa: E731
    res = {"labels": ["MODEL_ONLY", "HFSS_ANCHORED_TO_CURRENT_RUN", "NOT_HFSS_VALIDATED", "GAMMA_0_COLLINEAR_ONLY"],
           "anchor": {"hfss_gain_dbi": HFSS_GAIN, "launch": LAUNCH, "offset_db": float(off), "k_main": float(k_main),
                      "grazing_back_fraction_of_nonbeam": f_gb, "lhcp_fraction_of_nonbeam": f_lh}, "sizes": {}}
    for side, r in rows.items():
        g = r["gain_unanchored_dbi"] + off
        main = k_main * LAUNCH * r["eta_conv"] * r["within8"]
        non = 1.0 - main
        pin_direct = 1 - LAUNCH
        res["sizes"][str(side)] = {
            **r, "area_m2": float(area(side) * 1e-6), "gain_anchored_dbi": float(g),
            "aperture_efficiency_vs_4piAcos30": float(10 ** ((g - 10 * np.log10(4 * np.pi * area(side) * 1e-6 * np.cos(np.radians(30)) / (299792458.0 / 10.54e9) ** 2)) / 10)),
            "launch": LAUNCH, "pin_direct_share": pin_direct, "unconverted_share": float(LAUNCH * (1 - r["eta_conv"])),
            "main_beam_8deg_share": float(main), "grazing_plus_back_share": float(non * f_gb), "lhcp_share": float(non * f_lh),
            "four_hex_board_mm_extent": [float(2.5 * np.sqrt(3.0) * side), float(3.5 * side)],
        }
        print(side, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in res["sizes"][str(side)].items() if k in ("eta_conv", "d_ap_dbi", "gain_anchored_dbi", "main_beam_8deg_share", "grazing_plus_back_share", "lhcp_share", "xpd_db", "sll_db")})
    Path(out).write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main(sys.argv[1])
