"""I4: probe launcher optimisation with the repo's PinLaunch model: pin top height and top-hat disc radius. Launch fraction = P_sw/(P_sw+P_space); S11 from the calibrated probe
impedance Z = 1.30 R_rad + j (X_est + dX), dX = 17.8 - X_est(6 mm, no disc) (the Phase 26 HFSS calibration of the bare 6 mm probe), 50 ohm line. Whatever needs more than this (anti-pad, via, feed pin
dielectric) is an HFSS item, not a claim. MODEL_ONLY."""
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
warnings.simplefilter("ignore")
import numpy as np

from src.dispersion.pin_launch import Launcher, PinLaunch

F0 = 10.54e9


def one(top_mm, disk_mm, r=0.305):
    s = PinLaunch(F0, 2.2, 1.575e-3, Launcher(top_above_surface_m=top_mm * 1e-3, pin_radius_m=r * 1e-3, disk_radius_m=disk_mm * 1e-3), tan_delta=0.0009).split()
    return s


if __name__ == "__main__":
    ref = one(6.0, 0.0)
    dx = 17.8 - ref["x_est_ohm"]
    rows = []
    for top in (2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0):
        for disk in (0.0, 1.0, 1.5, 2.0, 3.0):
            s = one(top, disk)
            z = 1.30 * s["r_rad_ohm"] + 1j * (s["x_est_ohm"] + dx)
            s11 = 20 * np.log10(abs((z - 50) / (z + 50)))
            rows.append({"top_mm": top, "disk_mm": disk, "launch_fraction": s["launch_fraction"], "R_rad_cal": 1.30 * s["r_rad_ohm"], "X_cal": s["x_est_ohm"] + dx, "S11_db": s11, "beta_k0": s["beta_over_k0"]})
            print(rows[-1], flush=True)
    Path("results/sidequest_hybrid/launcher_opt.json").write_text(json.dumps(rows, indent=1), encoding="utf-8", newline="\n")
