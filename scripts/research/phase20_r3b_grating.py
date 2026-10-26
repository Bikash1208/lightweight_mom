"""Phase 20 R3b planning: grating-lobe / quantisation-lobe limits of the 15 mm CP-2 pitch at 10.5 GHz. MODEL_ONLY (array-factor arithmetic, isotropic elements, no coupling).

(a) Floquet (Wood) onset: largest plane-wave incidence angle for which no higher Floquet order propagates (margin >= 1 / WOOD_GUARD) vs frequency and azimuth.
(b) Steering limit at normal incidence: largest steering angle with no grating lobe in the principal plane, sin(theta_max) = lambda/P - 1.
(c) Quantisation lobes: a 64-cell row with the ideal phase gradient for steering angle theta_s quantised to 1 or 2 bits (nearest of 0/180 or 0/90/180/270 deg); array factor over -90..90 deg.

    python scripts/research/phase20_r3b_grating.py
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.ris.state_library import C0_MM_GHZ, PERIOD_MM, WOOD_GUARD, higher_order_margin  # noqa: E402

N_CELLS = 64


def onset_angles(f_ghz, phis=(0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0)):
    out = {}
    th = np.arange(0.0, 89.9, 0.1)
    for phi in phis:
        ok = [float(t) for t in th if higher_order_margin(f_ghz, t, phi, PERIOD_MM) >= WOOD_GUARD]
        out[f"{phi:g}"] = max(ok) if ok else 0.0
    return out


def lobes(f_ghz, steer_deg, bits):
    lam = C0_MM_GHZ / f_ghz
    k = 2 * np.pi / lam
    n = np.arange(N_CELLS)
    ideal = -k * PERIOD_MM * n * np.sin(np.radians(steer_deg))
    if bits == 0:
        phase = ideal
    else:
        step = 2 * np.pi / 2**bits
        phase = np.round(ideal / step) * step
    th = np.radians(np.arange(-90.0, 90.0001, 0.05))
    af = np.abs(np.exp(1j * (phase[None, :] + k * PERIOD_MM * n[None, :] * np.sin(th)[:, None])).sum(axis=1)) / N_CELLS
    db = 20 * np.log10(np.maximum(af, 1e-6))
    i0 = int(np.argmax(db))
    main = float(np.degrees(th[i0]))
    # strongest lobe more than 3 beamwidths (~3 deg) away from the main lobe
    far = np.abs(np.degrees(th) - main) > 6.0
    j = int(np.argmax(np.where(far, db, -999)))
    return {
        "steer_deg": steer_deg,
        "bits": bits,
        "main_lobe_deg": main,
        "main_level_db_re_ideal": float(db[i0]),
        "strongest_other_lobe_deg": float(np.degrees(th[j])),
        "strongest_other_lobe_db_re_main": float(db[j] - db[i0]),
    }


def main():
    out = {
        "labels": ["MODEL_ONLY", "ARRAY_FACTOR_ISOTROPIC_ELEMENTS_NO_COUPLING"],
        "period_mm": PERIOD_MM,
        "wood_guard": WOOD_GUARD,
    }
    out["onset_by_frequency_ghz"] = {f"{f:g}": onset_angles(f) for f in (10.0, 10.25, 10.5, 10.75, 11.0)}
    out["grating_free_steering_limit_deg"] = {
        f"{f:g}": float(np.degrees(np.arcsin(max(min(C0_MM_GHZ / f / PERIOD_MM - 1.0, 1.0), -1.0))))
        for f in (10.0, 10.25, 10.5, 10.75, 11.0)
    }
    lam = C0_MM_GHZ / 10.5
    out["at_10p5_ghz"] = {
        "wavelength_mm": lam,
        "period_over_wavelength": PERIOD_MM / lam,
        "phase_step_per_cell_for_60deg_deg": float(np.degrees(2 * np.pi * PERIOD_MM / lam * np.sin(np.radians(60.0)))),
        "cells_per_period_for_60deg": float(lam / PERIOD_MM / np.sin(np.radians(60.0))),
        "cells_per_period_for_45deg": float(lam / PERIOD_MM / np.sin(np.radians(45.0))),
        "cells_per_period_for_30deg": float(lam / PERIOD_MM / np.sin(np.radians(30.0))),
        "max_f_ghz_60deg_no_grating_lobe": float(C0_MM_GHZ / PERIOD_MM / (1.0 + np.sin(np.radians(60.0)))),
        "max_f_ghz_60deg_with_guard": float(C0_MM_GHZ / PERIOD_MM / (WOOD_GUARD + np.sin(np.radians(60.0)))),
    }
    out["quantisation"] = [lobes(10.5, s, b) for s in (20.0, 30.0, 45.0, 60.0) for b in (0, 2, 1)]
    (ROOT / "results/phase20_r3b/grating_limits.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("grating_free_steering_limit_deg", "at_10p5_ghz")}, indent=1))
    print(json.dumps(out["onset_by_frequency_ghz"]["10.5"]))
    for r in out["quantisation"]:
        print(
            f"steer {r['steer_deg']:4.0f} bits {r['bits']} main {r['main_lobe_deg']:6.1f} deg level {r['main_level_db_re_ideal']:6.2f} dB | strongest other {r['strongest_other_lobe_deg']:6.1f} deg {r['strongest_other_lobe_db_re_main']:6.1f} dB"
        )


if __name__ == "__main__":
    main()
