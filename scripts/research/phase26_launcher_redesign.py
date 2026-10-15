"""Phase 26 item 5: launcher redesign on the LOADED surface -- cut the pin's direct space-wave radiation.

    python scripts/research/phase26_launcher_redesign.py <out.json>

Model: `dispersion.pin_launch.PinLaunch` (vertical filament + optional top-hat disk in the grounded slab, with an
optional uniform shunt sheet admittance `y_sheet_s` = the patterned surface around the feed). Launch fraction =
TM0 share of the radiated power; direct space wave = 1 - launch fraction.

Feed calibration (HFSS_CALIBRATED, single board, not validated physics). The HFSS S11 sweep of
`hfss_results/hex150_tensor_cp_single` (min -27.3 dB at 10.10 GHz, -15.0 dB at 10.54, band 8.8-11.6 GHz) is
matched by the thin-wire model for the user's pin (6.0 mm above the surface, 7.6 mm from the ground; NOT the
110 mil = 2.794 mm of the old README) only with an effective radius of ~1.5 mm (0.635 mm: resonance 5 % low) and an
input resistance multiplied by `R_SCALE` (fit below, ~1.3): the thin-wire radiated-power resistance is about
3/4 of the HFSS one. The redesigns are judged with the same R_SCALE and the same series-reactance model; S11 is an
ESTIMATE (tens of ohms), the HFSS run is the arbiter.
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.dispersion.pin_launch import Launcher, PinLaunch  # noqa: E402

F0, H, EPS_R = 10.54e9, 1.575e-3, 2.2
R_PIN = 1.5e-3  # EFFECTIVE radius of the fit (project pin 0.635 mm; see README: the connector/gap acts as a thicker pin)
TOP_BASE = 6.0e-3  # 6.0 mm above the surface (7.6 mm from the ground), as reported by the user
N_SEG = 10
S11_CSV = ROOT / "hfss_results/hex150_tensor_cp_single/S11_8_12GHz.csv"


def zin(f, top, disk, y_ms, r_pin=R_PIN):
    pl = PinLaunch(f, EPS_R, H, Launcher(top, r_pin, disk), n_seg=N_SEG, y_sheet_s=1j * y_ms * 1e-3)
    r = pl.split()
    return complex(r["r_rad_ohm"], r["x_est_ohm"]), float(r["launch_fraction"])


def s11_db(z, r_scale=1.0):
    z = complex(z.real * r_scale, z.imag)
    return 20 * np.log10(max(abs((z - 50.0) / (z + 50.0)), 1e-9))


def calibrate():
    a = np.array([[float(v) for v in r] for r in list(csv.reader(open(S11_CSV)))[1:]])
    f = a[::10, 0] * 1e9
    s = a[::10, 1]
    zs = np.array([zin(fi, TOP_BASE, 0.0, 0.0)[0] for fi in f])

    def cost(c):
        m = 20 * np.log10(np.abs((zs.real * c + 1j * zs.imag - 50) / (zs.real * c + 1j * zs.imag + 50)))
        return float(np.sqrt(np.mean((m - s) ** 2)))

    r = minimize_scalar(cost, bounds=(0.5, 6.0), method="bounded")
    c = float(r.x)
    m = 20 * np.log10(np.abs((zs.real * c + 1j * zs.imag - 50) / (zs.real * c + 1j * zs.imag + 50)))
    return {"r_scale": c, "rms_db": float(r.fun), "rms_db_unscaled": cost(1.0),
            "model_min_db": float(m.min()), "model_f_at_min_ghz": float(f[np.argmin(m)] / 1e9),
            "hfss_min_db": float(s.min()), "hfss_f_at_min_ghz": float(f[np.argmin(s)] / 1e9),
            "model_s11_at_10p54_db": float(np.interp(10.54e9, f, m)), "hfss_s11_at_10p54_db": float(np.interp(10.54e9, f, s))}


def main(out):
    cal = calibrate()
    c = cal["r_scale"]
    print("calibration", cal, flush=True)
    base_z, base_lf = zin(F0, TOP_BASE, 0.0, 0.0)
    res = {"labels": ["MODEL_ONLY", "HFSS_CALIBRATED_FEED_R_SCALE", "MATCH_ESTIMATE_THIN_WIRE_TL", "NOT_HFSS_VALIDATED"],
           "calibration": cal,
           "baseline": {"top_above_surface_mm": TOP_BASE * 1e3, "launch_fraction": base_lf, "z_model_ohm": [base_z.real, base_z.imag],
                        "s11_db_calibrated": s11_db(base_z, c)}, "rows": []}
    hfss_gain = 23.05
    tops = np.arange(1.0, 9.01, 0.5) * 1e-3
    for y_ms in (0.0, 4.0, 8.0, 12.0, 16.0):
        for disk in (0.0, 1.5e-3, 3.0e-3, 4.5e-3):
            # pin top for X = 0 (series resonance) at this (disk, sheet)
            def x_of(t):
                return zin(F0, t, disk, y_ms)[0].imag

            xs = [x_of(t) for t in tops]
            roots = []
            for t0, t1, x0, x1 in zip(tops, tops[1:], xs, xs[1:], strict=False):
                if x0 * x1 < 0 and abs(x0) < 400 and abs(x1) < 400:
                    roots.append(brentq(x_of, t0, t1))
            cand = [("resonant", t) for t in roots] + [("fixed", t) for t in (2.0e-3, 4.0e-3, TOP_BASE)]
            for kind, t in cand:
                z, lf = zin(F0, t, disk, y_ms)
                g = s11_db(z, c)
                res["rows"].append({
                    "kind": kind, "top_mm": t * 1e3, "disk_mm": disk * 1e3, "sheet_mS": y_ms, "launch_fraction": lf,
                    "direct_space_share": 1 - lf, "r_ohm": z.real, "x_ohm": z.imag, "s11_db_calibrated": g,
                    "match_ok": bool(g <= -10.0),
                    "gain_change_vs_baseline_db": float(10 * np.log10(lf / base_lf)),
                    "predicted_gain_dbi_hfss_anchored": hfss_gain + float(10 * np.log10(lf / base_lf)),
                })
        print("sheet", y_ms, "done", flush=True)
    ok = [r for r in res["rows"] if r["match_ok"]]
    res["best_matched"] = max(ok, key=lambda r: r["launch_fraction"]) if ok else None
    res["best_any"] = max(res["rows"], key=lambda r: r["launch_fraction"])
    Path(out).write_text(json.dumps(res, indent=1))
    b = res["best_matched"]
    print("baseline launch %.3f (direct %.3f), S11 est %.1f dB" % (base_lf, 1 - base_lf, res["baseline"]["s11_db_calibrated"]))
    if b:
        print("best matched:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in b.items()})
    print("best any:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in res["best_any"].items()})


if __name__ == "__main__":
    main(sys.argv[1])
