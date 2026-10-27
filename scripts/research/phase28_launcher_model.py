"""Phase 28: side (Yagi-type) surface-wave launcher for RT/duroid 5880 3.175 mm at 10.0 GHz: driven pin + reflector post + directors (grounded posts), as in [78]. MODEL_ONLY, NOT_HFSS_VALIDATED.

    python scripts/research/phase28_launcher_model.py <out.json>

Model (ESTIMATES, stated limits). Each pin is a vertical filament with the project's `dispersion.pin_launch.PinLaunch` (grounded slab, thin-wire current): per pin height it gives the radiated resistance split into the
TM0 surface-wave part R_sw and the space-wave part, the thin-wire reactance X and beta. The pins are coupled through the TM0 cylindrical wave only: Z_mn = sqrt(R_sw,m R_sw,n) H0^(2)(beta d) (real part J0, imaginary part -Y0);
the near-field/space-wave mutual terms are NOT modelled (the pins are 0.15 - 0.4 lambda apart: this is the weakest point). The diagonal is R_n + j X_n. Driven pin: V = 1 at the base; the posts are shorted to the ground (V = 0).
The azimuthal TM0 pattern is F(psi) = sum_n sqrt(R_sw,n) I_n exp(+j beta r_n . psi_hat); P_sw = (1/2)(1/2 pi) int |F|^2 dpsi. S11 is referred to 50 ohm at the pin base.
Two feed calibrations bracket the physics: (A) thin wire, r = 0.635 mm (the project pin), no scaling; (B) the Phase 26 HFSS calibration of the 1.575 mm board (effective radius 1.5 mm, radiated resistance x 1.3): it was fitted to ONE
HFSS run on another board and is not validated for 3.175 mm. The posts have the driven pin's radius. The surface around the pins is the BARE slab (a keep-out around the launcher; the loaded surface starts outside it).
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import differential_evolution
from scipy.special import hankel2

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.dispersion.pin_launch import Launcher, PinLaunch  # noqa: E402

F0, EPS_R, H = 10.0e9, 2.2, 3.175e-3
CAL = {"A_thin_wire_r0p635": (0.635e-3, 1.0), "B_phase26_calibrated_r1p5_x1p3": (1.5e-3, 1.3)}
T_GRID = np.arange(1.0e-3, 10.01e-3, 0.5e-3)  # pin top above the substrate surface
BETA_NOTE = "bare-slab TM0 constant"


def tables():
    out = {}
    for name, (r, scale) in CAL.items():
        rows = []
        for t in T_GRID:
            pl = PinLaunch(F0, EPS_R, H, Launcher(float(t), r, 0.0), n_seg=10)
            s = pl.split()
            rows.append((t, 2 * s["p_sw_w_per_a2"] * scale, s["r_rad_ohm"] * scale, s["x_est_ohm"], s["beta_over_k0"]))
        a = np.array(rows)
        out[name] = {
            "t": a[:, 0],
            "r_sw": CubicSpline(a[:, 0], a[:, 1]),
            "r_tot": CubicSpline(a[:, 0], a[:, 2]),
            "x": CubicSpline(a[:, 0], a[:, 3]),
            "beta_over_k0": float(a[0, 4]),
            "raw": a.tolist(),
        }
    return out


def yagi(params, tab, f=F0):
    """params: dict(Td, reflector=(s, T), directors=[(s_from_previous, T), ...]). Returns S11 [dB], pattern, metrics. Positions along +x (directors), reflector on -x."""
    k0 = 2 * math.pi * f / 299792458.0
    beta = tab["beta_over_k0"] * k0
    xs, ts = [0.0], [params["Td"]]
    xs.append(-params["reflector"][0]); ts.append(params["reflector"][1])
    pos = 0.0
    for s, t in params["directors"]:
        pos += s
        xs.append(pos); ts.append(t)
    n = len(xs)
    rsw = np.array([float(tab["r_sw"](t)) for t in ts])
    z = np.zeros((n, n), dtype=complex)
    for m in range(n):
        for q in range(n):
            if m == q:
                z[m, q] = float(tab["r_tot"](ts[m])) + 1j * float(tab["x"](ts[m]))
            else:
                d = abs(xs[m] - xs[q])
                z[m, q] = math.sqrt(rsw[m] * rsw[q]) * complex(hankel2(0, beta * d))
    e = np.zeros(n, complex)
    e[0] = 1.0
    cur = np.linalg.solve(z, e)
    zin = 1.0 / cur[0] if abs(cur[0]) > 0 else complex("nan")
    s11 = (zin - 50.0) / (zin + 50.0)
    psi = np.radians(np.arange(0.0, 360.0, 1.0))
    amp = np.sqrt(rsw) * cur
    ff = np.array([np.sum(amp * np.exp(1j * beta * np.array(xs) * np.cos(p))) for p in psi])
    pw = np.abs(ff) ** 2
    p_sw = 0.5 * np.mean(pw)
    p_in = 0.5 * zin.real * abs(cur[0]) ** 2
    pk = int(np.argmax(pw))
    half = pw >= 0.5 * pw[pk]
    lo = hi = pk
    while half[(lo - 1) % 360] and (pk - lo) < 180:
        lo -= 1
    while half[(hi + 1) % 360] and (hi - pk) < 180:
        hi += 1
    front = pw[np.abs(((np.degrees(psi) + 180) % 360) - 180) <= 90].sum() / pw.sum()
    return {
        "zin_ohm": [zin.real, zin.imag], "s11_db": float(20 * math.log10(max(abs(s11), 1e-12))), "hpbw_deg": float(hi - lo), "peak_azimuth_deg": float(pk),
        "front_to_back_db": float(10 * math.log10(pw[0] / max(pw[180], 1e-30))), "front_half_power_fraction": float(front), "launch_fraction": float(p_sw / p_in) if p_in > 0 else float("nan"),
        "currents_abs_rel": [float(abs(c / cur[0])) for c in cur], "currents_phase_deg": [float(np.degrees(np.angle(c / cur[0]))) for c in cur], "x_mm": [x * 1e3 for x in xs], "top_mm": [t * 1e3 for t in ts],
        "pattern_db": [float(10 * math.log10(max(p / pw.max(), 1e-6))) for p in pw],
    }


def unpack(v):
    return {"Td": v[0] * 1e-3, "reflector": (v[1] * 1e-3, v[2] * 1e-3), "directors": [(v[3] * 1e-3, v[4] * 1e-3), (v[5] * 1e-3, v[6] * 1e-3)]}


def cost1(v, tab):
    r = yagi(unpack(v), tab)
    c = ((r["hpbw_deg"] - 78.0) / 10.0) ** 2 + max(0.0, r["s11_db"] + 10.0) ** 2 / 4.0 - 0.1 * min(r["front_to_back_db"], 20.0) - 2.0 * r["front_half_power_fraction"] - 4.0 * min(r["launch_fraction"], 0.85)
    return c if math.isfinite(c) else 1e3


def cost(v, tabs):
    """Robust over the two feed calibrations: the worse of the two."""
    return max(cost1(v, t) for t in tabs)


def main(out):
    tabs = tables()
    bounds = [(2.0, 8.0), (3.0, 14.0), (2.0, 9.0), (3.0, 14.0), (1.5, 8.0), (3.0, 14.0), (1.5, 8.0)]  # Td, s_r, T_r, s_1, T_1, s_2, T_2 [mm]
    sol = differential_evolution(cost, bounds, args=(list(tabs.values()),), seed=5, maxiter=80, popsize=16, tol=1e-7, polish=False)
    v = sol.x
    names = ["Td", "s_reflector", "T_reflector", "s_dir1", "T_dir1", "s_dir2", "T_dir2"]
    res = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "ESTIMATE_COUPLED_MONOPOLES"], "frequency_ghz": F0 / 1e9, "substrate": "RT/duroid 5880 eps_r 2.2, 3.175 mm", "note": __doc__,
           "design_mm": dict(zip(names, [float(a) for a in v], strict=True)), "robust_cost": float(sol.fun), "evaluated_under": {}, "freq_sweep": {}}
    for n2 in CAL:
        res["evaluated_under"][n2] = yagi(unpack(v), tabs[n2])
        res["freq_sweep"][n2] = {f"{f / 1e9:g}": {k: yagi(unpack(v), tabs[n2], f)[k] for k in ("s11_db", "hpbw_deg", "front_to_back_db", "launch_fraction")} for f in (9.5e9, 9.75e9, 10e9, 10.25e9, 10.5e9)}
    # the bare driven pin for reference (no posts)
    res["bare_pin"] = {n2: yagi({"Td": v[0] * 1e-3, "reflector": (1.0, 1.0), "directors": []}, tabs[n2]) if False else None for n2 in CAL}
    res["tables"] = {n: tabs[n]["raw"] for n in tabs}
    res["table_columns"] = ["top_above_surface_m", "R_sw_ohm", "R_rad_total_ohm", "X_est_ohm", "beta_over_k0"]
    Path(out).write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    print({k: round(x, 2) for k, x in res["design_mm"].items()}, "cost", round(res["robust_cost"], 2))
    for n2, r in res["evaluated_under"].items():
        print("  under", n2[:2], "zin", [round(x, 1) for x in r["zin_ohm"]], "S11 %.1f HPBW %.0f F/B %.1f front %.2f launch %.2f" % (r["s11_db"], r["hpbw_deg"], r["front_to_back_db"], r["front_half_power_fraction"], r["launch_fraction"]))
        print("     sweep", {f: (round(d["s11_db"], 1), round(d["hpbw_deg"]), round(d["front_to_back_db"], 1)) for f, d in res["freq_sweep"][n2].items()})


if __name__ == "__main__":
    main(sys.argv[1])
