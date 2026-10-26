"""Phase 26: writes `hfss_exports/launcher_calibration/` -- the bare grounded RT/duroid 5880 board (same outline as
`four_port_hex150_gap0`, NO hologram) with ONE feed at port 3's position and a chosen launcher, as native AEDT
build scripts + README with the model's predictions to compare.

    python scripts/research/phase26_launcher_package.py
"""

import json
import sys
from decimal import Decimal
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.dispersion.leakage_model import surface_wave_power_per_h2  # noqa: E402
from src.dispersion.pin_launch import Launcher, PinLaunch, _Layer  # noqa: E402
from src.export.aedt_script import generate_aedt_script  # noqa: E402
from src.export.exact_board import ExactBoard  # noqa: E402
from src.frequency import free_space_wavenumber  # noqa: E402

OUT = ROOT / "hfss_exports" / "launcher_calibration"
F, H, EPS_R, TAND = 10.54e9, 1.575e-3, 2.2, 0.0009
FEED = (Decimal("-64.952"), Decimal("-112.5"))
R_PIN_MM = 0.3048  # 12 mil
VARIANTS = {
    # name: (mono_h from ground [mm], hat radius [mm], description)
    "A_current_pin": (2.794, 0.0, "the user's current pin read as PinLength 110 mil from the ground plane"),
    "B_matched_quarter_wave_pin": (6.575, 0.0, "bare pin 5.0 mm above the substrate: best launched share of accepted power in the model sweep"),
    "C_tall_pin_top_hat": (11.575, 4.0, "pin 10 mm above the substrate with a 4 mm top-hat: best launch fraction (0.82), unmatched in the model"),
}


def bare_board():
    ten = Decimal("326.4"), Decimal("268.8")
    return ExactBoard(
        patches=(),
        feeds_mm=(FEED,),
        port_rects_mm=((-ten[0], ten[0], -ten[1], ten[1]),),
        x_min=-ten[0],
        x_max=ten[0],
        y_min=-ten[1],
        y_max=ten[1],
        thickness_mm=Decimal("1.575"),
        period_mm=Decimal("6.4"),
        max_rounding_mm=0.0,
    )


_FIELD_LINE = '''
# radial line 1 mm above the surface, phi = 225 deg from the feed, rho = 10 .. 140 mm (export Mag_E, E_z along it)
import math
_fx, _fy = %s, %s
_p0 = (_fx + 10.0 * math.cos(math.radians(225.0)), _fy + 10.0 * math.sin(math.radians(225.0)))
_p1 = (_fx + 140.0 * math.cos(math.radians(225.0)), _fy + 140.0 * math.sin(math.radians(225.0)))
oEditor.CreatePolyline(["NAME:PolylineParameters", "IsPolylineCovered:=", False, "IsPolylineClosed:=", False,
    ["NAME:PolylinePoints",
     ["NAME:PLPoint", "X:=", mm(_p0[0]), "Y:=", mm(_p0[1]), "Z:=", H + "+1mm"],
     ["NAME:PLPoint", "X:=", mm(_p1[0]), "Y:=", mm(_p1[1]), "Z:=", H + "+1mm"]],
    ["NAME:PolylineSegments", ["NAME:PLSegment", "SegmentType:=", "Line", "StartIndex:=", 0, "NoOfPoints:=", 2]],
    ["NAME:PolylineXSection", "XSectionType:=", "None", "XSectionOrient:=", "Auto", "XSectionWidth:=", "0mm",
     "XSectionTopWidth:=", "0mm", "XSectionHeight:=", "0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=", "Corner"]],
    ["NAME:Attributes", "Name:=", "RadialLine225", "Flags:=", "NonModel#", "Color:=", "(255 0 0)", "Transparency:=", 0,
     "PartCoordinateSystem:=", "Global", "UDMId:=", "", "MaterialValue:=", '"vacuum"', "SurfaceMaterialValue:=", '""',
     "SolveInside:=", True, "IsMaterialEditable:=", True, "UseMaterialAppearance:=", False, "IsLightweight:=", False])
'''

_HAT = '''
HAT_R = "%smm"      # top-hat disk on the pin (0mm = none); 0.1 mm thick, centred on the pin axis, united with the pin
if HAT_R != "0mm":
    for k, (x, y) in enumerate(FEEDS):
        oEditor.CreateCylinder(["NAME:CylinderParameters", "XCenter:=", mm(x), "YCenter:=", mm(y),
            "ZCenter:=", "mono_h-0.1mm", "Radius:=", HAT_R, "Height:=", "0.1mm", "WhichAxis:=", "Z", "NumSides:=", "0"],
            attrs("Hat%%d" %% (k + 1), "pec", False, "(255 255 0)"))
        unite(["Monopole%%d" %% (k + 1), "Hat%%d" %% (k + 1)])
'''


def build_script(name, mono_h, hat):
    text = generate_aedt_script(bare_board(), "launcher_" + name, EPS_R, TAND, F, R_PIN_MM, "2.05", mono_h)
    a = text.index("# patches, united in blocks")
    b = text.index('oBnd = oDesign.GetModule("BoundarySetup")')
    text = text[:a] + text[b:]
    text = text.replace('oBnd.AssignPerfectE(["NAME:PerfE_Patches", "Objects:=", ["Patches"], "InfGroundPlane:=", False])\n', "")
    marker = "    subtract(\"Ground\", holes, False)"
    text = text.replace(marker, (_HAT % ("%g" % hat)).replace("\n", "\n    ").rstrip(" ") + "\n" + marker, 1) if hat > 0 else text
    text = text.replace("oEditor.FitAll()", (_FIELD_LINE % (FEED[0], FEED[1])).strip("\n") + "\noEditor.FitAll()", 1)
    text = text.replace('"Built: 0 patches', '"Built: bare board, no patches')
    # far field to 180 deg for the power integral
    text = text.replace('"ThetaStart:=", "0deg", "ThetaStop:=", "90deg", "ThetaStep:=", "0.5deg"', '"ThetaStart:=", "0deg", "ThetaStop:=", "180deg", "ThetaStep:=", "1deg"')
    text = text.replace('"PhiStart:=", "0deg", "PhiStop:=", "360deg", "PhiStep:=", "1deg"', '"PhiStart:=", "-180deg", "PhiStop:=", "180deg", "PhiStep:=", "2deg"')
    return text


def complex_beta(eps_r, tan_delta, f, h):
    k0 = free_space_wavenumber(f)
    eps_c = eps_r * (1 - 1j * tan_delta)
    from src.dispersion.analytic_tm0 import tm0_beta

    b = complex(tm0_beta(f, eps_r, h))
    d = lambda kt: _Layer(kt, k0, eps_c, h).d  # noqa: E731
    b0, b1 = b, b * (1 + 1e-6)
    f0, f1 = d(b0), d(b1)
    for _ in range(60):
        b2 = b1 - f1 * (b1 - b0) / (f1 - f0)
        b0, f0, b1 = b1, f1, b2
        f1 = d(b1)
        if abs(b1 - b0) < 1e-12 * abs(b1):
            break
    return b1


def predictions():
    k0 = free_space_wavenumber(F)
    from src.dispersion.analytic_tm0 import tm0_beta

    beta = tm0_beta(F, EPS_R, H)
    kappa = np.sqrt(beta**2 - k0**2)
    w = surface_wave_power_per_h2(beta, EPS_R, H, F)
    eta0 = 376.730313
    bc = complex_beta(EPS_R, TAND, F, H)
    out = {"beta_over_k0": beta / k0, "kappa_np_per_m": kappa, "alpha_dielectric_np_per_m": float(-bc.imag), "variants": {}}
    for name, (mono, hat, _desc) in VARIANTS.items():
        top = (mono - 1.575) * 1e-3
        pl = PinLaunch(F, EPS_R, H, Launcher(top, R_PIN_MM * 1e-3, hat * 1e-3), n_seg=24)
        r = pl.split()
        z = complex(r["r_rad_ohm"], r["x_est_ohm"])
        gam = abs((z - 50) / (z + 50))
        rho = np.array([10.0, 20.0, 50.0, 100.0, 140.0]) * 1e-3
        p_sw = r["launch_fraction"]  # per 1 W ACCEPTED (independent of the match estimate)
        h0 = np.sqrt(p_sw / (2 * np.pi * rho * w)) * np.exp(-(-bc.imag) * rho)
        ez = eta0 * (beta / k0) * h0 * np.exp(-kappa * 1e-3)
        out["variants"][name] = {
            "launch_fraction": r["launch_fraction"],
            "r_rad_ohm": r["r_rad_ohm"],
            "x_est_ohm": r["x_est_ohm"],
            "gamma_est_db": float(20 * np.log10(max(gam, 1e-9))),
            "ez_peak_v_per_m_per_sqrt_w_accepted_at_1mm": {f"{x * 1e3:g}": float(v) for x, v in zip(rho, ez, strict=True)},
        }
    return out


README = """# launcher_calibration — bare grounded board, no hologram (Phase 26)

Purpose: measure what the surface-wave launcher alone does, so the hologram's leakage and the launcher's
launch fraction are no longer entangled in one HFSS number. Same substrate (RT/duroid 5880, eps_r 2.2,
tan d 0.0009, 1.575 mm), same board outline as `four_port_hex150_gap0` (652.8 x 537.6 mm, x -326.4..326.4,
y -268.8..268.8), ONE feed at port 3's position (-64.952, -112.5) mm, **no patches**. 10.54 GHz setup,
+-5 % sweep, far-field sphere theta 0..180 deg (1 deg), phi -180..180 deg (2 deg), as in the main export.

## Which pin is which

| Script | Pin | Why |
|---|---|---|
| `build_A_current_pin.py` | `mono_h` = 2.794 mm from the ground (110 mil), radius 12 mil | your current pin, read as PinLength measured from the ground plane (the `mono_h` convention of the exported scripts). **If your PinLength is measured above the substrate surface, set `mono_h = 4.369mm` in HFSS** (design variable) -- please tell me which it is. |
| `build_B_matched_quarter_wave_pin.py` | `mono_h` = 6.575 mm (5.0 mm above the surface), radius 12 mil | the model's best *launched share of accepted power* (launch fraction x match) |
| `build_C_tall_pin_top_hat.py` | `mono_h` = 11.575 mm with a 4 mm-radius, 0.1 mm top-hat | the model's best launch fraction (0.82) -- estimated unmatched, so this is a physics probe, not a design |

Run in AEDT: Tools > Run Script. Each script builds substrate, ground (with the coax opening), pin, one 50 ohm
lumped port, air box/radiation boundary, setup, sweep, far-field sphere, and a non-model polyline
`RadialLine225` 1 mm above the surface from rho = 10 to 140 mm at phi = 225 deg from the feed.

## What I need back (per script you run; A first, then B)

1. **Antenna parameters** table (port 1 excited, 1 W): radiated / accepted power, efficiencies, peak
   directivity, as in `hfss_results/hex150_gap0/port3_antenna_params_*.txt`.
2. **Far-field CSV**: realized gain RHCP and LHCP (or Etheta/Ephi, whichever is easier) over the full sphere,
   at 10.5 or 10.54 GHz, same format as `port3_realized_gain_RHCP_LHCP_10p5GHz.csv`.
3. **E-field along `RadialLine225`** at the same frequency: Fields calculator > Input > Quantity > E; Complex
   > Mag (for Mag_E) and separately the **complex Ez** (Real and Imag, or Mag and Phase); Input > Geometry >
   Line > `RadialLine225`; Value (eval) at 141+ points; Results > Create Fields Report > Rectangular Plot with
   X = Distance, then File > Export to a CSV. Export both Mag_E and the complex Ez (Re, Im) vs distance.
4. S11 vs frequency (the touchstone `.s1p` or a CSV of dB(S11) and the real/imag input impedance, 8-12 GHz) --
   this is what checks the model's thin-wire match estimate.

## How I will use them

* The surface wave above a grounded slab: E_z falls off as exp(-kappa z) with kappa = %(kappa).1f Np/m
  (bound mode beta/k0 = %(b).4f), so |E_z| at 1 mm gives the field at the surface. The cylindrical wave then
  gives the power through a circle of radius rho: P_sw(rho) = 2 pi rho w |H0|^2 with w = %(w).3f ohm m and
  |H0| = |E_z0| omega eps0 / beta. Fitting |E_z|(rho) to rho^-1/2 exp(-alpha_d rho) separates the launched power
  from the (small, predicted alpha_d = %(ad).3f Np/m) dielectric loss, and `launch fraction = P_sw / accepted power`.
* The far-field integral gives the directly radiated power, a second, independent measurement of the same split.
* The far-field *pattern* of A (donut-like, azimuth-uniform, peak towards the horizon) vs the patterned
  board's HFSS pattern shows how much of the hex150 board's grazing power is the launcher itself.

## Model predictions to compare (MODEL_ONLY; thin-wire match estimate, good to tens of ohms)

%(table)s

|E_z| is a peak amplitude per sqrt(W) of ACCEPTED power (HFSS antenna-parameter table: Accepted Power; scale by sqrt(accepted [W])), 1 mm above the
surface, rho = 10, 20, 50, 100, 140 mm, including the predicted dielectric loss and cylindrical spreading.
A near-field radial line passes close to the pin at rho = 10 mm: the 10 mm value includes non-TM0 near field.
Model launch fraction is the TM0 share of the power radiated by the pin (`src/dispersion/pin_launch.py`).
"""


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    files = []
    for name, (mono, hat, _d) in VARIANTS.items():
        fn = f"build_{name}.py"
        (OUT / fn).write_text(build_script(name, mono, hat))
        files.append(fn)
    pred = predictions()
    rows = ["| Variant | launch fraction | R_rad [ohm] | X est [ohm] | est. |S11| [dB] | |E_z| 20 / 50 / 100 mm [V/m per sqrt(W accepted)] |", "|---|---|---|---|---|---|"]
    for n, v in pred["variants"].items():
        e = v["ez_peak_v_per_m_per_sqrt_w_accepted_at_1mm"]
        rows.append(
            "| %s | %.3f | %.1f | %.0f | %.1f | %.1f / %.1f / %.1f |"
            % (n, v["launch_fraction"], v["r_rad_ohm"], v["x_est_ohm"], v["gamma_est_db"], e["20"], e["50"], e["100"])
        )
    k0 = free_space_wavenumber(F)
    text = README % {
        "kappa": pred["kappa_np_per_m"],
        "b": pred["beta_over_k0"],
        "w": surface_wave_power_per_h2(pred["beta_over_k0"] * k0, EPS_R, H, F),
        "ad": pred["alpha_dielectric_np_per_m"],
        "table": "\n".join(rows),
    }
    (OUT / "README.md").write_text(text)
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1))
    print(json.dumps(pred, indent=1))


if __name__ == "__main__":
    main()
