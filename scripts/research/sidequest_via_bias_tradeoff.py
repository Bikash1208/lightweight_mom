"""Bias selection: R_on(I) of MA4AGP907 from the vendor curves, and the ballast / bus-voltage / RF-loss trade-off of the via GB cell (resistive feed).  MODEL_ONLY.
R_on(I): MACOM V4 sheet, insertion loss of a series diode (50 ohm line) at 10 GHz read from the plot: 5 mA -0.41 dB, 15 mA -0.33 dB, 50 mA -0.31 dB -> R_eff = 100 (10^(IL/20) - 1) = 4.8 / 3.9 / 3.6 ohm
(spec: 4.2 ohm at 10 mA, 10 GHz).  Fit R = R0 + A / I  (R0 3.45 ohm, A 6.75 ohm mA) reproduces all four (4.1 at 10 mA).  1 mA is OUTSIDE the vendor data (extrapolation).  The user's measured model (2.44 ohm at 10 mA) is 0.58 of the spec: scaled fit 2.0 + 3.9 / I.
    python sidequest_via_bias_tradeoff.py"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_via_hw as W
import sidequest_via_cell as V

RES = ROOT / "results/sidequest_hybrid"
VF10, VF10_MAX = 1.33, 1.45
NVT_LN10 = 0.105  # n V_T ln(10) for n about 1.75 (AlGaAs), V: V_F(1 mA) = V_F(10 mA) - 0.105 - (R_s drop difference)


def r_spec(i_mA):
    return 3.45 + 6.75 / i_mA


def r_meas(i_mA):
    return 0.58 * r_spec(i_mA)


def vf(i_mA, vf10=VF10):
    return vf10 - NVT_LN10 * math.log10(10.0 / i_mA) - 0.0042 * (10 - i_mA) * 0.0 - (4.2e-3 * 10 - 4.2e-3 * i_mA)


def main():
    P = W.patch()
    print("R_on(I): spec-scaled fit / measured-scaled: " + ", ".join("%.1f mA: %.1f / %.1f ohm" % (i, r_spec(i), r_meas(i)) for i in (0.5, 1, 2, 5, 10)))
    print("V_F(I) estimate (typ 1.33 V at 10 mA): " + ", ".join("%.1f mA: %.2f V" % (i, vf(i)) for i in (0.5, 1, 2, 5, 10)))
    out = []
    for Rb in (1000.0, 1500.0, 2200.0, 3300.0, 4700.0):
        rb = dict(W.RB, R_b=Rb, C_p=8e-15, C_A=10e-15, cathode="via")
        cs = W.knob2(P, rb=rb)
        if not cs:
            out.append({"R_b": Rb, "feasible": False})
            print("R_b %.0f: no C_s gives beta_mean 1.062" % Rb)
            continue
        C = cs[0]
        yon, yoff = W.branch2("ON", V.F0, C, rb), W.branch2("OFF", V.F0, C, rb)
        ex = {st: [x for x in W.extrema(C, st, rb) if abs(x - 10.5) < 1.5] for st in ("ON", "OFF")}
        r_off = V.root(P, yoff)
        row = {"R_b": Rb, "feasible": True, "C_s_pF": C * 1e12, "step_mS": (yoff - yon).imag * 1e3, "Re_Y_OFF_mS": yoff.real * 1e3, "alpha_OFF_Np_m": r_off[1] if r_off else None, "extrema_within_1p5GHz": {k: [float(x) for x in v] for k, v in ex.items()}, "by_current": {}}
        for I in (0.5, 1.0, 2.0, 5.0):
            for tag, rf in (("spec", r_spec(I)), ("meas", r_meas(I))):
                yon2 = W.branch2("ON", V.F0, C, rb, dR=rf / 2.44)
                ron = V.root(P, yon2)
                row["by_current"][f"{I:g}mA_{tag}"] = {"R_on": rf, "Re_Y_ON_mS": yon2.real * 1e3, "alpha_ON_Np_m": ron[1] if ron else None, "dstep_pct": ((yoff - yon2).imag / (yoff - yon).imag - 1) * 100}
        # board-level extra loss: 150 GB cells of 1430 carry the Re Y; mean-sheet alpha with and without
        ymean = 0.5 * (yon + yoff)
        r0 = V.root(P, complex(0.0, ymean.imag))
        r1 = V.root(P, complex(0.5 * (yon.real + yoff.real) * 150 / 1430, ymean.imag))
        row["board_extra_alpha_Np_m"] = (r1[1] - r0[1]) if (r0 and r1) else None
        row["board_extra_loss_pct_150mm_path"] = (1 - math.exp(-2 * row["board_extra_alpha_Np_m"] * 0.10)) * 100 if row["board_extra_alpha_Np_m"] is not None else None
        out.append(row)
        print("R_b %4.0f ohm: C_s %.3f pF step %.1f mS Re Y_OFF %.2f mS alpha_OFF %.2f Np/m | extrema within 1.5 GHz %s | board extra alpha %.3f Np/m -> %.1f %% over 100 mm" % (Rb, C * 1e12, row["step_mS"], row["Re_Y_OFF_mS"], row["alpha_OFF_Np_m"], ex, row["board_extra_alpha_Np_m"], row["board_extra_loss_pct_150mm_path"]))
        for k, v in row["by_current"].items():
            if k.startswith(("1mA", "2mA")):
                print("     %-9s R_on %.1f ohm: Re Y_ON %.2f mS alpha_ON %.2f Np/m step change %+.1f %%" % (k, v["R_on"], v["Re_Y_ON_mS"], v["alpha_ON_Np_m"], v["dstep_pct"]))
    (RES / "via_bias_tradeoff.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
    print()
    print("DC (150 diodes, 1 mA / 2 mA): P = parallel bus, S4 = strings of 4 cells")
    for Rb in (1500.0, 2200.0, 4700.0):
        for I in (1.0, 2.0):
            v = vf(I)
            vp = v + I * 1e-3 * Rb
            vs = 4 * vp
            print("R_b %4.0f ohm, I %.0f mA: V_F %.2f V -> P bus %.2f V, %.0f mA total, %.2f W | S4 line %.1f V, %.0f mA total, %.2f W | resistor %.1f mW | I spread for V_F +-0.12 V: P %.2f mA, S4 %.2f mA" % (
                Rb, I, v, vp, 150 * I, 150 * I * vp / 1e3 * 1, vs, math.ceil(150 / 4) * I, math.ceil(150 / 4) * I * vs / 1e3, I * I * Rb * 1e-3, 0.12 / Rb * 1e3, 0.48 / (4 * Rb) * 1e3))


if __name__ == "__main__":
    main()
