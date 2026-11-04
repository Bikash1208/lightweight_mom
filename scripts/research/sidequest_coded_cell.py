"""1-bit cell pair for the grouped coded hologram: split patch (two plates + one MA4AGP907 across a gap, current along x), measured diode ON 10 mA / OFF 0 V, Phase 26 stack, 10.54 GHz.
Scalar two-level sheet: Y_ON / Y_OFF from the MoM (normal-incidence tensor, mesh 32); mean sheet -> beta; first-order leakage alpha (repo `leakage_alpha_first_order`, c1 = |dY|/pi, duty 0.5,
beam theta 25 deg, zeta 0) and the ohmic absorption alpha_abs of the ON cells (half of the cells) from the complex mean sheet. MODEL_ONLY."""
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_gate_cell as gc
from src.dispersion import sheet_model as sm
from src.dispersion.leakage_model import leakage_alpha_first_order
from src.hologram.coded_cell_eval import sheet_tensor


def cell(a, h, g, wd, state):
    return gc.make_cell(a * 1e-3, h * 1e-3, g * 1e-3, wd * 1e-3, state)


def job(p):
    warnings.simplefilter("ignore")
    a, h, g, wd = p
    out = {"a": a, "h": h, "g": g, "wd": wd}
    for st in ("OFF", "ON"):
        y = sheet_tensor(cell(a, h, g, wd, st), gc.F0, 32)
        out[st] = [complex(y[0, 0]).real * 1e3, complex(y[0, 0]).imag * 1e3, complex(y[1, 1]).real * 1e3, complex(y[1, 1]).imag * 1e3]
    return out


def alpha_pair(y_off, y_on, theta_deg=25.0, zeta_deg=0.0):
    """y_off, y_on: complex sheet admittances [S] along the propagation direction. Returns beta/k0, alpha_beam, alpha_abs [Np/m] for a 50/50 pattern."""
    y_mean = 0.5 * (y_on + y_off)
    mode = sm.find_beta_from_sheet_admittance(1j * y_mean.imag, gc.SUB, gc.F0)
    beta = mode.beta_rad_per_m
    c1 = abs(y_on - y_off) / np.pi
    a_beam = leakage_alpha_first_order(beta, 1j * y_mean.imag, c1, np.radians(theta_deg), np.radians(zeta_deg), gc.SUB, gc.F0)
    # ohmic absorption of the mean conductance: d beta / d b * g_mean (transverse resonance, small loss)
    db = 1e-5
    b2 = sm.find_beta_from_sheet_admittance(1j * (y_mean.imag + db), gc.SUB, gc.F0).beta_rad_per_m
    a_abs = (b2 - beta) / db * y_mean.real
    return beta / gc.K0, a_beam, a_abs


if __name__ == "__main__":
    grid = [(a, h, 0.4, wd) for a in (5.2, 5.6, 6.0) for h in (4.8, 5.6, 6.0) for wd in (0.8, 1.2, 2.0)]
    with ProcessPoolExecutor(4) as ex:
        res = list(ex.map(job, grid))
    rows = []
    for r in res:
        yo = complex(r["OFF"][0], r["OFF"][1]) * 1e-3
        yn = complex(r["ON"][0], r["ON"][1]) * 1e-3
        try:
            bk, ab, aa = alpha_pair(yo, yn)
        except Exception as e:  # noqa: BLE001
            bk, ab, aa = None, None, None
        r.update({"beta_k0": bk, "alpha_beam": ab, "alpha_abs": aa, "dB_mS": (yn.imag - yo.imag) * 1e3, "g_on_mS": yn.real * 1e3})
        rows.append(r)
        print(r["a"], r["h"], r["wd"], "OFF %.2f%+.2fj ON %.2f%+.2fj" % (r["OFF"][0], r["OFF"][1], r["ON"][0], r["ON"][1]), "beta/k0", bk and round(bk, 3), "alpha_beam", ab and round(ab, 2), "alpha_abs", aa and round(aa, 2))
    Path("results/sidequest_hybrid/coded_cell_scan.json").write_text(json.dumps(rows, indent=1), encoding="utf-8", newline="\n")
