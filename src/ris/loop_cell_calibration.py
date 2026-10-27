"""HFSS calibration hook for the Phase 28 loop-cell family (label HFSS_CHECKED; additive module).

Source: ONE HFSS run of the loop cell R4 (user, 2026-10-27; RT/duroid 5880 3.175 mm, 9 mm period, normal incidence, ROT 0, the diode as R only 2.4 ohm / C only 28.2 fF), extracted by
`scripts/research/phase28_hfss_cell_compare.py` into `results/phase28/f10/loop_R4/hfss_calibration.json`. Floquet mode 2 is x (along the diode), mode 1 is y. The state names here are the PHYSICAL ones:
`R` = the closed, lossy diode (2.4 ohm), `C` = the open diode (28.2 fF); in the HFSS project they were labelled the other way round (OFF.csv holds R, ON.csv holds C: the loss and the isotropy identify them).

* `hfss_tensors(f_ghz)`: (Y_R, Y_C), diagonal 2 x 2 sheet admittances [S] of the two states, linear interpolation in frequency (9 - 11 GHz only; outside: ValueError).
* `delta_to_mom(f_ghz)`: HFSS minus MoM (R-only / C-only states, mesh 60) of each component, valid 9.5 - 10.5 GHz.
* `calibrate(y_mom, state, f_ghz)`: the MoM tensor of ANOTHER cell of the loop family plus the R4 deltas: an APPROXIMATE correction (about +4 % of the switching contrast and the mean susceptance in the one case checked); it carries no validation for other geometries.
* `alpha_exact(f_ghz)`: the sheet-model leakage with the HFSS tensors (script value at the nearest tabulated frequency, 9.5 - 10.5 GHz).
"""

import json
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[2] / "results/phase28/f10/loop_R4/hfss_calibration.json"
LABEL = "HFSS_CHECKED"


def _cal():
    return json.loads(DATA.read_text(encoding="utf-8"))


def _interp(f, fs, vals):
    v = np.array([complex(*c) for c in vals])
    return np.interp(f, fs, v.real) + 1j * np.interp(f, fs, v.imag)


def hfss_tensors(f_ghz):
    c = _cal()
    fs = np.array(c["f_ghz"])
    if not (fs.min() - 1e-9 <= f_ghz <= fs.max() + 1e-9):
        raise ValueError(f"{f_ghz} GHz outside the HFSS-checked band {fs.min():g} - {fs.max():g} GHz")
    out = []
    for s in ("Y_R", "Y_C"):
        xx, yy = _interp(f_ghz, fs, c[s]["xx_re_im_S"]), _interp(f_ghz, fs, c[s]["yy_re_im_S"])
        out.append(np.array([[xx, 0.0], [0.0, yy]], dtype=complex))
    return out[0], out[1]


def delta_to_mom(f_ghz):
    d = _cal()["delta_to_mom"]
    fs = np.array(sorted(float(k) for k in d))
    if not (fs.min() - 1e-9 <= f_ghz <= fs.max() + 1e-9):
        raise ValueError(f"{f_ghz} GHz outside the tabulated {fs.min():g} - {fs.max():g} GHz")
    res = {}
    for s in ("R", "C"):
        res[s] = np.diag([_interp(f_ghz, fs, [d[f"{f:g}"][s]["xx"] for f in fs]), _interp(f_ghz, fs, [d[f"{f:g}"][s]["yy"] for f in fs])])
    return res


def calibrate(y_mom, state, f_ghz):
    """`y_mom` + the R4 delta of `state` ('R' or 'C'): approximate (see the module docstring)."""
    return np.asarray(y_mom, dtype=complex) + delta_to_mom(f_ghz)[state]


def alpha_exact(f_ghz):
    v = _cal()["vs_mom"]
    k = min(v, key=lambda s: abs(float(s) - f_ghz))
    return float(v[k]["hfss"]["alpha_exact"])
