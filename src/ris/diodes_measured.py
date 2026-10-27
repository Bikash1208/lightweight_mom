"""MEASURED_VENDOR_SPARAM diode: the MACOM MA4AGP907 (= MADP-000907-14020) from its vendor-measured S-parameters (Phase 28, additive module).

The series impedance Z(f) = Z0 (S11 + S22)/S21 was extracted from the measured 2-port files by `scripts/research/phase28_diode_extract.py`; only the extracted values (4 - 20 GHz) are committed
(`results/phase28/diodes/ma4agp907_series_z.json`, with file names, sha256, bias and file date), the vendor .s2p files are not. States: ON at 5 / 10 / 15 mA, OFF at 0 V / -5 V.
`z_of_f(state)` returns a callable f [Hz] -> Z [ohm] (linear interpolation of the table, passive: the small negative real part of the OFF extraction is measurement noise and is set to zero);
`load(...)` is the `LumpedLoad` of a diode footprint (frequency-dependent load, `transverse_conducts=False`). `scale` = {"R": s} scales the real part of an ON state, {"C": s} scales the OFF capacitance (tolerance studies).
"""

import json
from pathlib import Path

import numpy as np

from ..loads import LumpedLoad

DATA = Path(__file__).resolve().parents[2] / "results/phase28/diodes/ma4agp907_series_z.json"
LABEL = "MEASURED_VENDOR_SPARAM"


def _table():
    return json.loads(DATA.read_text(encoding="utf-8"))["states"]


def provenance(state):
    e = _table()[state]
    return {"label": LABEL, "file": e["file"], "sha256": e["sha256"], "file_date": e["file_date"], "bias": e["bias_header"]}


def z_of_f(state, scale=None):
    e = _table()[state]
    f = np.array(e["f_ghz"]) * 1e9
    re, im = np.array(e["z_re_ohm"]), np.array(e["z_im_ohm"])
    sc = scale or {}
    if state.startswith("ON"):
        re = re * sc.get("R", 1.0)
    else:
        re = np.zeros_like(re)  # passive: the small negative OFF real part is noise
        im = im / sc.get("C", 1.0)  # C -> s C: Im Z = -1/(w C)

    def z(freq):
        fr = np.asarray(freq, dtype=float)
        return np.interp(fr, f, re) + 1j * np.interp(fr, f, im)

    return z


def load(state, x, y, length, width, direction="x", scale=None):
    return LumpedLoad(x=x, y=y, length=length, width=width, direction=direction, z_of_f=z_of_f(state, scale), transverse_conducts=False)
