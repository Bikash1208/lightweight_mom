"""Phase 26 (launch fraction): TM0 launch fraction of the user's coax pin and the item-A bracket re-run with it.

    python scripts/research/phase26_launch_fraction.py <out.json>

The pin (HFSS variables as given by the user: PinLength 110 mil, PinRadius 12 mil, through RT/duroid 5880
1.575 mm) is evaluated under both readings of "PinLength": measured from the ground plane (the convention of
`mono_h` in the exported AEDT scripts) and measured above the substrate surface.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.dispersion.pin_launch import Launcher, PinLaunch  # noqa: E402

MIL = 25.4e-6
H, EPS_R = 1.575e-3, 2.2
PIN_L, PIN_R = 110 * MIL, 12 * MIL


def split(f, top, a=0.0, n_seg=24):
    pl = PinLaunch(f, EPS_R, H, Launcher(top, PIN_R, a), n_seg=n_seg)
    r = pl.split()
    return {k: float(v) for k, v in r.items()}


def main(out):
    res = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "MATCH_ESTIMATE_THIN_WIRE_TL"], "pin": {}}
    f = 10.5e9
    readings = {
        "PINLENGTH_FROM_GROUND": PIN_L - H,
        "PINLENGTH_ABOVE_SURFACE": PIN_L,
    }
    for name, top in readings.items():
        r = split(f, top)
        r["top_above_surface_mm"] = top * 1e3
        res["pin"][name] = r
    res["frequency_ghz"] = f / 1e9
    # launch fraction vs protrusion (no disk)
    res["vs_protrusion_mm"] = {
        f"{t * 1e3:.2f}": split(f, t, n_seg=16)
        for t in (-0.8e-3, 0.0, 0.6e-3, 1.2e-3, 2.0e-3, 3.0e-3, 4.4e-3, 6.0e-3, 8.0e-3, 10.0e-3)
    }
    # item-A bracket with the launch fraction
    a = json.loads((ROOT / "results" / "phase26" / "hex150_leakage.json").read_text())
    hfss = a["hfss"]["peak_realized_gain_rhcp_dbi"]
    br = {}
    for rname, r in res["pin"].items():
        lf = r["launch_fraction"]
        rows = {}
        for case, c in a["cases"].items():
            g = c["gain_upper_bound_dbi"] + 10 * np.log10(lf)
            rows[case] = {"gain_dbi": float(g), "miss_vs_hfss_db": float(g - hfss)}
        br[rname] = {"launch_fraction": lf, "cases": rows}
    res["bracket_with_launch"] = br
    res["hfss_realized_gain_dbi"] = hfss
    Path(out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res["pin"], indent=1))
    print(json.dumps(br, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
