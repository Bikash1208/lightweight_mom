"""Four ON-cell variants (and OFF twins for two) bracketing a via-patch resonance near 10.5-11 GHz, from the CALIBRATED model (L_br = 6.0 nH from HFSS ON f(180) = 4.92 GHz). MODEL_ONLY."""
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_mushroom_calibrate as cal
import sidequest_mushroom_cell as mc

FREQS = [4, 6, 8, 9.5, 10.54, 11.5, 12.5, 14, 16, 18]
GEOMS = {f"w{w}_p{p}_g{g}": (w, g, p) for w in (3.2, 4.0, 4.8) for p in (0.8, 2.0) for g in (0.4, 0.8) if w / 2 - (p / 2 + g) >= 0.399}


def job(item):
    name, (w, g, p) = item
    return name, mc.sheet_table(w * 1e-3, FREQS, 32, g * 1e-3, p * 1e-3)


if __name__ == "__main__":
    with ProcessPoolExecutor(4) as ex:
        tabs = dict(ex.map(job, GEOMS.items()))
    Path("results/sidequest_hybrid/mushroom_variant_sheets.json").write_text(json.dumps(tabs), encoding="utf-8", newline="\n")
    import math
    H_MM = 1.575
    def lbr(scn, g, p, r):
        strap = 0.99 * (g / 0.4) * (0.8 / p)
        via = 0.63 * (math.log(2 * H_MM / r) - 0.75) / 2.01
        if scn == "S0":
            return 6.0
        if scn == "S1":
            return 4.38 + via + strap
        return 6.0 * (g / 0.4) * (0.8 / p)  # S2: all of L_br is strap-like
    out = {}
    for name, tab in tabs.items():
        w, g, p = GEOMS[name]
        for r in (0.2, 0.5):
            row = {}
            for scn in ("S0", "S1", "S2"):
                m = cal.ModelPB(tab, lbr(scn, g, p, r) * 1e-9)
                on = {ph: [round(x / 1e9, 2) for x in cal.surface_band(m, "ON", ph)][:2] for ph in (90, 180)}
                row[scn] = {"L_br_nH": round(lbr(scn, g, p, r), 2), "ON": on}
            m = cal.ModelPB(tab, 6e-9)
            row["OFF_beta"] = [round(b / mc.free_space_wavenumber(10.54e9), 3) for b in m.beta_at("OFF", 10.54e9, None) if b / mc.free_space_wavenumber(10.54e9) < 1.4][:1]
            out[f"{name}_r{r}"] = row
            print(name, "r", r, {k: (v["L_br_nH"], v["ON"][180]) if isinstance(v, dict) else v for k, v in row.items()}, flush=True)
    Path("results/sidequest_hybrid/mushroom_variant_scan.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline="\n")
