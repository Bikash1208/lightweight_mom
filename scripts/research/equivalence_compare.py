"""Compare two `equivalence_case.py` outputs (infra brief rule 5): max |delta phase| [deg], max |delta level| [dB], max relative complex difference over all states / angles / TE-TM entries.
Interchangeable if the relative difference is <= 1e-6 (last-digit CPU / BLAS differences); otherwise stop and report.

    python scripts/research/equivalence_compare.py <a.json> <b.json> [--tol 1e-6]
"""

import json
import math
import sys
from pathlib import Path

import numpy as np


def load(path):
    d = json.loads(Path(path).read_text())
    return d["provenance"], {
        (c["state"], c["theta_deg"], c["phi_deg"]): np.array(c["s11"])[..., 0] + 1j * np.array(c["s11"])[..., 1]
        for c in d["cases"]
    }


def compare(a, b):
    worst_ph = worst_db = worst_rel = 0.0
    for k, sa in a.items():
        sb = b[k]
        for i in range(2):
            for j in range(2):
                za, zb = sa[i, j], sb[i, j]
                if abs(za) < 1e-12 and abs(zb) < 1e-12:
                    continue
                if abs(za) > 1e-6:  # skip phases/levels of numerically zero cross terms
                    worst_ph = max(worst_ph, abs(math.degrees(np.angle(zb / za))))
                    worst_db = max(worst_db, abs(20 * math.log10(abs(zb) / abs(za))))
                worst_rel = max(worst_rel, abs(zb - za) / max(abs(za), abs(zb)))
    return {"max_abs_dphase_deg": worst_ph, "max_abs_dlevel_db": worst_db, "max_rel_complex_diff": worst_rel}


def main(pa, pb, tol=1e-6):
    (ma, a), (mb, b) = load(pa), load(pb)
    out = compare(a, b)
    out["machines"] = [ma["machine"], mb["machine"]]
    out["verdict"] = (
        "INTERCHANGEABLE" if out["max_rel_complex_diff"] <= tol else "DIFFERENCE_ABOVE_TOLERANCE_STOP_AND_REPORT"
    )
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    tol = float(sys.argv[sys.argv.index("--tol") + 1]) if "--tol" in sys.argv else 1e-6
    main(sys.argv[1], sys.argv[2], tol)
