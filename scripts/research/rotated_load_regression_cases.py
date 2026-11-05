"""Reference cases for the rotated-load regression (fcc6691): pre-existing behaviour that must be bit-identical before/after the fix.

    python scripts/research/rotated_load_regression_cases.py OUT.json     (run from the project root of the checkout under test)
Values are stored as sha256 fingerprints of the exact complex128 bytes."""

import hashlib
import json
import sys
import warnings
from dataclasses import replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.basis_functions import MeshSettings, RooftopBasis  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.loads import load_matrix  # noqa: E402
from src.ris.cells import OFF, ON, Cp2Cell, ideal_kit, madp_kit  # noqa: E402
from src.ris.diodes import DiodeModel, diode_load, hfss_chain, smp1345_on  # noqa: E402
from src.sparameters import compute_s_parameters_oblique  # noqa: E402

W = 2 * np.pi * 10.5e9
PLANES = None  # set per cell below


def hx(a):
    """Exact fingerprint: sha256 of the complex128 bytes + shape (bitwise comparison without storing the arrays)."""
    a = np.ascontiguousarray(np.asarray(a.toarray() if hasattr(a, "toarray") else a, dtype=np.complex128))
    return {"shape": list(a.shape), "sha256": hashlib.sha256(a.tobytes()).hexdigest(), "norm": float(np.linalg.norm(a))}


def compute():
    out = {}
    warnings.simplefilter("ignore")
    cell = Cp2Cell()
    uc = cell.unit_cell((OFF, ON, OFF, OFF), ideal_kit())
    b = RooftopBasis(uc, 40, 40)
    cases = {
        "x": diode_load(smp1345_on(), 0.0, 0.0, 1e-3, 0.5e-3, "x"),
        "y": diode_load(smp1345_on(), 0.0, 0.0, 1e-3, 0.5e-3, "y"),
        "d0": diode_load(smp1345_on(), 0.0, 0.0, 1e-3, 0.5e-3, 0),
        "d90": diode_load(smp1345_on(), 0.0, 0.0, 1e-3, 0.5e-3, 90),
        "d30_explicit_blocked": diode_load(smp1345_on(), 0.0, 0.0, 1e-3, 0.5e-3, 30.0, transverse_conducts=False),
        "d30_explicit_conducting": diode_load(smp1345_on(), 0.0, 0.0, 1e-3, 0.5e-3, 30.0, transverse_conducts=True),
        "d30_default": diode_load(smp1345_on(), 0.0, 0.0, 1e-3, 0.5e-3, 30.0),
    }
    for k, ld in cases.items():
        out["load_matrix_" + k] = hx(load_matrix(b, [ld], W))
    chain, _ = hfss_chain(smp1345_on(), 0.0, 0.0, 0.25e-3, 0.5e-3, 0.5e-3)
    out["load_matrix_hfss_chain_x"] = hx(load_matrix(b, list(chain), W))
    chain30, _ = hfss_chain(smp1345_on(), 0.0, 0.0, 0.25e-3, 0.5e-3, 0.5e-3, 30.0)
    out["load_matrix_hfss_chain_30_default"] = hx(load_matrix(b, list(chain30), W))
    for name, kit in (("ideal", ideal_kit()), ("madp", madp_kit())):
        for st in ((OFF, OFF, OFF, OFF), (OFF, ON, OFF, OFF), (ON, ON, OFF, OFF)):
            u = cell.unit_cell(st, kit)
            s = compute_s_parameters_oblique(u, 10.5e9, MeshSettings(40, 40, 20), ReferencePlanes(0.0, None, u.substrate.thickness), 0.0, 0.0)
            out[f"cp2_{name}_{''.join(str(int(x == ON)) for x in st)}_s11"] = hx(s.s11)
            out[f"cp2_{name}_{''.join(str(int(x == ON)) for x in st)}_s21"] = hx(s.s21)
    return out


def main():
    out = compute()
    Path(sys.argv[1]).write_text(json.dumps(out, indent=0), encoding="utf-8", newline="\n")
    print(len(out), "cases")


if __name__ == "__main__":
    main()
