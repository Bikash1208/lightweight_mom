"""Generator of tests/data/r1_theta0_reference.json (normal-incidence S-parameters at the Phase 20 R1 base commit, hex floats)."""
import json
import warnings

from src.basis_functions import MeshSettings
from src.excitation import ReferencePlanes
from src.geometry import Geometry, UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.sparameters import compute_s_parameters

CASES = [
    (backing, geo)
    for backing in ("pec", "air")
    for geo in (dict(kind="rectangle", length=6e-3, width=5e-3), dict(kind="rectangle", length=4e-3, width=4e-3, rotation=0.3))
]


def compute():
    warnings.simplefilter("ignore")
    sub = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
    out = {}
    for backing, geo in CASES:
        uc = UnitCell(RectangularLattice(10e-3, 8e-3), Geometry(**geo), sub, Metal(), backing=backing)
        pl = ReferencePlanes(2e-3, None if backing == "pec" else 2e-3, sub.thickness)
        s = compute_s_parameters(uc, 10e9, MeshSettings(16, 16, 8), pl)
        blocks = {"s11": s.s11} if backing == "pec" else {"s11": s.s11, "s12": s.s12, "s21": s.s21, "s22": s.s22}
        out[f"{backing}_{geo['length']}_{geo.get('rotation', 0)}"] = {
            k: [[(float(v.real).hex(), float(v.imag).hex()) for v in row] for row in m] for k, m in blocks.items()
        }
    return out


if __name__ == "__main__":
    json.dump(compute(), open("tests/data/r1_theta0_reference.json", "w"), indent=0)
