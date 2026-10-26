"""Phase 20 R2 item 3: does the HFSS representation of a series PIN-diode circuit (a chain of single-type lumped boundaries joined by
conductor islands, `src.ris.diodes.hfss_chain`) give the same MoM response as the ideal series model (ONE load of the total Z_L over
the same span)?

    python scripts/research/phase20_r2_hfss_representation.py [out_dir] [--quick]

Cell: R1 demo cell (10 mm, RT5880 1.575 mm, tan d 0.0009, PEC ground, 10 GHz): two 3.5 x 6 mm arms joined by a 0.5 mm wide strip; the
diode span (1.0 mm, along x) sits in the middle of the strip with `pad` mm of strip conductor on each side. Compared at theta = 0, 30
deg (phi = 0), TE and TM, on meshes 40 and 80. Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET.
"""

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.basis_functions import MeshSettings  # noqa: E402
from src.composite_geometry import AxisRect, CompositeGeometry  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.geometry import UnitCell  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402
from src.ris.diodes import diode_load, hfss_chain, smp1345_off, smp1345_on  # noqa: E402
from src.sparameters import compute_s_parameters_oblique_over_loads  # noqa: E402

SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
PLANES = ReferencePlanes(0.0, None, SUB.thickness)
F = 10e9
SPAN = 1.0e-3  # total diode span along x
STRIP_W = 0.5e-3
ELEMENT = 0.25e-3  # length of each lumped boundary in the chain
MODES = ("TE", "TM")


def models():
    return {
        "ON 10mA (Rs 1.5)": smp1345_on("10mA"),
        "ON 1mA (Rs 3.5)": smp1345_on("1mA"),
        "OFF": smp1345_off(),
        "OFF Rp 10k (assumed)": smp1345_off(10e3),
        "OFF Rp 3k (assumed)": smp1345_off(3e3),
    }


def cell(pad):
    half = SPAN / 2 + pad
    arms = (
        AxisRect(-(half + 1.75e-3), 0.0, 3.5e-3, 6e-3),
        AxisRect(half + 1.75e-3, 0.0, 3.5e-3, 6e-3),
        AxisRect(0.0, 0.0, 2 * half, STRIP_W),
    )
    return UnitCell(RectangularLattice(10e-3, 10e-3), CompositeGeometry(arms), SUB, Metal(), backing="pec")


def ideal_and_chain(model):
    ideal = diode_load(model, 0.0, 0.0, SPAN, STRIP_W)
    n = sum([model.R_s > 0, model.C_T is not None, model.L_s > 0])
    island = (SPAN - n * ELEMENT) / (n - 1)
    chain, _ = hfss_chain(model, 0.0, 0.0, ELEMENT, island, STRIP_W)
    return ideal, chain


def compare(pad, mesh_n, thetas=(0.0, 30.0)):
    out = []
    ms = models()
    sets = []
    for m in ms.values():
        i, c = ideal_and_chain(m)
        sets += [(i,), tuple(c)]
    base = cell(pad)
    for th in thetas:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = compute_s_parameters_oblique_over_loads(
                base, F, MeshSettings(mesh_n, mesh_n, max(24, mesh_n // 2)), PLANES, math.radians(th), 0.0, sets
            )
        for k, name in enumerate(ms):
            a, b = res[2 * k].s11, res[2 * k + 1].s11
            for j, mode in enumerate(MODES):
                out.append(
                    dict(
                        model=name,
                        theta_deg=th,
                        mode=mode,
                        mesh=mesh_n,
                        pad_mm=pad * 1e3,
                        ideal_abs=float(abs(a[j, j])),
                        ideal_phase_deg=float(np.degrees(np.angle(a[j, j]))),
                        chain_abs=float(abs(b[j, j])),
                        chain_phase_deg=float(np.degrees(np.angle(b[j, j]))),
                        d_abs=float(abs(b[j, j]) - abs(a[j, j])),
                        d_phase_deg=float(np.degrees(np.angle(b[j, j] / a[j, j]))),
                    )
                )
    return out


def main(out_dir, quick=False):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for pad in (0.0, 0.5e-3):
        for n in (40,) if quick else (40, 80):
            rows += compare(pad, n)
    (out / "hfss_representation.json").write_text(
        json.dumps({"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "PROVISIONAL_DATASHEET"], "rows": rows}, indent=1)
    )
    lines = ["| pad mm | mesh | model | theta | mode | d|S| | d phase deg |", "|---|---|---|---|---|---|---|"]
    lines += [
        f"| {r['pad_mm']:g} | {r['mesh']} | {r['model']} | {r['theta_deg']:g} | {r['mode']} | {r['d_abs']:+.4f} | {r['d_phase_deg']:+.2f} |"
        for r in rows
    ]
    (out / "hfss_representation.md").write_text("\n".join(lines) + "\n")
    return rows


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    main(a[0] if a else ROOT / "results" / "phase20_r2", quick="--quick" in sys.argv)
