"""Phase 28 (loop cell R4, RT/duroid 5880 3.175 mm, 9.0 mm period): MoM solves for (a) the rotation check (cell rotated by 20 deg inside the lattice) and (b) the exact HFSS-package representation
(diode as SEPARATE single-type lumped elements joined by perfect-conductor islands when L_s > 0, `hfss_chain`). MODEL_ONLY, NOT_HFSS_VALIDATED.

    python scripts/research/phase28_loop_rotation.py <out.jsonl> [--workers N]

Cell R4: closed square loop, outer side 6.6 mm, strip width 1.2 mm, diode region 0.6 mm (along the local x axis) x 1.2 mm in the middle of the top side. Rotation by psi turns the loop AND the diode axis about the cell centre
(psi = 20 deg keeps the loop corners inside the 9 mm cell: 3.3 (cos psi + sin psi) = 4.23 mm < 4.5 mm; 40 deg would overlap the neighbours: 4.65 mm).
Diode cases: L0 = MADP, L_s 0 (ON 5.2 ohm, OFF 0.025 pF, single element); L02 = MADP with L_s = 0.2 nH ASSUMED (ON: R then L, OFF: C then L; elements 0.25 mm, islands 0.10 mm).
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
import warnings  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from dataclasses import dataclass  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

PITCH, SO, WR, GAP, ISLAND = 9.0e-3, 6.6e-3, 1.2e-3, 0.6e-3, 0.1e-3
SUB = (2.2, 0.0009, 3.175e-3)


@dataclass(frozen=True)
class RotLoop:
    outer: float
    width: float
    angle_deg: float
    kind = "loop_patch_rot"
    polarity = "metal"
    rotation = 0.0

    @property
    def feature_area(self):
        return self.outer**2

    def corners(self):
        a = math.radians(self.angle_deg)
        h = self.outer / 2 * (abs(math.cos(a)) + abs(math.sin(a)))
        return np.array([[h, h], [-h, h], [-h, -h], [h, -h]])

    def in_feature(self, x, y):
        a = math.radians(self.angle_deg)
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        u, v = x * math.cos(a) + y * math.sin(a), -x * math.sin(a) + y * math.cos(a)
        h = self.outer / 2 + 1e-12
        inner = self.outer / 2 - self.width - 1e-12
        return (np.abs(u) <= h) & (np.abs(v) <= h) & ~((np.abs(u) < inner) & (np.abs(v) < inner))


def diode_model(case, on, ls_nh=None):
    from src.ris.diodes import DiodeModel

    ls = (0.2e-9 if case == "L02" else 0.0) if ls_nh is None else ls_nh * 1e-9
    return DiodeModel("ON", R_s=5.2, L_s=ls) if on else DiodeModel("OFF", L_s=ls, C_T=0.025e-12)


def segments(case, on):
    """Series elements of the diode along its axis, in HFSS order: [(kind, value)], kind in R, C, L."""
    segs = [("R", 5.2)] if on else [("C", 0.025e-12)]
    if case == "L02":
        segs.append(("L", 0.2e-9))
    return segs


def build(case, on, rot_deg):
    from src.geometry import UnitCell
    from src.loads import LumpedLoad
    from src.materials import Metal, Substrate
    from src.periodic_bc import RectangularLattice
    from src.ris.cells import RT5880  # noqa: F401

    a = math.radians(rot_deg)
    x0, y0 = 0.0, SO / 2 - WR / 2
    cx, cy = x0 * math.cos(a) - y0 * math.sin(a), x0 * math.sin(a) + y0 * math.cos(a)
    segs = segments(case, on)
    n = len(segs)
    elem = (GAP - (n - 1) * ISLAND) / n
    loads, pos = [], -GAP / 2
    direction = "x" if rot_deg == 0 else float(rot_deg)
    for kind, val in segs:
        c = pos + elem / 2
        ox, oy = c * math.cos(a), c * math.sin(a)
        kw = {"R": dict(R=val), "L": dict(L=val), "C": dict(C=val)}[kind]
        loads.append(LumpedLoad(cx + ox, cy + oy, elem, WR, direction, transverse_conducts=False, **kw))
        pos += elem + ISLAND
    sub = Substrate("rt5880_3p175", *SUB)
    return UnitCell(RectangularLattice(PITCH, PITCH), RotLoop(SO, WR, rot_deg), sub, Metal(), backing="pec", loads=tuple(loads))


def solve(job):
    case, on, rot, mesh, f_hz = job
    warnings.simplefilter("ignore")
    from src.hologram.coded_cell_eval import sheet_tensor

    y = sheet_tensor(build(case, on, rot), f_hz, mesh)
    return {"case": case, "on": on, "rot_deg": rot, "mesh": mesh, "f_hz": f_hz, "y_s": [[[float(v.real), float(v.imag)] for v in row] for row in y]}


def jobs():
    out = []
    for on in (True, False):
        out += [("L0", on, 20, 60, 10e9), ("L0", on, 20, 90, 10e9), ("L02", on, 20, 60, 10e9), ("L02", on, 0, 90, 10e9)]
        out += [("L02", on, 0, 60, f * 1e9) for f in (9.5, 9.75, 10.0, 10.25, 10.5)]
        out += [("L0", on, 0, 60, f * 1e9) for f in (9.5, 9.75, 10.0, 10.25, 10.5)]
    return out


def main(out, workers):
    out = Path(out)
    done = set()
    if out.exists():
        for line in out.read_text(encoding="utf-8").splitlines():
            r = json.loads(line)
            done.add((r["case"], r["on"], r["rot_deg"], r["mesh"], r["f_hz"]))
    todo = [j for j in jobs() if j not in done]
    print(len(todo), "solves", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex, out.open("a", encoding="utf-8") as fh:
        for r in ex.map(solve, todo):
            fh.write(json.dumps(r) + "\n")
            fh.flush()


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], int(a[a.index("--workers") + 1]) if "--workers" in a else 30)
