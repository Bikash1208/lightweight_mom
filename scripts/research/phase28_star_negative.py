"""Phase 28 item 3: the NEGATIVE result for the 'shortenable star' cell (`coded_cell_candidates.ShortenableStarCell`): a diode wall across a slot does not close the slot.

    python scripts/research/phase28_star_negative.py <out.json>

Four through-slots (0, 45, 90, 135 deg) in a 6 mm circular patch on a 6.4 mm cell at 10.0 GHz, mesh 32. State 0 shortens the arms of the perpendicular slot by an ideal-switch
wall (ON = wall conducts). Compared with (a) the same cell with the walls OFF (all arms long, isotropic) and (b) the geometry in which the perpendicular slot is genuinely
shorter (what the wall was meant to achieve). If the walls worked, (state 0) would approach (b).
"""

import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
warnings.simplefilter("ignore")

from src.geometry import UnitCell  # noqa: E402
from src.hologram.coded_cell_candidates import ANGLES_DEG, ShortenableStarCell, _rot  # noqa: E402
from src.hologram.coded_cell_eval import principal_metrics, sheet_tensor  # noqa: E402
from src.materials import Metal  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402
from src.ris.cells import OFF, RT5880, ideal_kit  # noqa: E402

F = 10.0e9


class _Geo:
    kind = "star_lengths"
    polarity = "metal"
    rotation = 0.0

    def __init__(self, lengths_m, width_m=1.6e-3, diameter_m=6.0e-3):
        self.lengths, self.width, self.diameter = lengths_m, width_m, diameter_m
        self.feature_area = 1.0

    def corners(self):
        r = self.diameter / 2
        return np.array([[r, r], [-r, r], [-r, -r], [r, -r]])

    def in_feature(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        m = np.hypot(x, y) <= self.diameter / 2
        for a, length in zip(np.radians(ANGLES_DEG), self.lengths, strict=True):
            u, v = _rot(x, y, a)
            m &= ~((np.abs(u) <= length / 2) & (np.abs(v) <= self.width / 2))
        return m


def _metrics(uc):
    m = principal_metrics(sheet_tensor(uc, F, 32), RT5880, F)
    return {k: m[k] for k in ("b_bar_ms", "delta_ms", "dx_aniso_ohm", "alpha_beam_np_per_m")}


def main(out):
    cell = ShortenableStarCell()
    kit = ideal_kit()
    res = {
        "label": "PHASE28_ITEM3_STAR_CELL_NEGATIVE_RESULT",
        "frequency_ghz": F / 1e9,
        "wall_state_0_ideal_switch": _metrics(cell.unit_cell(cell.states()[0], kit)),
        "all_walls_open_ideal_switch": _metrics(cell.unit_cell((OFF,) * cell.n_diodes, kit)),
    }
    long_m, short_m = cell.slot_length, 2 * cell.wall_radius
    lengths = [long_m, long_m, short_m, long_m]  # axis 0 long, axis 90 genuinely short
    res["perpendicular_slot_really_shortened"] = _metrics(
        UnitCell(RectangularLattice(cell.period, cell.period), _Geo(lengths), RT5880, Metal(), backing="pec")
    )
    Path(out).write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
