"""Phase 28 item 3: a consistency check of the lumped-diode treatment in the candidate cells. With an IDEAL switch (ON = 1 mohm, OFF = 1e-18 F) every diode neck is either a plain
metal bridge (ON) or an empty gap (OFF). The same states are built a second way WITHOUT any lumped load: the ON necks are plain metal in the pattern, the OFF necks are removed from
the pattern. The two sheet tensors must agree (they would not if the loads were mis-placed, mis-oriented or mis-weighted).

    python scripts/research/phase28_neck_load_check.py <out.json>
"""

import json
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
warnings.simplefilter("ignore")

from src.geometry import UnitCell  # noqa: E402
from src.hologram.coded_cell_candidates import RationalSectorCell, RationalSectorGeometry  # noqa: E402
from src.hologram.coded_cell_eval import principal_metrics, sheet_tensor  # noqa: E402
from src.materials import Metal  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402
from src.ris.cells import OFF, ON, RT5880, ideal_kit  # noqa: E402

F = 10.0e9


@dataclass(frozen=True)
class _GeoNoNecks(RationalSectorGeometry):
    """The same pattern with only the listed necks present (`keep`: set of (kind, k))."""

    keep: tuple = ()

    def in_feature(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        full = RationalSectorGeometry.in_feature(self, x, y)
        all_necks = np.zeros(x.shape, dtype=bool)
        kept = np.zeros(x.shape, dtype=bool)
        for kind in ("arm_neck", "corner_neck"):
            for k in range(4):
                x0, x1, y0, y1 = self.neck_rect("arm_neck" if kind == "arm_neck" else "diag", k)
                inside = (x >= x0 - 1e-12) & (x <= x1 + 1e-12) & (y >= y0 - 1e-12) & (y <= y1 + 1e-12)
                all_necks |= inside
                if (kind, k) in self.keep:
                    kept |= inside
        # remove the neck regions that are not kept, but not where another piece (hub / arm / diagonal) covers them
        base = _base_pieces(self, x, y)
        return (full & ~all_necks) | kept | (full & base)


def _base_pieces(geo, x, y):
    u, v = np.abs(x), np.abs(y)
    r_in, t = geo.r_in, geo.gap_t
    m = np.maximum(u, v)
    hub = m <= geo.hub_half
    arm_x = (u >= r_in) & (u <= geo.arm_outer) & (v <= u / 2 - t)
    arm_y = (v >= r_in) & (v <= geo.arm_outer) & (u <= v / 2 - t)
    diag = (m >= r_in) & (m <= geo.diag_outer) & (v >= u / 2 + t) & (u >= v / 2 + t)
    return hub | arm_x | arm_y | diag


def compare_state(idx, mesh=32):
    """Relative difference of the sheet tensor of state `idx` built with ideal-switch loads and built as a plain pattern (the second way, no loads), plus the two leakages."""
    cell = RationalSectorCell(
        period=6.4e-3,
        hub_half=1.2e-3,
        neck_gap=0.4e-3,
        neck_width=0.4e-3,
        gap_t=0.2e-3,
        arm_outer=3.0e-3,
        diag_outer=2.8e-3,
        clip=3.0e-3,
        mode="notch",
    )
    geo = cell.geometry()
    sites = cell._sites()
    st = cell.states()[idx]
    with_loads = sheet_tensor(cell.unit_cell(st, ideal_kit()), F, mesh)
    keep = tuple((s[0], s[1]) for s, v in zip(sites, st, strict=True) if v == ON)
    geo2 = _GeoNoNecks(
        geo.hub_half,
        geo.neck_gap,
        geo.neck_width,
        geo.gap_t,
        geo.arm_outer,
        geo.diag_outer,
        geo.clip,
        geo.neck_width_diag,
        keep,
    )
    plain = sheet_tensor(
        UnitCell(RectangularLattice(cell.period, cell.period), geo2, RT5880, Metal(), backing="pec"), F, mesh
    )
    m1, m2 = principal_metrics(with_loads, RT5880, F), principal_metrics(plain, RT5880, F)
    return {
        "state": idx,
        "mesh": mesh,
        "n_on": sum(v == ON for v in st),
        "n_off": sum(v == OFF for v in st),
        "relative_tensor_difference": float(np.linalg.norm(with_loads - plain) / np.linalg.norm(plain)),
        "delta_ms_with_loads": m1["delta_ms"],
        "delta_ms_plain_pattern": m2["delta_ms"],
        "alpha_with_loads": m1["alpha_beam_np_per_m"],
        "alpha_plain_pattern": m2["alpha_beam_np_per_m"],
    }


def main(out):
    res = {
        "label": "PHASE28_ITEM3_NECK_LOAD_CONSISTENCY",
        "frequency_ghz": F / 1e9,
        "states": [compare_state(i, 32) for i in (0, 1)],
    }
    Path(out).write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
