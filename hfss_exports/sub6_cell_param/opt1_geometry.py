"""Analytic geometry and local DC checks for the top-layer bias cell."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math

PERIOD = 11.6
DISC_R = 5.435
SLOT_L = 8.696
SLOT_W = 2.9
ISLAND_W = 1.3
DIODE_GAP = 0.3
DIODE_W = 0.6
CU_T = 0.035
RT_T = 3.175


@dataclass(frozen=True)
class Geometry:
    route: str
    island_l: float
    line_w: float
    channel_w: float

    @property
    def name(self) -> str:
        route = self.route if self.route == "REF" else f"{self.route}_w{self.line_w:.2f}"
        return f"{route}_L{self.island_l:.3f}"


def check(g: Geometry) -> dict:
    if g.route not in ("REF", "Y", "X"):
        raise ValueError(g.route)
    left = -SLOT_W / 2 + DIODE_GAP
    right = left + ISLAND_W
    half_l = g.island_l / 2
    center_x = (left + right) / 2
    line_clearance = (g.channel_w - g.line_w) / 2 if g.route != "REF" else None
    if g.route == "Y":
        rim = math.sqrt(DISC_R**2 - center_x**2)
        outbound = {"center_x": center_x, "from_y": half_l, "to_y": PERIOD / 2}
        channel = {"center_x": center_x, "from_y": SLOT_L / 2, "to_y": PERIOD / 2}
        receiver = {"center_x": center_x, "from_y": -PERIOD / 2, "to_y": -5.0,
                    "overlap_into_disc_mm": rim - 5.0}
        island_to_channel = SLOT_L / 2 - half_l
        receiver_to_island_gap = 5.0 - half_l
    elif g.route == "X":
        rim = DISC_R
        outbound = {"center_y": 0, "from_x": right, "to_x": PERIOD / 2}
        channel = {"center_y": 0, "from_x": SLOT_W / 2, "to_x": PERIOD / 2}
        receiver = {"center_y": 0, "from_x": -PERIOD / 2, "to_x": -5.0,
                    "overlap_into_disc_mm": rim - 5.0}
        island_to_channel = SLOT_W / 2 - right
        receiver_to_island_gap = left - (-5.0)
    else:
        outbound = channel = receiver = None
        island_to_channel = None
        receiver_to_island_gap = None
    local_copper_separated = (
        DIODE_GAP > 0
        and SLOT_W / 2 - right > 0
        and SLOT_L / 2 - half_l > 0
        and (g.route == "REF" or (line_clearance > 0 and receiver_to_island_gap > 0))
    )
    checks = {
        "diode_gap_mm": DIODE_GAP,
        "island_to_slot_y_clearance_mm": SLOT_L / 2 - half_l,
        "island_to_slot_right_clearance_mm": SLOT_W / 2 - right,
        "line_channel_clearance_each_side_mm": line_clearance,
        "outbound_line_separated_from_local_disc": g.route == "REF" or line_clearance > 0,
        "island_to_opposite_receiver_gap_mm": receiver_to_island_gap,
        "island_local_disc_copper_dc_contact": not local_copper_separated,
        "connection_to_neighbor_disc_requires_opposite_edge_receiver": g.route != "REF",
        "receiver_overlaps_local_disc": g.route == "REF" or receiver["overlap_into_disc_mm"] > 0,
    }
    conflicts = []
    if g.route != "REF" and line_clearance < 0.15 - 1e-9:
        conflicts.append(f"{g.line_w:.2f} mm line in {g.channel_w:.2f} mm channel leaves only {line_clearance:.2f} mm per side, not 0.15 mm")
    if g.route != "REF" and island_to_channel <= 0:
        conflicts.append("island reaches channel or disc")
    if not local_copper_separated:
        conflicts.append("local island and disc copper are not separated")
    return {"name": g.name, "geometry": asdict(g), "outbound_line": outbound,
            "channel": channel, "opposite_edge_receiver": receiver, "checks": checks,
            "status": "CONFLICT" if conflicts else "PASS", "conflicts": conflicts}


def all_geometries(wide_channel: float = 0.4) -> list[dict]:
    specs = [Geometry("REF", l, 0, 0) for l in (4.348, 2.5)]
    specs += [Geometry(route, l, w, 0.4 if w == 0.1 else wide_channel)
              for route in ("Y", "X") for l in (4.348, 2.5) for w in (0.1, 0.2)]
    return [check(g) for g in specs]
