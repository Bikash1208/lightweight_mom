"""Option 4 geometry and conservative analytic DC topology gate."""
from dataclasses import dataclass, asdict
import math

from opt1_geometry import PERIOD, DISC_R, SLOT_L, SLOT_W, ISLAND_W, DIODE_GAP

VIA_DRILL = 0.3
VIA_PAD = 0.6
VIA_DISC_XY = (-3.0, 0.0)
LINE_W = 0.1
CHANNEL_W = 0.4


@dataclass(frozen=True)
class Geometry:
    route: str
    island_l: float

    @property
    def name(self):
        return f"{self.route}_L{self.island_l:.3f}"


GEOMETRIES = [Geometry(route, length) for route in
              ("G-ISL-a", "G-ISL-b", "G-ISL-noline", "G-DISC")
              for length in (4.348, 2.5)]


def check(g: Geometry):
    left = -SLOT_W / 2 + DIODE_GAP
    right = left + ISLAND_W
    cx = (left + right) / 2
    half_l = g.island_l / 2
    is_disc = g.route == "G-DISC"
    via_xy = VIA_DISC_XY if is_disc else (cx, 0.0)
    radius = VIA_PAD / 2
    disc_slot_clearance = abs(via_xy[0]) - SLOT_W / 2 - radius if is_disc else None
    disc_rim_clearance = DISC_R - math.hypot(*via_xy) - radius if is_disc else None
    island_pad_clearance_x = min(via_xy[0] - left, right - via_xy[0]) - radius if not is_disc else None
    island_pad_clearance_y = half_l - radius if not is_disc else None
    line_disc_overlap = DISC_R - 5.0 if g.route in ("G-ISL-a", "G-ISL-b") else None
    channel_clearance = (CHANNEL_W - LINE_W) / 2 if is_disc else None
    conflicts = []
    if g.route not in ("G-ISL-a", "G-ISL-b", "G-ISL-noline", "G-DISC"):
        conflicts.append("unknown route")
    if not (DIODE_GAP > 0 and right < SLOT_W / 2 and half_l < SLOT_L / 2):
        conflicts.append("island touches disc copper outside diode")
    if is_disc:
        if disc_slot_clearance < 0.15 - 1e-9 or disc_rim_clearance < 0.15 - 1e-9:
            conflicts.append("G-DISC pad lacks 0.15 mm slot/rim clearance")
        if channel_clearance < 0.15 - 1e-9:
            conflicts.append("G-DISC line/channel clearance below 0.15 mm")
        if abs(via_xy[0] - cx) < ISLAND_W / 2 + radius and abs(via_xy[1]) < half_l + radius:
            conflicts.append("disc via pad overlaps island")
    else:
        if island_pad_clearance_x < 0 or island_pad_clearance_y < 0:
            conflicts.append("island via pad outside island")
        if g.route in ("G-ISL-a", "G-ISL-b") and line_disc_overlap <= 0:
            conflicts.append("disc line does not overlap disc rim")
    # Copper graph: D=disc, I=island, G=ground. The sheet is the sole D-I edge.
    # Every variant has exactly one direct via edge, either I-G or D-G.
    graph = {"disc": ["ground"] if is_disc else [],
             "island": [] if is_disc else ["ground"], "ground": []}
    anode, cathode = ("island", "disc") if is_disc else ("disc", "island")
    if graph[anode]:
        conflicts.append("anode has a direct copper path to ground")
    return {"name": g.name, "geometry": asdict(g), "status": "PASS" if not conflicts else "CONFLICT",
            "conflicts": conflicts, "via_xy_mm": via_xy, "via_drill_mm": VIA_DRILL,
            "via_pad_mm": VIA_PAD, "disc_via_slot_clearance_mm": disc_slot_clearance,
            "disc_via_rim_clearance_mm": disc_rim_clearance,
            "island_pad_clearance_x_mm": island_pad_clearance_x,
            "island_pad_clearance_y_mm": island_pad_clearance_y,
            "line_to_disc_overlap_mm": line_disc_overlap,
            "line_channel_clearance_each_side_mm": channel_clearance,
            "dc_graph_without_diode": graph, "anode": anode, "cathode": cathode,
            "anode_reaches_ground_only_through_diode": not graph[anode] and not conflicts}
