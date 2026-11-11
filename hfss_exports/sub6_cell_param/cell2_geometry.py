"""Cell 2 analytic geometry and periodic-wrap DC gate."""
from dataclasses import dataclass
import math

PERIOD = 11.6
DISC_R = 5.435
SLOT_L = 8.696
SLOT_W = 2.9
LINE_W = 0.1
DIODE_GAP = 0.3
DIODE_W = 0.6


@dataclass(frozen=True)
class Geometry:
    family: str
    split_w: float
    lines: bool

    @property
    def name(self):
        if self.family == "SD":
            return f"SD_w{self.split_w:.1f}_{'ylines' if self.lines else 'noline'}".replace(".", "p")
        return "ID_yline" if self.lines else "ID_noline"


GEOMETRIES = [Geometry("SD", w, lines) for w in (1.0, 2.9) for lines in (False, True)]
GEOMETRIES += [Geometry("ID", 0.0, lines) for lines in (False, True)]


def connected(graph, a, b):
    pending, seen = [a], set()
    while pending:
        x = pending.pop()
        if x == b:
            return True
        if x not in seen:
            seen.add(x)
            pending.extend(graph.get(x, []))
    return False


def check(g: Geometry):
    p = PERIOD / 2
    clearances = {"disc_to_cell_face_mm": p - DISC_R}
    conflicts = []
    if g.family == "SD":
        half = g.split_w / 2
        tab_span = half - DIODE_GAP / 2
        clearances.update({"tab_to_opposite_tab_mm": DIODE_GAP,
                           "cut_to_yline_edge_mm": 1.5 - LINE_W / 2 if g.lines else None,
                           "line_rim_overlap_min_mm": 0.15 if g.lines else None})
        if g.split_w not in (1.0, 2.9) or tab_span <= 0:
            conflicts.append("split/tab dimensions invalid")
        if g.lines and clearances["cut_to_yline_edge_mm"] < 0.15:
            conflicts.append("yline too close to opposite disc half")
        graph = {"LEFT(n)": ["LEFT(n+y)"] if g.lines else [],
                 "LEFT(n+y)": ["LEFT(n)"] if g.lines else [],
                 "RIGHT(n)": ["RIGHT(n+y)"] if g.lines else [],
                 "RIGHT(n+y)": ["RIGHT(n)"] if g.lines else []}
        diode = ["LEFT(n)", "RIGHT(n)"]
        wrap = ["LEFT(n)↔LEFT(n+y)", "RIGHT(n)↔RIGHT(n+y)"] if g.lines else []
        if connected(graph, *diode):
            conflicts.append("LEFT and RIGHT short through lattice wrap")
        if g.lines:
            x = half + 1.5
            outer_x = x + LINE_W / 2
            rim_at_outer = math.sqrt(DISC_R**2 - outer_x**2)
            if rim_at_outer - 0.15 >= p:
                conflicts.append("yline does not reach cell face")
            clearances["line_x_from_centre_mm"] = x
            clearances["line_start_y_mm"] = rim_at_outer - 0.15
    elif g.family == "ID":
        clearances.pop("disc_to_cell_face_mm")
        disc_x = -0.2
        right_rim = disc_x + DISC_R
        left_rim = disc_x - DISC_R
        clearances.update({"disc_left_rim_to_cell_face_mm": left_rim + p,
                           "disc_right_rim_to_cell_face_mm": p - right_rim,
                           "diode_to_x_face_mm": p - 5.60,
                           "landing_to_disc_rim_mm": 5.60 - right_rim,
                           "right_tab_disc_overlap_min_mm": disc_x + math.sqrt(DISC_R**2 - 0.3**2) - 5.10,
                           "left_tab_disc_overlap_min_mm": -5.50 - (disc_x - math.sqrt(DISC_R**2 - 0.3**2)),
                           "diode_span_mm": 0.3,
                           "line_rim_overlap_min_mm": 0.15 if g.lines else None,
                           "line_x_mm": disc_x if g.lines else None})
        if clearances["right_tab_disc_overlap_min_mm"] <= 0 or clearances["left_tab_disc_overlap_min_mm"] <= 0:
            conflicts.append("tab misses shifted disc")
        graph = {"DISC(n)": ["DISC(n+y)"] if g.lines else [],
                 "DISC(n+y)": ["DISC(n)"] if g.lines else [],
                 "LANDING(n)": ["DISC(n+x)"],
                 "DISC(n+x)": ["LANDING(n)"]}
        diode = ["DISC(n)", "LANDING(n)"]
        wrap = ["LANDING(n)↔DISC(n+x)"] + (["DISC(n)↔DISC(n+y)"] if g.lines else [])
        if connected(graph, "DISC(n)", "DISC(n+x)"):
            conflicts.append("adjacent x discs short without diode")
        if connected(graph, *diode):
            conflicts.append("disc shorts to landing without diode")
        if not math.isclose(p - 5.60, 0.20):
            conflicts.append("diode face clearance is not 0.20 mm")
    else:
        raise ValueError(g.family)
    for key, value in clearances.items():
        if value is not None and ("overlap" not in key and "span" not in key and
                                  "start_y" not in key and "from_centre" not in key and
                                  key != "line_x_mm" and
                                  value < 0.15 - 1e-9):
            conflicts.append(f"clearance {key} below 0.15 mm")
    return {"geometry": g.name, "family": g.family, "split_w_mm": g.split_w,
            "lines": g.lines, "status": "PASS" if not conflicts else "CONFLICT",
            "conflicts": conflicts, "clearances_mm": clearances,
            "wrapped_graph_without_diode": graph, "diode_edge": diode,
            "periodic_wrap_connections": wrap,
            "only_diode_connects_target_nets": not connected(graph, *diode) and not conflicts}
