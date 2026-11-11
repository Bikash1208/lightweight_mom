"""Cell 2 six-geometry, seven-phase PEC150 sweep."""
from cell2_geometry import GEOMETRIES

PHASES = (70, 75, 80, 85, 90, 95, 100)


def variant(g):
    return g.name + "_PEC150"


def cases(g):
    prefix = "CELL2_" + variant(g)
    return [{"*": f"{prefix}_{state}_{direction}_{phase}", "state": state,
             "direction": direction, "phase_deg": str(phase)}
            for state in ("OFF", "ON") for direction in ("along", "across")
            for phase in PHASES]
