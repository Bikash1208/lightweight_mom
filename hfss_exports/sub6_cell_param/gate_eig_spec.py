"""Across-only SD gate-confirm case matrix; 100 degrees occurs once per C state."""
from cell2_geometry import GEOMETRIES as CELL2_GEOMETRIES

GEOMETRIES = [next(g for g in CELL2_GEOMETRIES if g.name == name)
              for name in ("SD_w2p9_ylines", "SD_w1p0_ylines")]
STATES = ("ON", "OFF", "OFF_C013", "OFF_C018")
PHASES_HIGH = tuple(range(100, 181, 10))
PHASES_C_LOW = tuple(range(70, 101, 5))


def variant(g):
    return g.name + "_PEC150"


def cases(g):
    out = []
    for state in STATES:
        phases = sorted(set(PHASES_HIGH) | (set(PHASES_C_LOW) if state.startswith("OFF_C") else set()))
        for phase in phases:
            out.append({"*": f"GATEEIG_{g.name}_PEC150_{state}_across_{phase}",
                        "state": state, "direction": "across", "phase_deg": str(phase)})
    return out
