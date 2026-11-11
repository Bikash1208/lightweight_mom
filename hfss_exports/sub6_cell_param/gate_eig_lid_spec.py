"""ON/OFF 70..180/10 at 40 and 300 mm; six missing rows at 150 mm."""
from gate_eig_spec import GEOMETRIES

GEOMETRY = next(g for g in GEOMETRIES if g.name == "SD_w2p9_ylines")
LIDS = (40, 150, 300)
STATES = ("ON", "OFF")
PHASES = tuple(range(70, 181, 10))


def variant(lid):
    return f"{GEOMETRY.name}_PEC{lid}"


def cases(lid):
    phases = (70, 80, 90) if lid == 150 else PHASES
    return [{"*": f"GATEEIG_{GEOMETRY.name}_PEC{lid}_{state}_across_{phase}",
             "state": state, "direction": "across", "phase_deg": str(phase)}
            for state in STATES for phase in phases]
