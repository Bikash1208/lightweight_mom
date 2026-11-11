"""Option 4 baseline and lid-check phase cases."""
from opt4_geometry import GEOMETRIES

LIDS = {"PEC150": (150.0, "PEC"), "PEC300": (300.0, "PEC"),
        "Z377_150": (150.0, "IMPEDANCE_377")}
PHASES = (70, 75, 80, 85, 90, 95, 100)


def variant(g, lid):
    return f"{g.name}_{lid}"


def cases(g, lid):
    prefix = "OPT4_" + variant(g, lid).replace(".", "p")
    return [{"*": f"{prefix}_{state}_{direction}_{phase}", "state": state,
             "direction": direction, "phase_deg": str(phase)}
            for state in ("OFF", "ON") for direction in ("along", "across")
            for phase in PHASES]
