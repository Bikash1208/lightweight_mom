"""Option 1b geometry, top-boundary, and phase matrix."""
from opt1_geometry import Geometry, all_geometries

GEOMETRY_NAMES = ("REF_L4.348", "Y_w0.10_L4.348", "Y_w0.10_L2.500")
GEOMETRIES = [Geometry(**item["geometry"]) for item in all_geometries(0.5)
              if item["name"] in GEOMETRY_NAMES]
LIDS = (("PEC150", 150.0, "PEC"), ("PEC300", 300.0, "PEC"),
        ("Z377_150", 150.0, "IMPEDANCE_377"))
PHASES = (70, 75, 80, 85, 90, 95, 100)


def variant(g: Geometry, lid: str) -> str:
    return f"{g.name}_{lid}"


def cases(g: Geometry, lid: str) -> list[dict]:
    prefix = "OPT1B_" + variant(g, lid).replace(".", "p")
    return [{"*": f"{prefix}_{state}_{direction}_{phase}", "state": state,
             "direction": direction, "phase_deg": str(phase)}
            for state in ("OFF", "ON") for direction in ("along", "across")
            for phase in PHASES]
