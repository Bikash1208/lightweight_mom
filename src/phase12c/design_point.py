"""The Phase 12C working design point, loaded from ``config_phase12c_grounded_8mm.yaml`` (brief's
own instruction: "Record it as a config; don't hard-code it") rather than repeated as module-level
constants the way phase12a/phase12b's own study modules did. Every other Phase 12C module builds
its unit cells from :func:`unit_cell_at`, which reuses this config's period/substrate/metal/backing
and only ever varies the patch length/width -- the one thing every task in this phase sweeps.
"""

from pathlib import Path

from ..config import load_config
from ..geometry import Geometry, UnitCell

CONFIG_PATH = Path(__file__).resolve().parent.parent.parent / "configs" / "config_phase12c_grounded_8mm.yaml"

DESIGN_CONFIG = load_config(CONFIG_PATH)

LATTICE = DESIGN_CONFIG.unit_cell.lattice
SUBSTRATE = DESIGN_CONFIG.unit_cell.substrate
METAL = DESIGN_CONFIG.unit_cell.metal
BACKING = DESIGN_CONFIG.unit_cell.backing
REFERENCE_PLANES = DESIGN_CONFIG.reference_planes
DESIGN_FREQUENCY_HZ = DESIGN_CONFIG.excitation.frequency
PERIOD_M = LATTICE.period_x
BASE_MESH_CELLS = DESIGN_CONFIG.mesh.cells_x
BASE_FLOQUET_ORDER = DESIGN_CONFIG.mesh.floquet_order
REPRESENTATIVE_LENGTH_WIDTH_M = (DESIGN_CONFIG.unit_cell.geometry.length, DESIGN_CONFIG.unit_cell.geometry.width)

assert LATTICE.period_x == LATTICE.period_y, "Phase 12C's design point is a square lattice"
assert BACKING == "pec", "config_phase12c_grounded_8mm.yaml must set substrate.backing: pec"


def unit_cell_at(length_m, width_m=None):
    """A :class:`~src.geometry.UnitCell` at this design point's period/substrate/metal/backing,
    with a rectangular patch ``length_m`` x ``width_m`` (default: square, ``length_m`` x
    ``length_m``) -- every Phase 12C sweep's only varying quantity."""
    width_m = length_m if width_m is None else width_m
    geometry = Geometry(kind="rectangle", length=length_m, width=width_m)
    return UnitCell(LATTICE, geometry, SUBSTRATE, METAL, backing=BACKING)
