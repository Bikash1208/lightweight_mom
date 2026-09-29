"""Phase 8A: locked configuration, gate states and CLI exit codes for the scientific feasibility
gate (`feasibility_gate.py`). This module defines *schema*, not behaviour, mirroring
`library_schema.py`'s split for Phase 7.
"""

from dataclasses import dataclass
from pathlib import Path

from .library_schema import DEFAULT_PHASE_TOLERANCE_DEG

# The Phase 7 production library this gate reads by default (brief Sec. 1). A different
# ``library_dir`` may be passed for tests, exactly like `LibraryConfig` accepts overrides.
DEFAULT_LIBRARY_DIR = Path(__file__).resolve().parent.parent / "results" / "unit_cell_library"

# Req. 35: explicit states, never a generic boolean. Every value below is actually returned by
# `feasibility_gate.run_feasibility_gate`; none is aspirational.
STATES = (
    "LIBRARY_INVALID",
    "PHASE_SOURCE_UNVALIDATED",
    "REFERENCE_INVALID",
    "FREQUENCY_UNSUPPORTED",
    "POLARIZATION_UNSUPPORTED",
    "GRID_MAPPING_UNVALIDATED",
    "PHASE_COVERAGE_INSUFFICIENT",
    "PHASE_MATCHING_INSUFFICIENT",
    "GRATING_RISK",
    "DESIGN_FEASIBLE",
    "DESIGN_FEASIBLE_WITH_WARNINGS",
)

# Req. 47's exit codes, extended (documented, not silently reused) for the three brief states
# (REFERENCE_INVALID, PHASE_MATCHING_INSUFFICIENT, GRATING_RISK) that Sec. 47's own list omits.
EXIT_CODES = {
    "DESIGN_FEASIBLE": 0,
    "DESIGN_FEASIBLE_WITH_WARNINGS": 1,
    "PHASE_COVERAGE_INSUFFICIENT": 2,
    "LIBRARY_INVALID": 3,
    "PHASE_SOURCE_UNVALIDATED": 4,
    "INPUT_INVALID": 5,
    "FREQUENCY_UNSUPPORTED": 6,
    "POLARIZATION_UNSUPPORTED": 7,
    "GRID_MAPPING_UNVALIDATED": 8,
    "REFERENCE_INVALID": 9,
    "PHASE_MATCHING_INSUFFICIENT": 10,
    "GRATING_RISK": 11,
}

# Reused, not reinvented: the same default query tolerance Phase 7's `UnitCellLibrary` already
# uses for phase candidate matching (`library_schema.DEFAULT_PHASE_TOLERANCE_DEG`).
DEFAULT_PHASE_MATCH_TOLERANCE_DEG = DEFAULT_PHASE_TOLERANCE_DEG

# This solver only ever excites/records the *linear* x/y basis at the reference planes
# (`database.CONVENTIONS['linear_basis']`); circular is a derived basis (`polarization.py`), not
# an independently excited/stored one, and the Phase 7 phase-vs-geometry pipeline this gate reads
# is pinned to the co-polarized "_xx" component only (`phase_extraction.py`). "y" is listed as a
# recognized request so it can be explicitly rejected (Sec. 27: never silently substituted),
# rather than silently accepted and quietly answered with x-polarized data.
SUPPORTED_POLARIZATIONS = ("x", "y")
PROVIDED_POLARIZATION = "x"


@dataclass(frozen=True)
class FeasibilityConfig:
    """Everything one Phase 8A gate run needs. Mode A (library-only) needs nothing but
    ``library_dir``; Mode B additionally needs a reference field / beam set, passed separately to
    `feasibility_gate.run_feasibility_gate` (not stored here, since they are per-design inputs,
    not configuration)."""

    library_dir: Path = DEFAULT_LIBRARY_DIR
    frequency_hz: float | None = None          # None -> the library's own design frequency
    polarization: str = "x"
    phase_match_tolerance_deg: float = DEFAULT_PHASE_MATCH_TOLERANCE_DEG
    expected_n_geometries: int = 324
    expected_n_frequency_points: int = 201

    def __post_init__(self):
        object.__setattr__(self, "library_dir", Path(self.library_dir))
        if self.polarization not in SUPPORTED_POLARIZATIONS:
            raise ValueError(f"polarization must be one of {SUPPORTED_POLARIZATIONS}, got {self.polarization!r}")
        if self.phase_match_tolerance_deg <= 0.0:
            raise ValueError(f"phase_match_tolerance_deg must be > 0, got {self.phase_match_tolerance_deg!r}")
