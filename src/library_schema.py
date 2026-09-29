"""Phase 7: locked physical configuration, database/CSV column schema, and classification
constants for the production rectangular-PEC-patch unit-cell library.

Every physical value here is either the Phase 7 brief's own locked spec (period, frequency grid,
patch sweep) or read from this project's own authoritative config file, `config.yaml` (substrate
`eps_r`/`tan_delta`/`name`; the brief's substrate thickness, 1.575 mm, already matches it exactly).
No electrical value is invented; see `LibraryConfig.default()`'s docstring for the provenance of
each field. `LibraryConfig` also accepts overrides (smaller L/W lists, a smaller frequency grid,
a different mesh) so `tests/` can exercise the same code paths at a tractable scale without
claiming to be "the" library -- only `LibraryConfig.default()` is the locked Phase 7 spec.

This module defines *schema*, not behaviour: it has no solver calls. `library_builder.py` runs
the sweep; `phase_extraction.py`, `library_classifier.py`, `geometry_export.py`,
`library_plots.py`, `library_validation.py` and `unit_cell_library.py` consume it.
"""

from dataclasses import dataclass, field
from pathlib import Path

from . import units
from .config import config_from_dict

# --------------------------------------------------------------------- authoritative source

# config.yaml (this project's own root example config, unmodified) is the authoritative source
# for these three substrate values; the Phase 7 brief's substrate thickness (1.575 mm) already
# matches it exactly, so nothing here is invented or overridden.
_AUTHORITATIVE_CONFIG_PATH = Path(__file__).resolve().parent.parent / "configs" / "config.yaml"
_AUTHORITATIVE_SUBSTRATE_NAME = "Rogers RT/Duroid 5880"
_AUTHORITATIVE_EPS_R = 2.2
_AUTHORITATIVE_TAN_DELTA = 0.0009
_AUTHORITATIVE_THICKNESS_M = 1.575e-3  # matches the Phase 7 brief's locked h = 1.575 mm exactly

# ----------------------------------------------------------------------------- locked sweep

_PERIOD_M = 4.8e-3  # period_x = period_y, brief Sec. 1
_F_MIN_HZ = 10.0e9
_F_MAX_HZ = 15.0e9
_F_STEP_HZ = 25.0e6  # (15 - 10) GHz / 200 steps = 25 MHz -> 201 points inclusive
_F_DESIGN_HZ = 12.5e9
_N_FREQ = 201

_L_MIN_M, _L_MAX_M, _L_STEP_M = 1.0e-3, 4.4e-3, 0.2e-3  # 18 values
_N_GEOM_PER_AXIS = 18


def _linspace_mm(start_mm, stop_mm, step_mm, n):
    """18 values 1.0, 1.2, ..., 4.4 mm, built by exact-count `round` (not float accumulation:
    1.0 + 17*0.2 must land on 4.4 mm bit-for-bit-reproducibly, not `4.399999999999999`)."""
    return tuple(round(start_mm + i * step_mm, 10) for i in range(n))


L_VALUES_MM = _linspace_mm(1.0, 4.4, 0.2, _N_GEOM_PER_AXIS)
W_VALUES_MM = L_VALUES_MM
assert L_VALUES_MM[0] == 1.0 and L_VALUES_MM[-1] == 4.4 and len(L_VALUES_MM) == 18

# Mesh: cells_x = cells_y = 32 (dx = dy = 0.15 mm), floquet_order = 32. See LIBRARY_REPORT.md
# Sec. "Mesh and discretization quality" and tests/test_library_schema.py for how this was
# chosen: a 24x24 mesh (dx = 0.2 mm, exactly the L/W step) looked ideal for zero staircase error
# but empirically produces a genuine floating-point misclassification of an entire extra mesh row
# for many (L, W) pairs (metal_area_error up to 44%, not an ordinary staircase pattern -- this is
# the same class of pitfall documented in validation_cases.py's module docstring, generalized:
# an unmodified solver-internal boundary comparison is sensitive to a few ULP of rounding when a
# feature edge sits exactly on a mesh line). 32x32 avoids exact coincidence (dx = 0.15 mm shares
# no simple ratio with the 0.2 mm sweep step) and gives an ordinary, smoothly-varying staircase
# error (verified non-pathological); the a priori `metal_area_error`/`area_ok` diagnostic
# (`discretization.check_cell`, Phase 5, unmodified, run for every cell by `runner.run_sweep`) is
# stored per cell rather than hidden, using this project's own existing
# `discretization.DEFAULT_MAX_AREA_ERROR` (0.05) policy threshold -- not a new one invented here.
# floquet_order = 32 satisfies the rank rule (2*32+1 = 65 >= 32) with a large margin; it is below
# the "2-4x cells" practical accuracy rule of thumb (MOM_FORMULATION.md) purely for runtime cost
# (see LIBRARY_REPORT.md for the measured build-time trade-off across mesh/Floquet choices, and
# for a `discretization.refinement_study` sample run that quantifies -- not assumes -- the
# resulting truncation error on a handful of representative cells).
DEFAULT_MESH = {"cells_x": 32, "cells_y": 32, "floquet_order": 32}

# Reference planes: not locked by the Phase 7 brief (they do not affect Z or V -- excitation.py,
# S_PARAMETERS.md); reused unchanged from the Phase 6 validation cases for consistency across the
# whole project's exported artifacts (5 mm of air modelled on each side of the unit cell).
DEFAULT_REFERENCE_PLANES = {"top_air_thickness": 5.0e-3, "bottom_air_thickness": 5.0e-3}


@dataclass(frozen=True)
class LibraryConfig:
    """Everything that defines one library build. `.default()` is the Phase 7 locked spec;
    other constructions (smaller L/W lists, fewer frequencies) exist only for tests, and are
    never passed to `examples/build_unit_cell_library.py`."""

    period_m: float
    substrate: dict  # name, eps_r, tan_delta, thickness_m
    f_min_hz: float
    f_max_hz: float
    f_step_hz: float
    f_design_hz: float
    l_values_mm: tuple
    w_values_mm: tuple
    mesh: dict = field(default_factory=lambda: dict(DEFAULT_MESH))
    reference_planes: dict = field(default_factory=lambda: dict(DEFAULT_REFERENCE_PLANES))

    @classmethod
    def default(cls):
        """The locked Phase 7 configuration (brief Sec. 1): period 4.8 mm, substrate from
        `config.yaml` (this project's authoritative source; matches the brief's own locked
        values exactly), 10-15 GHz / 25 MHz / 201 points / 12.5 GHz design frequency, 18x18
        L/W patch sweep 1.0-4.4 mm in 0.2 mm steps."""
        return cls(
            period_m=_PERIOD_M,
            substrate={
                "name": _AUTHORITATIVE_SUBSTRATE_NAME,
                "eps_r": _AUTHORITATIVE_EPS_R,
                "tan_delta": _AUTHORITATIVE_TAN_DELTA,
                "thickness_m": _AUTHORITATIVE_THICKNESS_M,
            },
            f_min_hz=_F_MIN_HZ,
            f_max_hz=_F_MAX_HZ,
            f_step_hz=_F_STEP_HZ,
            f_design_hz=_F_DESIGN_HZ,
            l_values_mm=L_VALUES_MM,
            w_values_mm=W_VALUES_MM,
        )

    def __post_init__(self):
        if self.period_m <= 0.0:
            raise ValueError(f"period_m must be > 0, got {self.period_m!r}")
        if max(self.l_values_mm) * 1e-3 >= self.period_m or max(self.w_values_mm) * 1e-3 >= self.period_m:
            raise ValueError(
                f"largest L/W ({max(self.l_values_mm)}/{max(self.w_values_mm)} mm) must stay "
                f"strictly inside the {units.to_mm(self.period_m):g} mm period"
            )
        n_expected = round((self.f_max_hz - self.f_min_hz) / self.f_step_hz) + 1
        if not (self.f_min_hz <= self.f_design_hz <= self.f_max_hz):
            raise ValueError(f"f_design_hz={self.f_design_hz!r} must be within [f_min_hz, f_max_hz]")
        object.__setattr__(self, "_n_freq_expected", n_expected)

    @property
    def n_geometries(self):
        return len(self.l_values_mm) * len(self.w_values_mm)

    def base_config_dict(self):
        """The shared (non-swept) `config_from_dict` mapping: lattice, substrate, mesh, planes,
        and a placeholder geometry (overwritten per cell by the sweep's length_m/width_m axes)."""
        return {
            "period_x": self.period_m,
            "period_y": self.period_m,
            "substrate": {
                "name": self.substrate["name"],
                "eps_r": self.substrate["eps_r"],
                "tan_delta": self.substrate["tan_delta"],
                "thickness": self.substrate["thickness_m"],
            },
            "geometry": {
                "type": "rectangle",
                "length": self.l_values_mm[0] * 1e-3,
                "width": self.w_values_mm[0] * 1e-3,
            },
            "excitation": {"frequency": self.f_design_hz},
            "mesh": dict(self.mesh),
            "reference_planes": dict(self.reference_planes),
        }

    def build_config(self):
        """A validated `SimulationConfig` for this library (Phase 1 `config_from_dict`,
        unmodified)."""
        return config_from_dict(self.base_config_dict())


# ------------------------------------------------------------------------------ CSV / record schema

# Req. 20: the flat per-(cell, frequency) CSV. `library_builder.py` populates it from the
# `database.Database` this config's sweep produces -- same numbers, same exp(+jwt)/S-parameter
# conventions (database.CONVENTIONS), reformatted to these exact column names/units (mm, GHz,
# dB) rather than a second, independently-derived S-parameter path (req. 23).
_S_FAMILIES = ("S11", "S21", "S12", "S22")
_S_FIELDS = ("real", "imag", "mag", "mag_db", "phase_deg")

LIBRARY_CSV_COLUMNS = ("cell_id", "L_mm", "W_mm", "f_GHz") + tuple(
    f"{fam}_{field_}" for fam in _S_FAMILIES for field_ in _S_FIELDS
)

# Req. 5's minimum per-record field list (JSON, one entry per (cell, frequency)); the JSON export
# additionally nests re/im/mag/mag_db/phase_deg for readability instead of flattening them.
LIBRARY_RECORD_FIELDS = (
    "cell_id",
    "L_mm",
    "W_mm",
    "period_x_mm",
    "period_y_mm",
    "substrate_thickness_mm",
    "epsilon_r",
    "tan_delta",
    "frequency_GHz",
    "S11",
    "S21",
    "S12",
    "S22",
    "polarization_convention",
    "time_convention",
    "geometry_status",
    "solver_status",
)

# ----------------------------------------------------------------------------- classification

# Req. 13: every classification field, its exact criterion, and why -- all consumed by
# library_classifier.py, none hidden or "just because". Independent flags (not one exclusive
# label): a cell can be RESONANT and PHASE_USABLE at once, so a single enum would either lose
# information or force an arbitrary priority order that Sec. 13 never asked for.
GEOMETRY_STATUS = ("VALID", "INVALID")
SOLVER_STATUS = ("OK", "FAILED")
RESONANCE_STATUS = ("RESONANT", "NON_RESONANT", "UNCLASSIFIED")
PHASE_STATUS = ("PHASE_USABLE", "PHASE_UNUSABLE", "UNCLASSIFIED")
GRATING_STATUS = ("GRATING_RISK", "OK")

# Resonance: a local maximum of |S11_xx(f)| at a strictly interior frequency sample, with
# prominence (peak value minus the lower of its two flanking minima -- phase_extraction.py's
# `_prominence`) at or above this. |S11_xx| in [0, 1] (a passive, normalized-incidence
# reflection magnitude), so 0.10 means "the peak stands at least 0.10 above the surrounding
# response" -- an explicit, documented, and tested policy value, not a theorem; see
# tests/test_library_classifier.py for the boundary check on synthetic curves and
# LIBRARY_REPORT.md for the measured distribution over the actual library.
RESONANCE_MIN_PROMINENCE = 0.10

# Phase usability: below this |S11_xx| magnitude (dB) at the design frequency, the phase of an
# almost-zero complex number is dominated by floating-point/model noise rather than physical
# information (the same reasoning already applied to the -300 dB CSV floor in
# validation_export.py/hfss_comparison.py, generalized into an actual usability threshold here).
PHASE_USABLE_MIN_MAG_DB = -40.0

# Phase-collision bin width (req. 11): geometries within this many degrees of each other at the
# design frequency are reported as candidates for the same phase region. A policy value (not a
# theorem); chosen as 1/36 of the full circle (10 degrees) for a readable number of bins.
PHASE_COLLISION_BIN_DEG = 10.0

# Default tolerance for library.find_phase_candidates / invert_phase (req. 24-25) when the
# caller does not specify one explicitly.
DEFAULT_PHASE_TOLERANCE_DEG = 5.0
