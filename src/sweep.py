"""Sweep definition and planning: frequency grid, geometry-parameter axes, cell expansion.

Nothing here solves anything. ``plan_sweep`` expands the definition into the full list of unit
cells and validates *all* of them up front, so a bad parameter is reported before the first
(potentially long) solve rather than halfway through.

Frequency grid  f_start -> f_stop -> step
    f_i = f_start + i * step,  i = 0 .. N-1,  N = (f_stop - f_start)/step + 1,
    both end points included, the last point set exactly to f_stop. (f_stop - f_start) must be
    an integer multiple of step (relative tolerance 1e-9 of the step); anything else is an
    error, so a requested end point is never silently dropped or overshot. A single point is
    f_start == f_stop. Every frequency is an independent solve: no state is shared between them.

Geometry parameters (the only swept quantities)
    length_m, width_m, center_x_m, center_y_m [m] and rotation_deg [deg, counter-clockwise]:
    the fields of ``Geometry``. Lattice periods, substrate, metal, mesh and reference planes
    are fixed for a sweep (a different lattice or stack is a different database), and so is the
    geometry kind. Axes combine as a Cartesian product; the first axis varies slowest, so the
    cell order is deterministic. Non-swept parameters keep the base configuration's value.

Cell identity
    ``cell_id`` = short SHA-256 of the canonical JSON of (sweep context, cell parameters); the
    context holds lattice, substrate, metal, geometry kind, mesh, reference planes and the code
    version. The frequency is not part of it: the checkpoint key of one result is
    (cell_id, frequency, a-priori policy).
"""

import itertools
import math
import warnings
from dataclasses import dataclass, field

import numpy as np

from . import __version__, units
from ._canonical import cell_id as make_cell_id
from ._validate import finite_real
from .discretization import (
    DEFAULT_MAX_AREA_ERROR,
    DEFAULT_MIN_CELLS_PER_WAVELENGTH,
    DEFAULT_REFINEMENT_FACTOR,
    DEFAULT_REFINEMENT_TOLERANCE,
    check_cell,
    check_frequency,
    duplicate_geometry_groups,
    metal_mask_key,
)
from .frequency import validate_frequency
from .geometry import AnnularSlotPatch, Geometry, UnitCell

PARAMETER_NAMES = (
    "length_m", "width_m", "center_x_m", "center_y_m", "rotation_deg",
    # PHASE15 Part E: `geometry.AnnularSlotPatch`'s own four shape parameters -- only meaningful
    # when `sweep_context(config)["geometry_kind"] == AnnularSlotPatch.kind`
    # (`geometry_from_parameters`/`base_parameters` dispatch on that; `ParameterAxis` itself has no
    # notion of geometry kind, so both families' names must coexist in this one flat tuple).
    "a_m", "b_m", "R_m", "r_m",
)
_ANNULAR_SLOT_PARAMETER_NAMES = ("a_m", "b_m", "R_m", "r_m", "center_x_m", "center_y_m", "rotation_deg")
ON_INVALID = ("error", "skip")
REFINEMENT_FREQUENCIES = ("max", "all")

MAX_POINTS = 100_000  # one axis / frequency grid
MAX_CELLS = 1_000_000  # Cartesian product of the axes
_RANGE_RTOL = 1e-9  # of the step, for "span is an integer number of steps"


class SweepError(ValueError):
    """Invalid sweep definition or a cell/frequency that fails validation."""


# ------------------------------------------------------------------------- ranges and grids


def inclusive_range(start, stop, step, what="range"):
    """start, start+step, ..., stop  (both included); see the module docstring."""
    try:
        start = finite_real(f"{what}.start", start)
        stop = finite_real(f"{what}.stop", stop)
        step = finite_real(f"{what}.step", step)
    except ValueError as e:
        raise SweepError(str(e)) from None
    if step <= 0.0:
        raise SweepError(f"{what}: step must be > 0, got {step}")
    if stop < start:
        raise SweepError(f"{what}: stop ({stop}) must be >= start ({start})")
    span = stop - start
    n = int(round(span / step))
    if abs(span - n * step) > _RANGE_RTOL * step:
        raise SweepError(
            f"{what}: (stop - start) = {span} is not an integer multiple of step = {step}; "
            f"the end point would be skipped or overshot"
        )
    if n + 1 > MAX_POINTS:
        raise SweepError(f"{what}: {n + 1} points exceeds the limit of {MAX_POINTS}")
    values = start + np.arange(n + 1) * step
    values[-1] = stop
    return values


@dataclass(frozen=True)
class FrequencyGrid:
    start: float  # [Hz]
    stop: float  # [Hz]
    step: float  # [Hz]; irrelevant (but must be > 0) when start == stop

    def __post_init__(self):
        try:
            start = validate_frequency(self.start)
            stop = validate_frequency(self.stop)
            step = finite_real("frequency.step", self.step)
        except ValueError as e:
            raise SweepError(str(e)) from None
        object.__setattr__(self, "start", start)
        object.__setattr__(self, "stop", stop)
        object.__setattr__(self, "step", step)
        inclusive_range(start, stop, step, "frequency")  # validates step > 0, multiple, size

    @classmethod
    def single(cls, frequency):
        f = validate_frequency(frequency)
        return cls(f, f, f)

    def frequencies(self):
        """1-D array of the grid frequencies [Hz]."""
        return inclusive_range(self.start, self.stop, self.step, "frequency")

    @property
    def count(self):
        return int(round((self.stop - self.start) / self.step)) + 1


@dataclass(frozen=True)
class ParameterAxis:
    """One swept geometry parameter: an ordered tuple of distinct finite values."""

    name: str
    values: tuple

    def __post_init__(self):
        if self.name not in PARAMETER_NAMES:
            raise SweepError(f"unknown geometry parameter {self.name!r}; sweepable: {list(PARAMETER_NAMES)}")
        try:
            vals = tuple(finite_real(f"{self.name} value", v) for v in self.values)
        except (ValueError, TypeError) as e:
            raise SweepError(str(e)) from None
        if not vals:
            raise SweepError(f"parameter {self.name!r} has no values")
        if len(vals) > MAX_POINTS:
            raise SweepError(f"parameter {self.name!r}: {len(vals)} values exceeds {MAX_POINTS}")
        if len(set(vals)) != len(vals):
            raise SweepError(f"parameter {self.name!r} has duplicate values")
        object.__setattr__(self, "values", vals)

    @classmethod
    def from_range(cls, name, start, stop, step):
        return cls(name, tuple(inclusive_range(start, stop, step, name)))


# ------------------------------------------------------------------------------- policies


@dataclass(frozen=True)
class RefinementPolicy:
    """Opt-in refinement study (discretization.py). ``frequencies``: 'max' = highest sweep
    frequency only (the smallest wavelength for a given mesh), 'all' = every frequency."""

    enabled: bool = False
    factor: int = DEFAULT_REFINEMENT_FACTOR
    tolerance: float = DEFAULT_REFINEMENT_TOLERANCE
    frequencies: str = "max"

    def __post_init__(self):
        if not isinstance(self.enabled, bool):
            raise SweepError(f"refinement.enabled must be a boolean, got {self.enabled!r}")
        if isinstance(self.factor, bool) or not isinstance(self.factor, (int, np.integer)) or self.factor < 2:
            raise SweepError(f"refinement.factor must be an integer >= 2, got {self.factor!r}")
        try:
            tol = finite_real("refinement.tolerance", self.tolerance)
        except ValueError as e:
            raise SweepError(str(e)) from None
        if tol <= 0.0:
            raise SweepError(f"refinement.tolerance must be > 0, got {tol}")
        object.__setattr__(self, "tolerance", tol)
        object.__setattr__(self, "factor", int(self.factor))
        if self.frequencies not in REFINEMENT_FREQUENCIES:
            raise SweepError(
                f"refinement.frequencies must be one of {list(REFINEMENT_FREQUENCIES)}, got {self.frequencies!r}"
            )


@dataclass(frozen=True)
class ValidationPolicy:
    """Thresholds of the a priori checks; ``strict`` turns a failed check into an error before
    any solve. Without ``strict`` failures are recorded and the solve proceeds."""

    min_cells_per_wavelength: float = DEFAULT_MIN_CELLS_PER_WAVELENGTH
    max_area_error: float = DEFAULT_MAX_AREA_ERROR
    strict: bool = False
    refinement: RefinementPolicy = field(default_factory=RefinementPolicy)

    def __post_init__(self):
        if not isinstance(self.strict, bool):
            raise SweepError(f"validation.strict must be a boolean, got {self.strict!r}")
        try:
            object.__setattr__(
                self,
                "min_cells_per_wavelength",
                finite_real("min_cells_per_wavelength", self.min_cells_per_wavelength),
            )
            object.__setattr__(self, "max_area_error", finite_real("max_area_error", self.max_area_error))
        except ValueError as e:
            raise SweepError(str(e)) from None
        if self.min_cells_per_wavelength <= 0.0 or self.max_area_error < 0.0:
            raise SweepError("min_cells_per_wavelength must be > 0 and max_area_error >= 0")

    def as_dict(self):
        r = self.refinement
        return {
            "min_cells_per_wavelength": self.min_cells_per_wavelength,
            "max_area_error": self.max_area_error,
            "strict": self.strict,
            "refinement": {
                "enabled": r.enabled,
                "factor": r.factor,
                "tolerance": r.tolerance,
                "frequencies": r.frequencies,
            },
        }


@dataclass(frozen=True)
class SweepDefinition:
    frequency: FrequencyGrid | None = None  # None: the single excitation.frequency
    axes: tuple = ()
    on_invalid: str = "error"
    validation: ValidationPolicy = field(default_factory=ValidationPolicy)

    def __post_init__(self):
        object.__setattr__(self, "axes", tuple(self.axes))
        names = [a.name for a in self.axes]
        if len(set(names)) != len(names):
            raise SweepError(f"duplicate sweep parameter(s) in {names}")
        if self.on_invalid not in ON_INVALID:
            raise SweepError(f"on_invalid must be one of {list(ON_INVALID)}, got {self.on_invalid!r}")
        size = math.prod(len(a.values) for a in self.axes)
        if size > MAX_CELLS:
            raise SweepError(f"{size} cells exceeds the limit of {MAX_CELLS}")


# ------------------------------------------------------------------ cells, context, planning


def base_parameters(geometry):
    """The sweepable parameters of a `Geometry` (five) or an `AnnularSlotPatch` (seven; PHASE15
    Part E) -- I/O units: m and degrees."""
    if isinstance(geometry, AnnularSlotPatch):
        return {
            "a_m": geometry.a, "b_m": geometry.b, "R_m": geometry.R, "r_m": geometry.r,
            "center_x_m": geometry.center_x, "center_y_m": geometry.center_y,
            "rotation_deg": units.to_deg(geometry.rotation),
        }
    return {
        "length_m": geometry.length,
        "width_m": geometry.width,
        "center_x_m": geometry.center_x,
        "center_y_m": geometry.center_y,
        "rotation_deg": units.to_deg(geometry.rotation),
    }


def geometry_from_parameters(kind, parameters):
    """Inverse of ``base_parameters`` (the only place degrees become radians for a sweep). PHASE15
    Part E: `kind == AnnularSlotPatch.kind` ("patch_annular_slot") builds an `AnnularSlotPatch`
    instead of a `Geometry` -- the swept parameter NAMES (`a_m`/`b_m`/`R_m`/`r_m` vs `length_m`/
    `width_m`) are what select the class, since `sweep_context`'s own `geometry_kind` already
    records which family a given sweep is for."""
    if kind == AnnularSlotPatch.kind:
        missing = set(_ANNULAR_SLOT_PARAMETER_NAMES) - set(parameters)
        if missing:
            raise ValueError(f"missing geometry parameter(s) {sorted(missing)}")
        return AnnularSlotPatch(
            a=parameters["a_m"], b=parameters["b_m"], R=parameters["R_m"], r=parameters["r_m"],
            center_x=parameters["center_x_m"], center_y=parameters["center_y_m"],
            rotation=units.from_deg(parameters["rotation_deg"]),
        )
    base_names = ("length_m", "width_m", "center_x_m", "center_y_m", "rotation_deg")
    missing = set(base_names) - set(parameters)
    if missing:
        raise ValueError(f"missing geometry parameter(s) {sorted(missing)}")
    return Geometry(
        kind=kind,
        length=parameters["length_m"],
        width=parameters["width_m"],
        center_x=parameters["center_x_m"],
        center_y=parameters["center_y_m"],
        rotation=units.from_deg(parameters["rotation_deg"]),
    )


def sweep_context(config):
    """Everything that defines a database except the swept parameters and frequencies."""
    u, m, rp = config.unit_cell, config.mesh, config.reference_planes
    return {
        "code_version": __version__,
        "geometry_kind": u.geometry.kind,
        "lattice": {"period_x_m": u.lattice.period_x, "period_y_m": u.lattice.period_y},
        "substrate": {
            "name": u.substrate.name,
            "eps_r": u.substrate.eps_r,
            "tan_delta": u.substrate.tan_delta,
            "thickness_m": u.substrate.thickness,
        },
        "metal": {"model": u.metal.model},
        # Phase 12A: which stack this was solved for. Included in the cell_id hash (below), so a
        # grounded and an ungrounded sweep of otherwise-identical geometry are never confused as
        # "the same cell" by the checkpoint/database machinery.
        "backing": u.backing,
        "mesh": {"cells_x": m.cells_x, "cells_y": m.cells_y, "floquet_order": m.floquet_order},
        "reference_planes": {
            "top_air_thickness_m": rp.top_air_thickness,
            "bottom_air_thickness_m": rp.bottom_air_thickness,
        },
    }


def unit_cell_from_context(context, parameters):
    """Rebuild the UnitCell that a database cell was solved for (context + parameters only).

    ``context.get("backing", "air")``: older contexts (written before Phase 12A) have no
    'backing' key at all; they were always the ungrounded ("air") stack, so absence means "air",
    exactly like the config-file default (CONVENTIONS.md).
    """
    from .materials import Metal, Substrate
    from .periodic_bc import RectangularLattice

    lat, sub = context["lattice"], context["substrate"]
    return UnitCell(
        lattice=RectangularLattice(lat["period_x_m"], lat["period_y_m"]),
        geometry=geometry_from_parameters(context["geometry_kind"], parameters),
        substrate=Substrate(
            name=sub["name"],
            eps_r=sub["eps_r"],
            tan_delta=sub["tan_delta"],
            thickness=sub["thickness_m"],
        ),
        metal=Metal(**context["metal"]),
        backing=context.get("backing", "air"),
    )


@dataclass(frozen=True, eq=False)
class SweepCell:
    index: int
    cell_id: str
    parameters: dict
    status: str  # 'ok' | 'invalid'
    reason: str | None = None  # why it is invalid
    unit_cell: UnitCell | None = None
    checks: dict | None = None  # frequency-independent a priori checks (discretization.check_cell)


@dataclass(frozen=True, eq=False)
class SweepPlan:
    config: object  # the base SimulationConfig
    definition: SweepDefinition
    frequencies: np.ndarray
    context: dict
    cells: tuple
    # Phase 12B: tuple of tuples of cell ``index`` sharing an identical discretized metal mask on
    # config.mesh (discretization.duplicate_geometry_groups) -- () for every pre-Phase-12B sweep
    # construction (the field has a default so nothing else changes).
    mesh_duplicate_groups: tuple = ()

    @property
    def valid_cells(self):
        return tuple(c for c in self.cells if c.status == "ok")

    @property
    def n_solves(self):
        return len(self.valid_cells) * len(self.frequencies)

    @property
    def n_mesh_distinct(self):
        """Number of OK cells that are geometrically distinct *on this mesh*: duplicates within
        each group of size k contribute k-1 "extra" cells that discretize to the same problem."""
        extra = sum(len(g) - 1 for g in self.mesh_duplicate_groups)
        return len(self.valid_cells) - extra


def _describe(index, parameters):
    return f"cell {index} " + ", ".join(f"{k}={v:.6g}" for k, v in parameters.items())


def plan_sweep(config, definition=None, *, warn_duplicates=True):
    """Expand ``definition`` (default ``config.sweep``) into validated cells; see module docstring.

    ``warn_duplicates`` (Phase 12B, default True for a direct/interactive call): issue a
    UserWarning if any two OK cells discretize to an identical metal mask on ``config.mesh``
    (``discretization.duplicate_geometry_groups``) -- they will solve to bit-identical S-parameters
    despite differing continuous geometry. The group data itself (``SweepPlan.mesh_duplicate_groups``)
    is always computed and returned regardless of this flag ("warning + record flag", Phase 12B
    brief). ``runner.run_sweep`` passes ``warn_duplicates=False`` (it records the condition
    per-plan instead), the same convention as ``compute_s_parameters``'s ``warn_grating``/
    ``warn_surface_wave``, so no pre-Phase-12B sweep gains a new warning.
    """
    definition = definition if definition is not None else getattr(config, "sweep", None)
    if definition is None:
        raise SweepError("no sweep definition: pass one or add a 'sweep:' section to the config")
    policy = definition.validation
    grid = definition.frequency or FrequencyGrid.single(config.excitation.frequency)
    frequencies = grid.frequencies()
    context = sweep_context(config)
    base = base_parameters(config.unit_cell.geometry)
    names = [a.name for a in definition.axes]

    cells, invalid = [], []
    for index, combo in enumerate(itertools.product(*[a.values for a in definition.axes])):
        parameters = dict(base)
        parameters.update(zip(names, combo))
        cid = make_cell_id(context, parameters)
        try:
            unit_cell = unit_cell_from_context(context, parameters)
            checks = check_cell(unit_cell, config.mesh, policy.max_area_error)
        except ValueError as e:
            invalid.append((index, parameters, str(e)))
            cells.append(SweepCell(index, cid, parameters, "invalid", str(e)))
            continue
        cells.append(SweepCell(index, cid, parameters, "ok", None, unit_cell, checks))

    if invalid and definition.on_invalid == "error":
        shown = "; ".join(f"{_describe(i, p)}: {why}" for i, p, why in invalid[:5])
        more = f" (+{len(invalid) - 5} more)" if len(invalid) > 5 else ""
        raise SweepError(
            f"{len(invalid)} of {len(cells)} sweep cell(s) are invalid: {shown}{more}. "
            f"Fix the ranges or set on_invalid: skip to record them as skipped."
        )

    if policy.strict:
        failures = []
        for cell in cells:
            if cell.status != "ok":
                continue
            if not cell.checks["area_ok"]:
                failures.append(
                    f"{_describe(cell.index, cell.parameters)}: metal area error "
                    f"{cell.checks['metal_area_error']:.3g} > {policy.max_area_error:g}"
                )
            for f in frequencies:
                chk = check_frequency(cell.unit_cell, config.mesh, f, policy.min_cells_per_wavelength)
                if not chk["wavelength_ok"]:
                    failures.append(
                        f"{_describe(cell.index, cell.parameters)} at {f:.6g} Hz: "
                        f"{chk['cells_per_wavelength']:.3g} cells/wavelength < "
                        f"{policy.min_cells_per_wavelength:g}"
                    )
                    break
                if not chk["specular_only"]:
                    failures.append(
                        f"{_describe(cell.index, cell.parameters)} at {f:.6g} Hz: "
                        f"a grating order propagates (period >= wavelength)"
                    )
                    break
        if failures:
            more = f" (+{len(failures) - 5} more)" if len(failures) > 5 else ""
            raise SweepError(f"strict validation failed for {len(failures)} check(s): {'; '.join(failures[:5])}{more}")

    duplicate_groups = duplicate_geometry_groups(
        (cell.index, metal_mask_key(cell.unit_cell, config.mesh)) for cell in cells if cell.status == "ok"
    )
    if warn_duplicates and duplicate_groups:
        n_extra = sum(len(g) - 1 for g in duplicate_groups)
        warnings.warn(
            f"{len(duplicate_groups)} group(s) of mesh-identical geometries on this mesh "
            f"(cells_x={config.mesh.cells_x}, cells_y={config.mesh.cells_y}): {n_extra} of "
            f"{sum(1 for c in cells if c.status == 'ok')} OK cell(s) discretize to the exact same "
            f"metal mask as another cell in this sweep, and will solve to bit-identical "
            f"S-parameters despite differing continuous geometry (discretization.py)",
            UserWarning,
            stacklevel=2,
        )

    return SweepPlan(config, definition, frequencies, context, tuple(cells), duplicate_groups)
