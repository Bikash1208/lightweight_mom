"""Configuration loading and validation.

Config values are SI unless a unit suffix is given (see units.py). Unknown keys
are rejected so typos cannot be silently ignored, and no electrical value has a
built-in default: substrate parameters must be in the file.

Note: PyYAML (YAML 1.1) parses ``8e-3`` as a *string* (a float needs a dot:
``8.0e-3``). Bare-number strings are therefore accepted as SI numbers.
"""

import math
from dataclasses import dataclass
from pathlib import Path

import yaml

from . import units
from .basis_functions import MeshSettings
from .excitation import PlaneWaveExcitation, ReferencePlanes
from .geometry import Geometry, UnitCell
from .materials import Metal, Substrate
from .periodic_bc import RectangularLattice
from .sweep import (
    FrequencyGrid,
    ParameterAxis,
    RefinementPolicy,
    SweepDefinition,
    ValidationPolicy,
)

# Lengths >= this are almost certainly a unit mistake (e.g. "8" meant as mm).
MAX_PLAUSIBLE_LENGTH = 1.0  # [m]


class ConfigError(ValueError):
    """Invalid or unsupported configuration."""


@dataclass(frozen=True)
class SimulationConfig:
    unit_cell: UnitCell
    excitation: PlaneWaveExcitation
    mesh: MeshSettings
    reference_planes: ReferencePlanes
    sweep: SweepDefinition | None = None  # optional 'sweep:' section (Phase 5)


def _mapping(data, name):
    if not isinstance(data, dict):
        raise ConfigError(f"'{name}' must be a mapping, got {type(data).__name__}")
    return data


def _check_keys(section, data, required, optional=()):
    missing = set(required) - set(data)
    if missing:
        raise ConfigError(f"'{section}': missing required key(s) {sorted(missing)}")
    unknown = set(data) - set(required) - set(optional)
    if unknown:
        raise ConfigError(
            f"'{section}': unknown key(s) {sorted(unknown)}; allowed: {sorted(set(required) | set(optional))}"
        )


def _length(section, key, value):
    what = f"{section}.{key}"
    try:
        v = units.parse_length(value, what)
    except ValueError as e:
        raise ConfigError(str(e)) from None
    if abs(v) >= MAX_PLAUSIBLE_LENGTH:
        raise ConfigError(
            f"{what} = {v} m is >= {MAX_PLAUSIBLE_LENGTH} m; plain numbers are metres (write '8 mm' or 0.008, not 8)"
        )
    return v


def _wrap(section, fn, *args, **kwargs):
    """Call a validating constructor, re-raising ValueError as ConfigError."""
    try:
        return fn(*args, **kwargs)
    except ValueError as e:
        raise ConfigError(f"'{section}': {e}") from None


def config_from_dict(data):
    """Build a validated SimulationConfig from an already-parsed mapping."""
    data = _mapping(data, "<root>")
    _check_keys(
        "<root>",
        data,
        required=(
            "period_x",
            "period_y",
            "substrate",
            "geometry",
            "excitation",
            "mesh",
            "reference_planes",
        ),
        optional=("metal", "sweep"),
    )

    # --- lattice -----------------------------------------------------------
    lattice = _wrap(
        "lattice",
        RectangularLattice,
        _length("<root>", "period_x", data["period_x"]),
        _length("<root>", "period_y", data["period_y"]),
    )

    # --- substrate ---------------------------------------------------------
    sub = _mapping(data["substrate"], "substrate")
    _check_keys("substrate", sub, required=("eps_r", "tan_delta", "thickness"), optional=("name", "backing"))
    substrate = _wrap(
        "substrate",
        Substrate,
        name=str(sub.get("name", "unnamed")),
        eps_r=sub["eps_r"],
        tan_delta=sub["tan_delta"],
        thickness=_length("substrate", "thickness", sub["thickness"]),
    )
    # Phase 12A: an optional PEC ground plane at z = -h. Absence means "air" (the default,
    # documented in CONVENTIONS.md): every pre-Phase-12A config has no 'backing' key and loads
    # with identical behaviour. Kept in the 'substrate' section (Sec. 1's own suggestion) rather
    # than a new top-level section, since it is a property of the stack's bottom termination, the
    # same place substrate.thickness (the stack's other z-extent) already lives.
    backing = sub.get("backing", "air")
    if backing not in ("air", "pec"):
        raise ConfigError(f"'substrate.backing' must be 'air' or 'pec', got {backing!r}")

    # --- metal -------------------------------------------------------------
    met = _mapping(data.get("metal", {"model": "pec"}), "metal")
    _check_keys("metal", met, required=(), optional=("model",))
    metal = _wrap("metal", Metal, **met)

    # --- geometry ----------------------------------------------------------
    geo = _mapping(data["geometry"], "geometry")
    _check_keys(
        "geometry",
        geo,
        required=("type", "length", "width"),
        optional=("center_x", "center_y", "rotation_deg"),
    )
    geometry = _wrap(
        "geometry",
        Geometry,
        kind=geo["type"],
        length=_length("geometry", "length", geo["length"]),
        width=_length("geometry", "width", geo["width"]),
        center_x=_length("geometry", "center_x", geo.get("center_x", 0.0)),
        center_y=_length("geometry", "center_y", geo.get("center_y", 0.0)),
        rotation=units.from_deg(_number("geometry.rotation_deg", geo.get("rotation_deg", 0.0))),
    )

    unit_cell = _wrap("unit_cell", UnitCell, lattice, geometry, substrate, metal, backing=backing)

    # --- excitation --------------------------------------------------------
    exc = _mapping(data["excitation"], "excitation")
    _check_keys(
        "excitation",
        exc,
        required=("frequency",),
        optional=("theta_deg", "phi_deg", "polarization"),
    )
    theta_deg = _number("excitation.theta_deg", exc.get("theta_deg", 0.0))
    if theta_deg != 0.0:
        raise ConfigError("excitation.theta_deg: Phase 1 supports normal incidence only (theta_deg must be 0)")
    try:
        frequency = units.parse_frequency(exc["frequency"], "excitation.frequency")
    except ValueError as e:
        raise ConfigError(str(e)) from None
    excitation = _wrap(
        "excitation",
        PlaneWaveExcitation,
        frequency=frequency,
        theta=units.from_deg(theta_deg),
        phi=units.from_deg(_number("excitation.phi_deg", exc.get("phi_deg", 0.0))),
        polarization=exc.get("polarization", "x"),
    )

    # --- discretization ----------------------------------------------------
    msh = _mapping(data["mesh"], "mesh")
    _check_keys("mesh", msh, required=("cells_x", "cells_y", "floquet_order"))
    mesh = _wrap("mesh", MeshSettings, **msh)

    # --- S-parameter reference planes (explicit; not used by Z or V) ---------
    ref = _mapping(data["reference_planes"], "reference_planes")
    if backing == "air":
        _check_keys("reference_planes", ref, required=("top_air_thickness", "bottom_air_thickness"))
        bottom_air_thickness = _length("reference_planes", "bottom_air_thickness", ref["bottom_air_thickness"])
    else:  # "pec": there is no bottom port, so bottom_air_thickness is meaningless (Sec. 1)
        _check_keys("reference_planes", ref, required=("top_air_thickness",))
        bottom_air_thickness = None
    reference_planes = _wrap(
        "reference_planes",
        ReferencePlanes,
        top_air_thickness=_length("reference_planes", "top_air_thickness", ref["top_air_thickness"]),
        bottom_air_thickness=bottom_air_thickness,
        substrate_thickness=substrate.thickness,
    )

    sweep = _parse_sweep(data["sweep"]) if "sweep" in data else None

    return SimulationConfig(
        unit_cell=unit_cell,
        excitation=excitation,
        mesh=mesh,
        reference_planes=reference_planes,
        sweep=sweep,
    )


# YAML parameter name -> (internal name, kind). Lengths accept unit suffixes; rotation is degrees.
_SWEEP_PARAMETERS = {
    "length": ("length_m", "length"),
    "width": ("width_m", "length"),
    "center_x": ("center_x_m", "length"),
    "center_y": ("center_y_m", "length"),
    "rotation_deg": ("rotation_deg", "degrees"),
}


def _sweep_value(section, key, kind, value):
    return _length(section, key, value) if kind == "length" else _number(f"{section}.{key}", value)


def _parse_axis(yaml_name, spec):
    if yaml_name not in _SWEEP_PARAMETERS:
        raise ConfigError(
            f"'sweep.parameters': unknown geometry parameter {yaml_name!r}; sweepable: {sorted(_SWEEP_PARAMETERS)}"
        )
    name, kind = _SWEEP_PARAMETERS[yaml_name]
    section = f"sweep.parameters.{yaml_name}"
    if isinstance(spec, list):
        spec = {"values": spec}
    spec = _mapping(spec, section)
    has_values, has_range = "values" in spec, bool({"start", "stop", "step"} & set(spec))
    if has_values == has_range:
        raise ConfigError(f"'{section}': give either 'values' or 'start'/'stop'/'step', not both or neither")
    if has_values:
        _check_keys(section, spec, required=("values",))
        if not isinstance(spec["values"], list):
            raise ConfigError(f"'{section}.values' must be a list")
        return _wrap(
            section, ParameterAxis, name, tuple(_sweep_value(section, "values", kind, v) for v in spec["values"])
        )
    _check_keys(section, spec, required=("start", "stop", "step"))
    return _wrap(
        section,
        ParameterAxis.from_range,
        name,
        *(_sweep_value(section, k, kind, spec[k]) for k in ("start", "stop", "step")),
    )


def _parse_sweep(data):
    sw = _mapping(data, "sweep")
    _check_keys("sweep", sw, required=(), optional=("frequency", "parameters", "on_invalid", "validation"))

    frequency = None
    if "frequency" in sw:
        fr = _mapping(sw["frequency"], "sweep.frequency")
        _check_keys("sweep.frequency", fr, required=("start", "stop", "step"))
        try:
            values = [units.parse_frequency(fr[k], f"sweep.frequency.{k}") for k in ("start", "stop", "step")]
        except ValueError as e:
            raise ConfigError(str(e)) from None
        frequency = _wrap("sweep.frequency", FrequencyGrid, *values)

    params = _mapping(sw.get("parameters", {}), "sweep.parameters")
    axes = tuple(_parse_axis(name, spec) for name, spec in params.items())

    kwargs = {}
    if "validation" in sw:
        val = _mapping(sw["validation"], "sweep.validation")
        _check_keys(
            "sweep.validation",
            val,
            required=(),
            optional=("min_cells_per_wavelength", "max_area_error", "strict", "refinement"),
        )
        kwargs = {k: v for k, v in val.items() if k != "refinement"}
        if "refinement" in val:
            ref = _mapping(val["refinement"], "sweep.validation.refinement")
            _check_keys(
                "sweep.validation.refinement",
                ref,
                required=(),
                optional=("enabled", "factor", "tolerance", "frequencies"),
            )
            kwargs["refinement"] = _wrap("sweep.validation.refinement", RefinementPolicy, **ref)
    validation = _wrap("sweep.validation", ValidationPolicy, **kwargs)

    extra = {"on_invalid": sw["on_invalid"]} if "on_invalid" in sw else {}
    return _wrap("sweep", SweepDefinition, frequency=frequency, axes=axes, validation=validation, **extra)


def _number(what, value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigError(f"{what}: expected a finite number, got {value!r}")
    return float(value)


def load_config(path):
    """Read a YAML file and return a validated SimulationConfig."""
    path = Path(path)
    try:
        with path.open("r", encoding="utf-8") as fh:
            data = yaml.safe_load(fh)
    except yaml.YAMLError as e:
        raise ConfigError(f"{path}: invalid YAML: {e}") from None
    return config_from_dict(data)
