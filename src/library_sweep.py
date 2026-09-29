"""Generalizes `library_builder.py`/`unit_cell_library.py`'s locked-to-Phase-7 L/W x frequency
sweep to an *arbitrary* `config.SimulationConfig` -- any period, substrate, or backing, not only
the one frozen ungrounded spec `library_schema.LibraryConfig` represents. This is what the app's
Library page needs: a grounded design point (e.g. Phase 12C's, `configs/config_phase12c_grounded_8mm.yaml`)
has its own phase-vs-size sweep to run, and CLAUDE.md is explicit that `LibraryConfig` deliberately
has no `backing` field ("Only one library is actually built, still ungrounded").

Reuses `runner.run_sweep` (Phase 5, unmodified) exactly as `library_builder.build_library_database`
does, and `database.Record.component`/`.phase_deg` (unmodified) exactly as
`unit_cell_library._s_entry` does. The only new logic is choosing which S-parameter families exist
to serialize: a grounded (`port_topology="one_port"`) record only has S11 (S12/S21/S22 are
undefined, not absent-by-omission -- `sparameters.SParameters`'s own docstring), so
`unit_cell_library.write_library_json`'s unconditional 4-family loop would raise on one; this
module checks `record.port_topology` first instead.
"""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .runner import run_sweep
from .sweep import FrequencyGrid, ParameterAxis, SweepDefinition, ValidationPolicy
from .validation_export import _db

LIBRARY_FORMAT_NAME = "lightweight_mom.app_unit_cell_library"
LIBRARY_FORMAT_VERSION = 1

_S_FAMILIES_TWO_PORT = ("S11", "S21", "S12", "S22")
_S_FAMILIES_ONE_PORT = ("S11",)


@dataclass(frozen=True)
class SweepSpec:
    """Just enough to drive a sweep and to duck-type as `phase_extraction`'s `config` argument
    (`f_design_hz`, `l_values_mm`, `w_values_mm`), without `library_schema.LibraryConfig`'s
    ungrounded-only `base_config_dict()`."""

    f_min_hz: float
    f_max_hz: float
    f_step_hz: float
    f_design_hz: float
    l_values_mm: tuple
    w_values_mm: tuple


def build_sweep_definition(spec: SweepSpec):
    l_values = tuple(v * 1e-3 for v in spec.l_values_mm)
    w_values = tuple(v * 1e-3 for v in spec.w_values_mm)
    frequency = FrequencyGrid(spec.f_min_hz, spec.f_max_hz, spec.f_step_hz)
    axes = (ParameterAxis("length_m", l_values), ParameterAxis("width_m", w_values))
    return SweepDefinition(frequency=frequency, axes=axes, on_invalid="skip", validation=ValidationPolicy())


def run_custom_sweep(sim_config, spec: SweepSpec, out_dir=None, *, resume=True, progress=None):
    """``sim_config``: an already-built `config.SimulationConfig` (any period/substrate/backing --
    the caller decides via `config.config_from_dict`, unmodified). Returns the same
    `runner.SweepResult` `library_builder.build_library_database` returns."""
    definition = build_sweep_definition(spec)
    return run_sweep(sim_config, out_dir, definition=definition, resume=resume, progress=progress)


def _s_entry(record, name):
    z = record.component(name)
    mag = abs(z)
    return {
        "real": float(z.real),
        "imag": float(z.imag),
        "mag": float(mag),
        "mag_db": float(_db(np.array([mag]))[0]),
        "phase_deg": float(record.phase_deg(name)),
    }


def _cell_doc(cell, database):
    l_mm = round(cell.parameters["length_m"] * 1e3, 9)
    w_mm = round(cell.parameters["width_m"] * 1e3, 9)
    geometry_status = "VALID" if cell.status == "ok" else "INVALID"
    if cell.status != "ok":
        return {
            "cell_id": cell.cell_id,
            "L_mm": l_mm,
            "W_mm": w_mm,
            "geometry_status": geometry_status,
            "solver_status": "FAILED",
            "invalid_reason": cell.reason,
            "frequency_GHz": [],
            "S11": [],
        }
    records = sorted(database.records_of(cell.index), key=lambda r: r.frequency_index)
    families = _S_FAMILIES_ONE_PORT if records[0].port_topology == "one_port" else _S_FAMILIES_TWO_PORT
    ctx = database.metadata["context"]
    doc = {
        "cell_id": cell.cell_id,
        "L_mm": l_mm,
        "W_mm": w_mm,
        "period_x_mm": ctx["lattice"]["period_x_m"] * 1e3,
        "period_y_mm": ctx["lattice"]["period_y_m"] * 1e3,
        "substrate_thickness_mm": ctx["substrate"]["thickness_m"] * 1e3,
        "epsilon_r": ctx["substrate"]["eps_r"],
        "tan_delta": ctx["substrate"]["tan_delta"],
        "backing": ctx["backing"],
        "geometry_status": geometry_status,
        "solver_status": "OK",
        "frequency_GHz": [r.frequency / 1e9 for r in records],
    }
    for fam in families:
        doc[fam] = [_s_entry(r, f"{fam}_xx") for r in records]
    return doc


@dataclass(frozen=True)
class DesignFrequencyRow:
    """Same shape as `phase_extraction.DesignFrequencyRow` (so `phase_extraction.phase_vs_geometry_table`
    and `library_classifier.classify_cell`, which only read `.l_mm`/`.w_mm`/`.s11`, accept it
    unchanged) -- but built without `phase_extraction.extract_design_frequency`'s unconditional
    S21/S12/S22 lookup, which raises for a grounded (one-port) record."""

    cell_id: str
    l_mm: float
    w_mm: float
    geometry_status: str
    solver_status: str
    s11: "dict | None"
    s21: "dict | None"
    s12: "dict | None"
    s22: "dict | None"


def extract_design_frequency(database, spec: SweepSpec):
    """Grounded-safe equivalent of `phase_extraction.extract_design_frequency`, for the same
    `f_design_hz`/`l_values_mm`/`w_values_mm` a `SweepSpec` already carries."""
    from .phase_extraction import design_frequency_index

    idx = design_frequency_index(database, spec.f_design_hz)
    rows = []
    for l_mm in spec.l_values_mm:
        for w_mm in spec.w_values_mm:
            cell = next(
                (
                    c
                    for c in database.cells
                    if abs(c.parameters["length_m"] - l_mm * 1e-3) < 1e-15
                    and abs(c.parameters["width_m"] - w_mm * 1e-3) < 1e-15
                ),
                None,
            )
            if cell is None:
                raise ValueError(f"no cell for L={l_mm} mm, W={w_mm} mm in this database -- incomplete sweep")
            if cell.status != "ok":
                rows.append(DesignFrequencyRow(cell.cell_id, l_mm, w_mm, "INVALID", "FAILED", None, None, None, None))
                continue
            record = next((r for r in database.records_of(cell.index) if r.frequency_index == idx), None)
            if record is None:
                rows.append(DesignFrequencyRow(cell.cell_id, l_mm, w_mm, "VALID", "FAILED", None, None, None, None))
                continue
            two_port = record.port_topology == "two_port"
            rows.append(
                DesignFrequencyRow(
                    cell.cell_id,
                    l_mm,
                    w_mm,
                    "VALID",
                    "OK",
                    _s_entry(record, "S11_xx"),
                    _s_entry(record, "S21_xx") if two_port else None,
                    _s_entry(record, "S12_xx") if two_port else None,
                    _s_entry(record, "S22_xx") if two_port else None,
                )
            )
    return rows


def write_app_library_json(database, path):
    """Same shape `unit_cell_library.UnitCellLibrary`'s query API reads (`cell["frequency_GHz"]`,
    `cell[family][i]`), except a grounded cell's document only ever has an `"S11"` key -- callers
    (e.g. `UnitCellLibrary.find_phase_candidates`) that default to `family="S11"` work unchanged;
    anything that asks a grounded library for S21 gets a normal `KeyError` on a missing dict key,
    which is the correct "undefined, not zero" behaviour for this format too."""
    doc = {
        "format": LIBRARY_FORMAT_NAME,
        "version": LIBRARY_FORMAT_VERSION,
        "cells": [_cell_doc(c, database) for c in database.cells],
    }
    path = Path(path)
    path.write_text(json.dumps(doc, indent=1, allow_nan=False) + "\n")
    return path
