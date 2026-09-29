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
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from .library_classifier import ResonanceResult, classify_cell
from .library_schema import RESONANCE_MIN_PROMINENCE
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


RESONANCE_PHASE_JUMP_MIN_DEG = 180.0  # a jump of exactly this much is an exact +-180deg crossing


def detect_resonance_phase_crossing(frequencies_hz, phase_deg, min_jump_deg=RESONANCE_PHASE_JUMP_MIN_DEG):
    """PHASE13B A6: `library_classifier.detect_resonance` looks for a magnitude PEAK -- correct
    for the ungrounded library it was written for, but measured on this project's own grounded
    design point (L=W=6.0mm, 8mm period, 12.5 GHz), |S11_xx| stays within 0.9957-0.9959 across
    12.0-13.0 GHz (prominence ~2e-4, far below `RESONANCE_MIN_PROMINENCE`): a low-loss grounded
    reflector barely absorbs anything even at resonance, so a magnitude-dip criterion can
    structurally never cross that threshold either -- it is not just "the wrong sign", it is the
    wrong physical quantity. The standard resonance signature for a low-loss reflective element
    (patch/reflectarray design curves, e.g. Huang & Encinar, "Reflectarray Antennas") is instead
    the reflection PHASE passing through +-180deg -- and that is exactly what the same measured
    sweep shows (-134.8deg at 12.0 GHz falling smoothly through -175.4deg at 12.6 GHz to +178.1deg
    at 12.7 GHz: a genuine crossing of the (-180,180] branch cut). Detected directly off the raw
    (wrapped) phase sequence: a jump of >= `min_jump_deg` between two adjacent frequency samples
    means the true continuous phase crossed an odd multiple of 180deg somewhere between them."""
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    phase_deg = np.asarray(phase_deg, dtype=float)
    criterion = (
        f"the raw ((-180,180]) reflection phase jumps by >= {min_jump_deg:g}deg between two "
        f"adjacent frequency samples -- the continuous phase crossed +-180deg between them (the "
        f"standard low-loss reflective-element resonance signature; no magnitude dip required, "
        f"unlike the ungrounded library's detect_resonance)"
    )
    if phase_deg.size < 2:
        return ResonanceResult(False, None, None, None, criterion)

    jumps = np.abs(np.diff(phase_deg))
    crossing_indices = np.flatnonzero(jumps >= min_jump_deg)
    if crossing_indices.size == 0:
        return ResonanceResult(False, None, None, None, criterion)

    i = int(crossing_indices[np.argmax(jumps[crossing_indices])])
    return ResonanceResult(
        is_resonant=True,
        resonance_frequency_hz=float(0.5 * (frequencies_hz[i] + frequencies_hz[i + 1])),
        resonance_metric=float(jumps[i]),
        peak_magnitude=None,  # not a magnitude criterion -- see docstring
        criterion=criterion,
    )


def classify_cell_backing_aware(cell_entry, database, design_row, unit_cell, min_prominence=RESONANCE_MIN_PROMINENCE):
    """Wraps `library_classifier.classify_cell` unmodified for an ungrounded cell; for a grounded
    one (PHASE13B A6), re-runs resonance detection with `detect_resonance_phase_crossing` instead
    -- every other flag (phase usability, grating risk) is backing-agnostic and left as
    `classify_cell` computed it."""
    classification = classify_cell(cell_entry, database, design_row, min_prominence)
    if unit_cell.backing != "pec" or classification.solver_status != "OK":
        return classification
    records = sorted(database.records_of(cell_entry.index), key=lambda r: r.frequency_index)
    freqs = np.array([r.frequency for r in records])
    phase = np.array([r.phase_deg("S11_xx") for r in records])
    resonance = detect_resonance_phase_crossing(freqs, phase)
    return replace(classification, resonance=resonance)


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


def build_square_library_cells(base_config_dict, l_values_mm, frequency_hz, progress=None):
    """PHASE13C: one real solve per L value, SQUARE cells only (L=W) -- the wizard's auto-built
    board library never exposes an independent L,W grid the way the Library page does. `base_
    config_dict`: `config.config_from_dict`'s own shape, with `geometry`/`excitation.frequency`
    as placeholders (overwritten per L here, same convention `routers/library.py`'s
    `LibrarySweepRequest` already documents). Deliberately one `run_custom_sweep` call per L
    (each a single-cell sweep) rather than `l_values_mm x w_values_mm` through
    `build_sweep_definition`'s two independent axes, which would solve the full cross product
    (e.g. 25 cells for 5 L values) when only the 5 diagonal (square) cells are wanted. Returns the
    list of per-cell JSON docs (the same shape `write_app_library_json`'s `doc["cells"]` has)."""
    from .config import config_from_dict

    cells = []
    for i, l_mm in enumerate(l_values_mm):
        if progress is not None:
            progress({"done": i, "total": len(l_values_mm)})
        cfg = config_from_dict(
            {
                **base_config_dict,
                "geometry": {"type": "rectangle", "length": l_mm * 1e-3, "width": l_mm * 1e-3},
                "excitation": {"frequency": frequency_hz},
            }
        )
        spec = SweepSpec(
            f_min_hz=frequency_hz,
            f_max_hz=frequency_hz,
            f_step_hz=1.0e9,
            f_design_hz=frequency_hz,
            l_values_mm=(l_mm,),
            w_values_mm=(l_mm,),
        )
        result = run_custom_sweep(cfg, spec, out_dir=None)
        cells.append(_cell_doc(result.database.cells[0], result.database))
    if progress is not None:
        progress({"done": len(l_values_mm), "total": len(l_values_mm)})
    return cells


def write_square_library_json(cells, path):
    """Writes a `build_square_library_cells` cell list in the same on-disk shape
    `write_app_library_json` produces."""
    doc = {"format": LIBRARY_FORMAT_NAME, "version": LIBRARY_FORMAT_VERSION, "cells": cells}
    path = Path(path)
    path.write_text(json.dumps(doc, indent=1, allow_nan=False) + "\n")
    return path
