"""Unit-cell database: in-memory model, JSON and CSV export/import.

One *record* = the 4x4 complex two-port S-matrix of one unit cell at one frequency, in the
linear x/y basis at the reference planes stored in the metadata (matrix layout, indexing and
sign conventions: S_PARAMETERS.md; circular basis and handedness: polarization.py, DATABASE.md).
Everything in a file is deterministic: no timestamps, no host names, no absolute paths.

JSON (archival, complete)
    Full metadata (conventions, lattice/substrate/metal/mesh/planes, sweep axes, frequency grid,
    validation policy), per-cell parameters/status/validation, per-record S. Each of the 16
    linear entries is stored as {re, im, mag, phase_deg}: ``re``/``im`` are the value (Python
    ``repr`` floats, so the round trip is bit-exact); ``mag``/``phase_deg`` are derived, checked
    against re/im on import. The 16 circular entries (mag, phase_deg) are derived likewise.
CSV (flat table, one row per record)
    Columns: identifiers, the five geometry parameters (SI/deg), diagnostics, then mag and
    phase_deg of every linear entry and of every circular entry. Only mag and phase_deg are
    stored, so the complex value read back differs from the JSON value by a few ulp; the mag and
    phase_deg columns themselves round-trip exactly (``repr`` floats). Metadata is not in the
    CSV: pass it to ``Database.from_csv(path, metadata=...)`` if needed. Skipped (invalid)
    cells appear only in the JSON.

Phase is the principal value in (-180, 180] degrees (phase.py).
"""

import csv
import json
import platform
from dataclasses import dataclass

import numpy as np

from ._canonical import cell_id as make_cell_id
from .phase import (
    STEP_WARN_DEG,
    max_wrapped_step_deg,
    principal_phase_deg,
    unwrap_deg,
    wrapped_difference_deg,
)
from .polarization import (
    CIRCULAR_COMPONENTS,
    LINEAR_COMPONENTS,
    block_to_circular,
    matrix4_to_circular,
)

FORMAT_NAME = "lightweight_mom.unit_cell_database"
FORMAT_VERSION = 1

FAMILIES = ("S11", "S12", "S21", "S22")
PARAMETER_NAMES = ("length_m", "width_m", "center_x_m", "center_y_m", "rotation_deg")

# name -> python type of every per-record diagnostic (strict schema for both file formats)
DIAGNOSTIC_TYPES = {
    "n_unknowns": int,
    "incoming_residual": float,
    "reciprocity_error": float,
    "passivity_min_eig": float,
    "cells_per_wavelength": float,
    "specular_only": bool,
    "checks_passed": bool,
}

CONVENTIONS = {
    "time_dependence": "exp(+j*omega*t)",
    "units": {
        "frequency": "Hz", "length": "m", "angle": "deg", "phase": "deg",
        "magnitude": "linear (|S|, not dB)",
    },
    "phase": "arg(S), principal value in (-180, 180] deg; a delay is a negative phase",
    "s_matrix": (
        "4x4 [[S11, S12], [S21, S22]] over (port1 x, port1 y, port2 x, port2 y); "
        "entry name S<out><in>_<received><incident>; port 1 = top reference plane (wave "
        "incident travelling -z), port 2 = bottom reference plane (wave incident travelling +z)"
    ),
    "linear_basis": "global x and y unit vectors, identical at both ports",
    "circular_basis": {
        "R": "(x + j y)/sqrt(2)",
        "L": "(x - j y)/sqrt(2)",
        "rotation_in_time": "R rotates clockwise seen from +z (x toward -y); L counter-clockwise",
        "handedness": (
            "right-hand rule, thumb along propagation: R is RHCP when travelling -z and LHCP "
            "when travelling +z; L is the reverse"
        ),
        "components": "S<out><in>_<received><incident> with labels R/L, e.g. S21_LR",
        "transform": "S_circ = T^H S_lin T, T = [e_R e_L] (unitary)",
    },
    "co_cross": {
        "linear": {"co": ["xx", "yy"], "cross": ["yx", "xy"]},
        "circular": {"co": ["RR", "LL"], "cross": ["LR", "RL"]},
    },
}


class DatabaseError(ValueError):
    """Malformed, inconsistent or unsupported database content."""


# ----------------------------------------------------------------------- S-matrix accessors

def _split_name(name):
    family, _, comp = name.partition("_")
    if family not in FAMILIES:
        raise ValueError(f"unknown S-parameter family in {name!r}; use one of {list(FAMILIES)}")
    return family, comp


def family_block(s4, family):
    """(2, 2) block [received pol, incident pol] of the 4x4 matrix for 'S11'/'S12'/'S21'/'S22'."""
    out, inc = int(family[1]) - 1, int(family[2]) - 1
    return np.asarray(s4)[2 * out : 2 * out + 2, 2 * inc : 2 * inc + 2]


def component_value(s4, name):
    """One complex entry by name: linear 'S21_yx' or circular 'S21_LR' ([received][incident]).

    ``s4`` may also be a (2, 2) one-port (grounded, Phase 12A) matrix -- S11 only, ``family``
    must then be 'S11' (any other family is undefined for a one-port record, not zero)."""
    family, comp = _split_name(name)
    s4 = np.asarray(s4)
    if s4.shape == (2, 2):
        if family != "S11":
            raise ValueError(f"{name}: only S11 is defined for a one-port (grounded) record")
        block = s4
    else:
        block = family_block(s4, family)
    if comp in LINEAR_COMPONENTS:
        return complex(block[LINEAR_COMPONENTS[comp]])
    if comp in CIRCULAR_COMPONENTS:
        return complex(block_to_circular(block)[CIRCULAR_COMPONENTS[comp]])
    raise ValueError(
        f"unknown component in {name!r}; linear {list(LINEAR_COMPONENTS)}, "
        f"circular {list(CIRCULAR_COMPONENTS)}"
    )


def component_names(basis):
    """All 16 entry names of 'linear' or 'circular', in file order (family-major)."""
    comps = {"linear": LINEAR_COMPONENTS, "circular": CIRCULAR_COMPONENTS}.get(basis)
    if comps is None:
        raise ValueError("basis must be 'linear' or 'circular'")
    return [f"{fam}_{c}" for fam in FAMILIES for c in comps]


def _csv_columns():
    cols = ["cell_id", "cell_index", "frequency_index", "frequency_hz"]
    cols += list(PARAMETER_NAMES) + list(DIAGNOSTIC_TYPES)
    for basis in ("linear", "circular"):
        for name in component_names(basis):
            cols += [f"{name}_mag", f"{name}_phase_deg"]
    return tuple(cols)


CSV_COLUMNS = _csv_columns()


# ------------------------------------------------------------------------------ data model

PORT_TOPOLOGIES = ("two_port", "one_port")


@dataclass(frozen=True, eq=False)
class Record:
    """S-matrix of one cell at one frequency.

    ``port_topology`` (Phase 12A): ``"two_port"`` (default; every pre-Phase-12A record) stores
    the full 4x4 two-port matrix, as always. ``"one_port"`` (a grounded stack, ``backing="pec"``:
    there is no port 2) stores just the (2, 2) S11 block -- S12/S21/S22 are undefined, not zero,
    and are simply absent, not stored as anything (SParameters, S_PARAMETERS.md).
    """

    cell_id: str
    frequency_index: int
    frequency: float          # [Hz]
    s: np.ndarray             # (4, 4) two_port, or (2, 2) one_port (S11 only): linear x/y basis
    diagnostics: dict
    port_topology: str = "two_port"

    def __post_init__(self):
        if self.port_topology not in PORT_TOPOLOGIES:
            raise DatabaseError(
                f"port_topology must be one of {PORT_TOPOLOGIES}, got {self.port_topology!r}"
            )
        s = np.array(self.s, dtype=complex)
        expected_shape = (2, 2) if self.port_topology == "one_port" else (4, 4)
        if s.shape != expected_shape or not np.all(np.isfinite(s)):
            raise DatabaseError(
                f"record S must be a finite {expected_shape} complex matrix "
                f"(port_topology={self.port_topology!r}), got shape {s.shape}"
            )
        s.setflags(write=False)
        object.__setattr__(self, "s", s)
        object.__setattr__(self, "frequency_index", int(self.frequency_index))
        object.__setattr__(self, "frequency", float(self.frequency))
        object.__setattr__(self, "diagnostics", _coerce_diagnostics(self.diagnostics))

    def component(self, name):
        return component_value(self.s, name)

    def magnitude(self, name):
        return abs(self.component(name))

    def phase_deg(self, name):
        return float(principal_phase_deg(self.component(name)))

    def same_as(self, other):
        return (
            self.cell_id == other.cell_id
            and self.frequency_index == other.frequency_index
            and self.frequency == other.frequency
            and self.port_topology == other.port_topology
            and np.array_equal(self.s, other.s)
            and self.diagnostics == other.diagnostics
        )


@dataclass(frozen=True)
class CellEntry:
    index: int
    cell_id: str
    parameters: dict           # PARAMETER_NAMES -> float
    status: str                # 'ok' | 'invalid'
    reason: str | None = None
    validation: dict | None = None      # cell-level checks and refinement result (free-form, JSON-able)


def _coerce_diagnostics(diag):
    if set(diag) != set(DIAGNOSTIC_TYPES):
        raise DatabaseError(
            f"diagnostics keys {sorted(diag)} differ from the schema {sorted(DIAGNOSTIC_TYPES)}"
        )
    out = {}
    for key, typ in DIAGNOSTIC_TYPES.items():
        v = diag[key]
        if typ is bool:
            if not isinstance(v, (bool, np.bool_)):
                raise DatabaseError(f"diagnostic {key!r} must be a boolean, got {v!r}")
            out[key] = bool(v)
        elif typ is int:
            if isinstance(v, bool) or not isinstance(v, (int, np.integer)):
                raise DatabaseError(f"diagnostic {key!r} must be an integer, got {v!r}")
            out[key] = int(v)
        else:
            if isinstance(v, bool) or not isinstance(v, (int, float, np.floating, np.integer)):
                raise DatabaseError(f"diagnostic {key!r} must be a number, got {v!r}")
            if not np.isfinite(v):
                raise DatabaseError(f"diagnostic {key!r} must be finite, got {v!r}")
            out[key] = float(v)
    return out


class Database:
    """Cells + records + metadata. Immutable by convention; build one from a sweep
    (runner.run_sweep) or import one (``from_json`` / ``from_csv``)."""

    def __init__(self, metadata, cells, records, environment=None):
        self.metadata = dict(metadata)
        self.cells = tuple(cells)
        self.environment = dict(environment) if environment is not None else default_environment()
        self._by_id, self._by_index = {}, {}
        for c in self.cells:
            if c.cell_id in self._by_id:
                raise DatabaseError(f"duplicate cell_id {c.cell_id!r}")
            if c.index in self._by_index:
                raise DatabaseError(f"duplicate cell index {c.index}")
            self._by_id[c.cell_id] = self._by_index[c.index] = c
        if [c.index for c in self.cells] != sorted(self._by_index):
            raise DatabaseError("cells must be listed in increasing index order")
        records = list(records)
        seen = set()
        for r in records:
            if r.cell_id not in self._by_id:
                raise DatabaseError(f"record refers to unknown cell_id {r.cell_id!r}")
            if self._by_id[r.cell_id].status != "ok":
                raise DatabaseError(f"record for cell {r.cell_id!r} whose status is not 'ok'")
            key = (r.cell_id, r.frequency_index)
            if key in seen:
                raise DatabaseError(f"duplicate record {key}")
            seen.add(key)
        self.records = tuple(sorted(records, key=lambda r: (self._by_id[r.cell_id].index, r.frequency_index)))

    # ------------------------------------------------------------------------ equality
    def __eq__(self, other):
        if not isinstance(other, Database):
            return NotImplemented
        return (
            self.metadata == other.metadata
            and self.cells == other.cells
            and len(self.records) == len(other.records)
            and all(a.same_as(b) for a, b in zip(self.records, other.records))
        )

    __hash__ = None

    # ------------------------------------------------------------------------- queries
    def cell(self, key):
        """A CellEntry by cell_id (str) or by its ``index`` (int)."""
        table = self._by_id if isinstance(key, str) else self._by_index
        try:
            return table[key]
        except KeyError:
            raise KeyError(f"no cell with {'cell_id' if isinstance(key, str) else 'index'} {key!r}") from None

    def records_of(self, key):
        cid = self.cell(key).cell_id
        return [r for r in self.records if r.cell_id == cid]

    def record(self, key, frequency_index):
        cid = self.cell(key).cell_id
        for r in self.records:
            if r.cell_id == cid and r.frequency_index == frequency_index:
                return r
        raise KeyError(f"no record for cell {cid!r} at frequency index {frequency_index}")

    @property
    def frequencies(self):
        """Sorted distinct record frequencies [Hz]."""
        return np.array(sorted({r.frequency for r in self.records}))

    def series(self, key, name):
        """(frequencies [Hz], complex values) of one component of one cell, by frequency index."""
        recs = self.records_of(key)
        return np.array([r.frequency for r in recs]), np.array([r.component(name) for r in recs])

    def phase_series_deg(self, key, name, unwrap=False, warn_step_deg=STEP_WARN_DEG):
        """Phase [deg] vs frequency for one cell/component: principal, or unwrapped along frequency.

        Unwrapping needs an unbroken run of frequency indices (a gap makes the step between
        neighbours meaningless) and a fine enough grid; a UserWarning is issued if the largest
        wrapped step exceeds ``warn_step_deg`` (phase.py explains why this cannot be detected
        reliably after the fact)."""
        recs = self.records_of(key)
        freqs = np.array([r.frequency for r in recs])
        phase = np.array([r.phase_deg(name) for r in recs])
        if not unwrap:
            return freqs, phase
        idx = [r.frequency_index for r in recs]
        if idx != list(range(idx[0], idx[0] + len(idx))):
            raise DatabaseError(f"cannot unwrap: frequency indices {idx} are not consecutive")
        step = max_wrapped_step_deg(phase)
        if step > warn_step_deg:
            import warnings

            warnings.warn(
                f"largest phase step between frequency samples is {step:.1f} deg "
                f"(> {warn_step_deg:g}); the unwrapped phase may be aliased",
                UserWarning,
                stacklevel=2,
            )
        return freqs, unwrap_deg(phase)

    def s_array(self, basis="linear"):
        """S as an array (n_cells, n_frequencies, 4, 4); row i is ``cells[i]``; NaN where a cell
        has no record. ``basis='circular'`` applies the change of basis of polarization.py."""
        if basis not in ("linear", "circular"):
            raise ValueError("basis must be 'linear' or 'circular'")
        n_f = 1 + max((r.frequency_index for r in self.records), default=-1)
        row = {c.cell_id: i for i, c in enumerate(self.cells)}
        out = np.full((len(self.cells), n_f, 4, 4), np.nan + 0j)
        for r in self.records:
            out[row[r.cell_id], r.frequency_index] = r.s
        return out if basis == "linear" else matrix4_to_circular(out)

    def parameter_array(self, name):
        """Values of one geometry parameter over all cells (in cell order)."""
        if name not in PARAMETER_NAMES:
            raise ValueError(f"unknown parameter {name!r}; use one of {list(PARAMETER_NAMES)}")
        return np.array([c.parameters[name] for c in self.cells])

    def unit_cell(self, key):
        """Rebuild the UnitCell this cell was solved for (needs the JSON metadata)."""
        from .sweep import unit_cell_from_context

        context = self.metadata.get("context")
        if not context:
            raise DatabaseError("metadata has no 'context' (was this database imported from CSV?)")
        cell = self.cell(key)
        if cell.status != "ok":
            raise DatabaseError(f"cell {cell.cell_id} is invalid: {cell.reason}")
        return unit_cell_from_context(context, cell.parameters)

    def failed_validation(self):
        """Cells whose a priori checks or refinement study did not pass (see discretization.py)."""
        bad = []
        for c in self.cells:
            v = c.validation or {}
            ref = v.get("refinement")
            if c.status == "ok" and (
                not v.get("area_ok", True)
                or (ref is not None and not ref.get("passed", True))
                or any(not r.diagnostics["checks_passed"] for r in self.records_of(c.index))
            ):
                bad.append(c)
        return bad

    # ------------------------------------------------------------------------ JSON
    def to_dict(self):
        cells_json = []
        for c in self.cells:
            cells_json.append(
                {
                    "index": c.index,
                    "cell_id": c.cell_id,
                    "status": c.status,
                    "reason": c.reason,
                    "parameters": {k: c.parameters[k] for k in PARAMETER_NAMES},
                    "validation": c.validation,
                    "records": [_record_to_json(r) for r in self.records_of(c.index)],
                }
            )
        return {
            "format": FORMAT_NAME,
            "version": FORMAT_VERSION,
            "metadata": self.metadata,
            "environment": self.environment,
            "cells": cells_json,
        }

    def to_json(self, path):
        text = json.dumps(self.to_dict(), indent=1, allow_nan=False) + "\n"
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)

    @classmethod
    def from_dict(cls, doc):
        if not isinstance(doc, dict) or doc.get("format") != FORMAT_NAME:
            raise DatabaseError(f"not a {FORMAT_NAME} document")
        if doc.get("version") != FORMAT_VERSION:
            raise DatabaseError(f"unsupported database version {doc.get('version')!r} (expected {FORMAT_VERSION})")
        for key in ("metadata", "cells"):
            if key not in doc:
                raise DatabaseError(f"missing top-level key {key!r}")
        metadata = doc["metadata"]
        context = metadata.get("context") if isinstance(metadata, dict) else None
        cells, records = [], []
        for c in doc["cells"]:
            try:
                params = {k: float(c["parameters"][k]) for k in PARAMETER_NAMES}
                cid = c["cell_id"]
                if context is not None and make_cell_id(context, params) != cid:
                    raise DatabaseError(
                        f"cell_id {cid!r} does not match its parameters and the metadata context "
                        f"(file edited, or produced by different code/settings)"
                    )
                cells.append(
                    CellEntry(int(c["index"]), cid, params, c["status"], c.get("reason"), c.get("validation"))
                )
                for r in c["records"]:
                    records.append(_record_from_json(cid, r))
            except (KeyError, TypeError) as e:
                raise DatabaseError(f"malformed cell entry: missing/invalid {e}") from None
        return cls(metadata, cells, records, doc.get("environment"))

    @classmethod
    def from_json(cls, path):
        try:
            with open(path, "r", encoding="utf-8") as fh:
                doc = json.loads(fh.read(), parse_constant=_reject_constant)
        except json.JSONDecodeError as e:
            raise DatabaseError(f"{path}: invalid JSON: {e}") from None
        return cls.from_dict(doc)

    # ------------------------------------------------------------------------- CSV
    def to_csv(self, path):
        """Write the flat CSV table. Raises ``DatabaseError`` if any record is one-port
        (grounded, Phase 12A): the fixed ``CSV_COLUMNS`` layout (16 linear + 16 circular columns)
        assumes a two-port S-matrix, and a mixed or one-port-only file would have to leave most
        columns meaninglessly empty rather than "undefined" (Sec. "S12/S21/S22 undefined" of
        SParameters) -- a minimal, explicit, documented restriction (Phase 12A brief Sec. 4's own
        "if a module can't support one-port records cleanly, reject with a clear error").
        Use ``to_json`` (which supports both topologies) for a grounded database instead.
        """
        if any(r.port_topology == "one_port" for r in self.records):
            raise DatabaseError(
                "CSV export does not support one-port (grounded, Phase 12A) records: the CSV "
                "schema's 16 linear + 16 circular columns assume a two-port S-matrix. Use "
                "Database.to_json() instead."
            )
        with open(path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, lineterminator="\n")
            writer.writerow(CSV_COLUMNS)
            for r in self.records:
                cell = self._by_id[r.cell_id]
                row = [r.cell_id, cell.index, r.frequency_index, repr(r.frequency)]
                row += [repr(cell.parameters[k]) for k in PARAMETER_NAMES]
                for key, typ in DIAGNOSTIC_TYPES.items():
                    v = r.diagnostics[key]
                    row.append(int(v) if typ in (bool, int) else repr(v))
                for basis in ("linear", "circular"):
                    for name in component_names(basis):
                        z = r.component(name)
                        row += [repr(abs(z)), repr(float(principal_phase_deg(z)))]
                writer.writerow(row)

    @classmethod
    def from_csv(cls, path, metadata=None):
        try:
            with open(path, "r", encoding="utf-8", newline="") as fh:
                rows = list(csv.reader(fh))
        except UnicodeDecodeError as e:
            raise DatabaseError(f"{path}: not a UTF-8 text file: {e}") from None
        if not rows or tuple(rows[0]) != CSV_COLUMNS:
            raise DatabaseError(f"{path}: header does not match the expected {len(CSV_COLUMNS)} columns")
        cells, records = {}, []
        for lineno, row in enumerate(rows[1:], start=2):
            if len(row) != len(CSV_COLUMNS):
                raise DatabaseError(f"{path}:{lineno}: expected {len(CSV_COLUMNS)} fields, got {len(row)}")
            try:
                rec, entry = _record_from_csv_row(dict(zip(CSV_COLUMNS, row)))
            except (ValueError, KeyError) as e:
                raise DatabaseError(f"{path}:{lineno}: {e}") from None
            prev = cells.setdefault(entry.cell_id, entry)
            if prev != entry:
                raise DatabaseError(f"{path}:{lineno}: inconsistent cell data for {entry.cell_id}")
            records.append(rec)
        ordered = sorted(cells.values(), key=lambda c: c.index)
        return cls(metadata or {}, ordered, records)


# ----------------------------------------------------------------------- (de)serialisation

def _reject_constant(name):
    raise DatabaseError(f"non-finite JSON constant {name!r} is not allowed")


def default_environment():
    """Provenance that may legitimately differ between machines; excluded from equality."""
    return {"python": platform.python_version(), "numpy": np.__version__}


def _entry_json(z):
    z = complex(z)
    return {
        "re": z.real,
        "im": z.imag,
        "mag": abs(z),
        "phase_deg": float(principal_phase_deg(z)),
    }


def _record_to_json(r):
    families = ("S11",) if r.port_topology == "one_port" else FAMILIES
    linear, circular = {}, {}
    for fam in families:
        block = r.s if r.port_topology == "one_port" else family_block(r.s, fam)
        linear[fam] = {c: _entry_json(block[idx]) for c, idx in LINEAR_COMPONENTS.items()}
        cblock = block_to_circular(block)
        circular[fam] = {
            c: {k: v for k, v in _entry_json(cblock[idx]).items() if k in ("mag", "phase_deg")}
            for c, idx in CIRCULAR_COMPONENTS.items()
        }
    return {
        "frequency_index": r.frequency_index,
        "frequency_hz": r.frequency,
        # Written for every new record (Phase 12A); absent in files written before Phase 12A,
        # where it is read back as "two_port" (the only topology that existed then) -- see
        # _record_from_json. Not a FORMAT_VERSION bump: this key is purely additive and optional,
        # so every pre-Phase-12A file still loads unchanged (bumping FORMAT_VERSION, which is
        # compared with strict equality, would instead break loading those old two_port-only
        # files, the opposite of what a version bump is for).
        "port_topology": r.port_topology,
        "diagnostics": dict(r.diagnostics),
        "S": linear,
        "S_circular": circular,
    }


_DERIVED_TOL = 1e-9
# Below this magnitude a phase is not reproducible from mag/phase columns (a difference of two
# nearly equal entries, e.g. S_LR of a symmetric cell, is rounding noise), so it is not checked.
_PHASE_MIN_MAG = 1e-6


def _check_derived(name, z, mag, phase_deg):
    if abs(abs(z) - mag) > _DERIVED_TOL * max(1.0, abs(z)):
        raise DatabaseError(f"{name}: stored magnitude {mag} disagrees with the value ({abs(z)})")
    if abs(z) > _PHASE_MIN_MAG and abs(wrapped_difference_deg(principal_phase_deg(z), phase_deg)) > 1e-6:
        raise DatabaseError(f"{name}: stored phase {phase_deg} deg disagrees with the value")


def _record_from_json(cid, r):
    port_topology = r.get("port_topology", "two_port")   # absent: written before Phase 12A
    families = ("S11",) if port_topology == "one_port" else FAMILIES
    s_mat = np.zeros((2, 2) if port_topology == "one_port" else (4, 4), dtype=complex)
    for fam in families:
        out, inc = (0, 0) if port_topology == "one_port" else (int(fam[1]) - 1, int(fam[2]) - 1)
        for comp, (row, col) in LINEAR_COMPONENTS.items():
            e = r["S"][fam][comp]
            z = complex(float(e["re"]), float(e["im"]))
            _check_derived(f"{fam}_{comp}", z, float(e["mag"]), float(e["phase_deg"]))
            s_mat[2 * out + row, 2 * inc + col] = z
    rec = Record(
        cid, int(r["frequency_index"]), float(r["frequency_hz"]), s_mat, r["diagnostics"],
        port_topology=port_topology,
    )
    for fam in families:
        for comp in CIRCULAR_COMPONENTS:
            e = r["S_circular"][fam][comp]
            _check_derived(f"{fam}_{comp}", rec.component(f"{fam}_{comp}"), float(e["mag"]), float(e["phase_deg"]))
    return rec


def _record_from_csv_row(row):
    s4 = np.zeros((4, 4), dtype=complex)
    for fam in FAMILIES:
        out, inc = int(fam[1]) - 1, int(fam[2]) - 1
        for comp, (r, c) in LINEAR_COMPONENTS.items():
            mag = float(row[f"{fam}_{comp}_mag"])
            phase = np.deg2rad(float(row[f"{fam}_{comp}_phase_deg"]))
            if mag < 0.0:
                raise ValueError(f"negative magnitude in {fam}_{comp}")
            s4[2 * out + r, 2 * inc + c] = mag * complex(np.cos(phase), np.sin(phase))
    diag = {}
    for key, typ in DIAGNOSTIC_TYPES.items():
        raw = row[key]
        diag[key] = bool(int(raw)) if typ is bool else (int(raw) if typ is int else float(raw))
    rec = Record(row["cell_id"], int(row["frequency_index"]), float(row["frequency_hz"]), s4, diag)
    for name in component_names("circular"):
        _check_derived(name, rec.component(name), float(row[f"{name}_mag"]), float(row[f"{name}_phase_deg"]))
    params = {k: float(row[k]) for k in PARAMETER_NAMES}
    entry = CellEntry(int(row["cell_index"]), row["cell_id"], params, "ok")
    return rec, entry
