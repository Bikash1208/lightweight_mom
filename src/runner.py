"""Run a planned sweep: solve every (cell, frequency), checkpoint each result, assemble a Database.

The solver is used unchanged: each result is ``sparameters.compute_s_parameters`` for one unit
cell at one frequency. Nothing is shared between results (no warm start, no cache), so the
value of a record depends only on its own inputs.

Resume / restart
    With ``out_dir`` every finished result is appended to ``checkpoint.jsonl`` and flushed to
    disk immediately (one JSON object per line). A checkpoint entry is identified by a *key*,
    the hash of everything that determines its value:

        record key      = hash(cell_id, frequency, a-priori thresholds)
        refinement key  = hash(cell_id, refinement frequencies, refinement factor)

    where ``cell_id`` already contains the lattice, substrate, metal, geometry parameters, mesh,
    reference planes and code version (sweep.py). On restart, every planned result whose key is
    in the checkpoint is reused verbatim and only the missing ones are computed, so

    * an interrupted run continues where it stopped (also after a hard kill: a torn last line is
      detected, dropped and recomputed; the file is rewritten atomically from the valid lines);
    * a sweep may be extended (more frequencies or parameter values) and only the new results
      are computed;
    * changing anything that changes a result (mesh, material, geometry, code version) changes
      the key, so stale results are never reused; they stay in the file, unused.

    ``resume=False`` discards an existing checkpoint. A file that is not a checkpoint of this
    format is never overwritten (``CheckpointError``).

The assembled database (and the exported database.json / database.csv) does not depend on the
order in which results were computed, on whether the run was interrupted, or on timing: the
same inputs give byte-identical files (on the same platform and BLAS thread count).
"""

import json
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._canonical import short_hash
from .database import CONVENTIONS, CellEntry, Database, Record
from .discretization import check_frequency, refinement_study
from .green_function import LayeredGreenFunction
from .sparameters import compute_s_parameters, surface_wave_proximity
from .sweep import plan_sweep

CHECKPOINT_NAME = "checkpoint.jsonl"
DATABASE_JSON_NAME = "database.json"
DATABASE_CSV_NAME = "database.csv"
CHECKPOINT_FORMAT = "lightweight_mom.checkpoint"
CHECKPOINT_VERSION = 1


class CheckpointError(RuntimeError):
    """The checkpoint path holds something that is not a usable checkpoint."""


@dataclass(frozen=True, eq=False)
class SweepResult:
    database: Database
    computed: int  # (cell, frequency) results and refinement studies solved in this call
    reused: int  # taken from the checkpoint
    dropped: int  # unreadable checkpoint lines discarded (and recomputed if needed)
    checkpoint_path: Path | None = None
    json_path: Path | None = None
    csv_path: Path | None = None


# ------------------------------------------------------------------------------ keys


def record_key(cell_id, frequency, policy):
    return short_hash(
        {
            "type": "record",
            "cell_id": cell_id,
            "frequency_hz": float(frequency),
            "min_cells_per_wavelength": policy.min_cells_per_wavelength,
            "max_area_error": policy.max_area_error,
        },
        length=32,
    )


def refinement_key(cell_id, frequencies, refinement):
    return short_hash(
        {
            "type": "refinement",
            "cell_id": cell_id,
            "frequencies_hz": [float(f) for f in frequencies],
            "factor": refinement.factor,
        },
        length=32,
    )


# ------------------------------------------------------------------------ checkpoint file


def _matrix_to_json(m):
    m = np.asarray(m, dtype=complex)
    return {"re": m.real.tolist(), "im": m.imag.tolist()}


def _matrix_from_json(d, expected_shape=(4, 4)):
    m = np.array(d["re"], dtype=float) + 1j * np.array(d["im"], dtype=float)
    if m.shape != expected_shape:
        raise ValueError(f"matrix is not {expected_shape}")
    return m


def _valid_entry(entry):
    """Structural validation of one parsed checkpoint line; raises on any defect."""
    kind, key, data = entry["type"], entry["key"], entry["data"]
    if not isinstance(key, str):
        raise ValueError("key is not a string")
    if kind == "record":
        # "port_topology" (Phase 12A): absent in a checkpoint written before Phase 12A, which
        # only ever wrote two-port (4x4) records.
        port_topology = data.get("port_topology", "two_port")
        shape = (2, 2) if port_topology == "one_port" else (4, 4)
        float(data["frequency_hz"])
        Record(
            data["cell_id"],
            0,
            data["frequency_hz"],
            _matrix_from_json(data["s"], shape),
            data["diagnostics"],
            port_topology=port_topology,
        )
    elif kind == "refinement":
        float(data["floquet_error"])
        float(data["mesh_error"])
        [float(f) for f in data["frequencies_hz"]]
    else:
        raise ValueError(f"unknown entry type {kind!r}")


class Checkpoint:
    """Append-only, fsync'd JSONL store keyed by result key. ``path=None``: in memory only."""

    def __init__(self, path, resume=True, sync=True):
        self.path = Path(path) if path is not None else None
        self.sync = sync
        self.entries = {}  # key -> (type, data)
        self.dropped = 0
        self._fh = None
        if self.path is None:
            return
        if self.path.exists() and self.path.stat().st_size > 0:
            self._read(load_entries=resume)  # raises CheckpointError for a foreign file
        self._write_all()  # header + valid entries (compaction), atomically
        self._fh = open(self.path, "a", encoding="utf-8", newline="\n")

    # -- file handling
    def _header(self):
        return json.dumps({"format": CHECKPOINT_FORMAT, "version": CHECKPOINT_VERSION}) + "\n"

    def _write_all(self):
        """Atomically (re)write header + the entries currently held."""
        tmp = self.path.with_name(self.path.name + ".tmp")
        with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(self._header())
            for key, (kind, data) in self.entries.items():
                fh.write(self._line(key, kind, data))
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, self.path)

    @staticmethod
    def _line(key, kind, data):
        return json.dumps({"type": kind, "key": key, "data": data}, allow_nan=False) + "\n"

    def _read(self, load_entries):
        with open(self.path, encoding="utf-8", errors="replace") as fh:
            lines = fh.read().split("\n")
        try:
            head = json.loads(lines[0])
            ok = head.get("format") == CHECKPOINT_FORMAT and head.get("version") == CHECKPOINT_VERSION
        except (json.JSONDecodeError, AttributeError):
            ok = False
        if not ok:
            raise CheckpointError(
                f"{self.path} exists but is not a version-{CHECKPOINT_VERSION} checkpoint; "
                f"refusing to overwrite it (move it away or choose another out_dir)"
            )
        if not load_entries:
            return
        for line in lines[1:]:
            if not line.strip():
                continue
            try:
                entry = json.loads(line)
                _valid_entry(entry)
            except (ValueError, KeyError, TypeError, AttributeError):
                self.dropped += 1
                continue
            self.entries[entry["key"]] = (entry["type"], entry["data"])

    # -- access
    def get(self, key, kind):
        hit = self.entries.get(key)
        return hit[1] if hit is not None and hit[0] == kind else None

    def put(self, key, kind, data):
        self.entries[key] = (kind, data)
        if self._fh is not None:
            self._fh.write(self._line(key, kind, data))
            self._fh.flush()
            if self.sync:
                os.fsync(self._fh.fileno())

    def close(self):
        if self._fh is not None:
            self._fh.close()
            self._fh = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


# ------------------------------------------------------------------------------- solving


def _diagnostics(s, n_unknowns, chk, area_ok):
    """Reciprocity/passivity from the (2, 2) S11 block for a grounded (one-port) result, or the
    full (4, 4) matrix otherwise -- the same two properties (S = S^T; I - S^H S is PSD), just
    generalized to however many ports actually exist (Phase 12A). Same 7 diagnostic keys either
    way: the on-disk ``DIAGNOSTIC_TYPES`` schema (database.py) is unchanged."""
    m = s.s11 if s.grounded else s.matrix()
    n = m.shape[0]
    scale = float(np.max(np.abs(m)))
    reciprocity = float(np.max(np.abs(m - m.T)) / scale) if scale > 0.0 else 0.0
    passivity = float(np.min(np.linalg.eigvalsh(np.eye(n) - m.conj().T @ m)))
    return {
        "n_unknowns": int(n_unknowns),
        "incoming_residual": float(s.incoming_residual),
        "reciprocity_error": reciprocity,
        "passivity_min_eig": passivity,
        "cells_per_wavelength": chk["cells_per_wavelength"],
        "specular_only": chk["specular_only"],
        "checks_passed": bool(area_ok and chk["wavelength_ok"] and chk["specular_only"]),
    }


def _solve_record(plan, cell, frequency, policy):
    config = plan.config
    s = compute_s_parameters(
        cell.unit_cell,
        float(frequency),
        config.mesh,
        config.reference_planes,
        warn_grating=False,
        warn_surface_wave=False,
    )
    chk = check_frequency(cell.unit_cell, config.mesh, frequency, policy.min_cells_per_wavelength)
    return {
        "cell_id": cell.cell_id,
        "frequency_hz": float(frequency),
        "port_topology": "one_port" if s.grounded else "two_port",
        "s": _matrix_to_json(s.s11 if s.grounded else s.matrix()),
        "diagnostics": _diagnostics(s, cell.checks["n_unknowns"], chk, cell.checks["area_ok"]),
    }


def _frequency_grid_metadata(plan):
    f = plan.frequencies
    grid = plan.definition.frequency
    return {
        "start_hz": float(f[0]),
        "stop_hz": float(f[-1]),
        "step_hz": float(grid.step) if grid is not None else float(f[0]),
        "count": int(f.size),
        "frequencies_hz": [float(x) for x in f],
    }


def _metadata(plan):
    d = plan.definition
    return {
        "conventions": CONVENTIONS,
        "context": plan.context,
        "frequency_grid": _frequency_grid_metadata(plan),
        "parameter_axes": [{"name": a.name, "values": list(a.values)} for a in d.axes],
        "on_invalid": d.on_invalid,
        "validation_policy": d.validation.as_dict(),
        # Phase 12B: which OK cell indices are mesh-identical to each other on this sweep's own
        "mesh_duplicate_groups": [list(g) for g in plan.mesh_duplicate_groups],
        "n_mesh_distinct": plan.n_mesh_distinct,
    }


def run_sweep(config, out_dir=None, *, definition=None, resume=True, progress=None, sync=True):
    """Solve the sweep of ``definition`` (default ``config.sweep``) and return a ``SweepResult``.

    ``out_dir``: directory for checkpoint.jsonl, database.json and database.csv (created if
    needed); ``None`` runs in memory without resume. ``progress(event)`` is called after each
    (cell, frequency) result with a dict {done, total, cell_index, frequency_index, reused};
    an exception raised by it aborts the run and leaves the checkpoint valid.
    """
    # warn_duplicates=False: Phase 12B's mesh-duplicate-geometry warning is suppressed here (the
    # group data is still on plan.mesh_duplicate_groups); same convention as compute_s_parameters's
    # warn_grating/warn_surface_wave -- a sweep records the condition instead of warning about it.
    plan = plan_sweep(config, definition, warn_duplicates=False)  # validates before the first solve
    policy = plan.definition.validation
    ref = policy.refinement
    out = Path(out_dir) if out_dir is not None else None
    if out is not None:
        out.mkdir(parents=True, exist_ok=True)

    computed = reused = 0
    records, cells = [], []
    total = plan.n_solves
    done = 0

    with Checkpoint(out / CHECKPOINT_NAME if out else None, resume=resume, sync=sync) as ckpt:
        for cell in plan.cells:
            if cell.status != "ok":
                cells.append(CellEntry(cell.index, cell.cell_id, cell.parameters, "invalid", cell.reason))
                continue
            matrices = []
            for fi, f in enumerate(plan.frequencies):
                key = record_key(cell.cell_id, f, policy)
                data = ckpt.get(key, "record")
                hit = data is not None
                if not hit:
                    data = _solve_record(plan, cell, f, policy)
                    ckpt.put(key, "record", data)
                    computed += 1
                else:
                    reused += 1
                port_topology = data.get("port_topology", "two_port")  # Phase 12A; see _valid_entry
                shape = (2, 2) if port_topology == "one_port" else (4, 4)
                matrices.append(_matrix_from_json(data["s"], shape))
                records.append(
                    Record(
                        cell.cell_id,
                        fi,
                        float(f),
                        matrices[-1],
                        data["diagnostics"],
                        port_topology=port_topology,
                    )
                )
                done += 1
                if progress is not None:
                    progress(
                        {"done": done, "total": total, "cell_index": cell.index, "frequency_index": fi, "reused": hit}
                    )

            validation = dict(cell.checks)
            validation["refinement"] = None
            # Phase 12A: per-cell (not per-record, to keep the strict per-record diagnostic
            # schema unchanged for old-file compatibility, database.py Sec. 6A) surface-wave-pole
            # proximity, at the sweep's own highest frequency (smallest wavelength, largest k0 --
            # the same "worst case" convention RefinementPolicy.frequencies="max" already uses).
            # A no-op (None) for an ungrounded cell.
            if getattr(cell.unit_cell, "backing", "air") == "pec":
                worst_f = float(plan.frequencies[-1])
                green = LayeredGreenFunction(cell.unit_cell.substrate, worst_f, backing="pec")
                validation["surface_wave"] = surface_wave_proximity(cell.unit_cell, green, config.mesh.floquet_order)
                validation["surface_wave"]["frequency_hz"] = worst_f
            else:
                validation["surface_wave"] = None
            if ref.enabled:
                idx = [len(plan.frequencies) - 1] if ref.frequencies == "max" else list(range(len(plan.frequencies)))
                freqs = [float(plan.frequencies[i]) for i in idx]
                rkey = refinement_key(cell.cell_id, freqs, ref)
                rdata = ckpt.get(rkey, "refinement")
                if rdata is None:
                    errors = [
                        refinement_study(
                            cell.unit_cell,
                            config.reference_planes,
                            config.mesh,
                            plan.frequencies[i],
                            matrices[i],
                            ref.factor,
                        )
                        for i in idx
                    ]
                    rdata = {
                        "frequencies_hz": freqs,
                        "floquet_error": max(e[0] for e in errors),
                        "mesh_error": max(e[1] for e in errors),
                    }
                    ckpt.put(rkey, "refinement", rdata)
                    computed += 1
                else:
                    reused += 1
                validation["refinement"] = {
                    "frequencies_hz": [float(f) for f in rdata["frequencies_hz"]],
                    "factor": ref.factor,
                    "floquet_error": float(rdata["floquet_error"]),
                    "mesh_error": float(rdata["mesh_error"]),
                    "tolerance": ref.tolerance,
                    "passed": bool(rdata["floquet_error"] <= ref.tolerance and rdata["mesh_error"] <= ref.tolerance),
                }
            cells.append(CellEntry(cell.index, cell.cell_id, cell.parameters, "ok", None, validation))
        dropped = ckpt.dropped

    database = Database(_metadata(plan), cells, records)
    # Phase 12A: a grounded (one-port) sweep has no CSV representation (Database.to_csv's fixed
    # two-port column layout; database.py). Write JSON only, rather than let to_csv's ValueError
    # abort the whole run -- the JSON file is unaffected and remains the complete record.
    grounded_sweep = any(r.port_topology == "one_port" for r in records)
    json_path = csv_path = None
    if out is not None:
        json_path = out / DATABASE_JSON_NAME
        writers = [(json_path, database.to_json)]
        if not grounded_sweep:
            csv_path = out / DATABASE_CSV_NAME
            writers.append((csv_path, database.to_csv))
        for path, writer in writers:
            tmp = path.with_name(path.name + ".tmp")
            writer(tmp)
            os.replace(tmp, path)
    return SweepResult(
        database,
        computed,
        reused,
        dropped,
        out / CHECKPOINT_NAME if out else None,
        json_path,
        csv_path,
    )
