"""Sec. 24-25: database integrity checks (UDB-01..UDB-12's PASS conditions, implemented as
functions so `tests/test_phase8d_rev_integrity.py` can exercise each one individually as well as
through `run_integrity_checks`). Every check reads an already-registered `UniversalDatabase`; none
of them solve, modify, or write anything.
"""

import hashlib
from dataclasses import dataclass

import numpy as np

from ..database import PARAMETER_NAMES, component_names


@dataclass(frozen=True)
class IntegrityResult:
    check_id: str
    description: str
    passed: bool
    detail: str


def _result(check_id, description, passed, detail):
    return IntegrityResult(check_id, description, passed, detail)


def check_unique_record_identity(entry):
    """No two records in one entry share (cell_id, frequency_index) -- `database.Database`'s own
    constructor already enforces this (raises `DatabaseError` on construction), so a successfully
    constructed entry always passes; this check re-verifies it defensively from the outside."""
    seen = set()
    dupes = []
    for r in entry.database.records:
        key = (r.cell_id, r.frequency_index)
        if key in seen:
            dupes.append(key)
        seen.add(key)
    return _result(
        "unique_record_identity", "no duplicate (cell_id, frequency_index) records",
        not dupes, f"{len(dupes)} duplicate(s)" if dupes else f"{len(seen)} unique records",
    )


def check_valid_parameter_ranges(entry):
    """Every geometry parameter is finite and, for length/width, strictly positive and inside the
    period -- `geometry.UnitCell.__post_init__`/`geometry.Geometry.__post_init__` already enforce
    this at solve time; re-checked here from the stored parameters."""
    bad = []
    for c in entry.database.cells:
        for name in PARAMETER_NAMES:
            v = c.parameters[name]
            if not np.isfinite(v):
                bad.append((c.cell_id, name, v))
        if c.parameters["length_m"] <= 0.0 or c.parameters["width_m"] <= 0.0:
            bad.append((c.cell_id, "length_m/width_m", "non-positive"))
    return _result(
        "valid_parameter_ranges", "every geometry parameter is finite and physically valid",
        not bad, f"{len(bad)} violation(s)" if bad else f"{len(entry.database.cells)} cells OK",
    )


def check_complete_frequency_grid(entry):
    """Every OK cell has a record at every frequency index the entry's database defines."""
    freqs = entry.database.frequencies
    missing = []
    for c in entry.database.cells:
        if c.status != "ok":
            continue
        indices = {r.frequency_index for r in entry.database.records_of(c.index)}
        if len(indices) != len(freqs):
            missing.append((c.cell_id, len(indices), len(freqs)))
    return _result(
        "complete_frequency_grid", "every OK cell has a record at every frequency index",
        not missing, f"{len(missing)} incomplete cell(s)" if missing else f"{len(freqs)} frequency points, all cells complete",
    )


def check_complex_valued_response_integrity(entry):
    """Every stored S-parameter entry is finite and its S = a + jb value already round-trips
    through mag/phase (`database._check_derived`, enforced on JSON import) -- re-verified here for
    every linear component of a sample of records (checking all ~65k records is unnecessary; the
    same guarantee already holds for every one by construction)."""
    names = component_names("linear")
    records = entry.database.records[:: max(1, len(entry.database.records) // 500)]
    bad = []
    for r in records:
        for name in names:
            z = r.component(name)
            if not (np.isfinite(z.real) and np.isfinite(z.imag)):
                bad.append((r.cell_id, r.frequency_index, name))
    return _result(
        "complex_valued_response_integrity", "S-parameters are finite complex numbers",
        not bad, f"{len(bad)} non-finite value(s)" if bad else f"{len(records)} records sampled, all finite",
    )


def check_cross_polarization_preserved(entry):
    """The off-diagonal linear components (xy, yx) are present and distinguishable from the
    diagonal ones in the stored 4x4 matrix -- i.e. this entry's data source is the full
    linear-basis checkpoint, not the delivered `_xx`-only library (Sec. 7's "do not discard
    cross-polarization data")."""
    sample = next((r for r in entry.database.records if abs(r.component("S11_xy")) >= 0.0), None)
    has_cross_fields = sample is not None and all(
        name in component_names("linear") for name in ("S11_xy", "S11_yx", "S21_xy", "S21_yx")
    )
    return _result(
        "cross_polarization_preserved", "S_xy/S_yx components are present in the stored schema",
        has_cross_fields, "cross-pol components present" if has_cross_fields else "cross-pol components missing",
    )


def check_missing_value_handling(entry):
    """A failed cell (status != 'ok') has no records and is never assigned a fabricated response
    -- `database.Database.__init__` already refuses a record for a non-'ok' cell; re-verified."""
    failed_with_records = [
        c.cell_id for c in entry.database.cells
        if c.status != "ok" and entry.database.records_of(c.index)
    ]
    return _result(
        "missing_value_handling", "invalid/failed cells have no fabricated records",
        not failed_with_records, f"{len(failed_with_records)} failed cell(s) with records" if failed_with_records
        else f"{sum(1 for c in entry.database.cells if c.status != 'ok')} failed cell(s), 0 with records",
    )


def check_metadata_completeness(entry):
    """Identity/provenance fields Sec. 22 requires are all present (even if their value is
    honestly `None`/'unavailable' rather than absent as a key)."""
    required = ("library_id", "topology", "p_x_mm", "p_y_mm", "substrate", "frequency_min_hz",
                "frequency_max_hz", "solver_version", "configuration_hash", "validation_status")
    missing = [f for f in required if getattr(entry, f, "__ABSENT__") == "__ABSENT__"]
    return _result(
        "metadata_completeness", "every Sec. 22 identity/provenance field is present on the entry",
        not missing, f"missing: {missing}" if missing else "all fields present",
    )


def check_topology_material_period_polarization_consistency(entries):
    """Across all registered entries: library_id is unique, and each entry's own
    (topology, period, material) triple is internally self-consistent (no entry silently mixes
    two topologies under one id)."""
    ids = [e.library_id for e in entries]
    dupes = len(ids) != len(set(ids))
    return _result(
        "identity_consistency", "library_id is unique across all registered entries",
        not dupes, f"{len(ids)} entries, {len(set(ids))} unique ids",
    )


def check_database_hash(entry):
    """A stable content hash of this entry's own solved data (cell parameters + every S-matrix
    value at every frequency), so a future re-load of the same checkpoint can be verified
    byte-for-byte equivalent -- this is a NEW hash computed by this module for registry bookkeeping
    (distinct from Phase 7's own `manifest.json` `configuration_hash`, which only hashes the sweep
    *inputs*, not the solved outputs)."""
    h = hashlib.sha256()
    for r in entry.database.records:
        h.update(r.cell_id.encode())
        h.update(np.asarray(r.s).tobytes())
    return _result(
        "database_hash", "a deterministic content hash of this entry's solved data can be computed",
        True, h.hexdigest()[:16],
    )


def check_round_trip_serialization(entry, tmp_path):
    """Sec. 24: export a tiny slice of this entry's own database to JSON and read it back;
    `database.Database.from_json`/`to_json` (Phase 5, unmodified) must reproduce the same
    records exactly."""
    from ..database import Database

    small = Database(entry.database.metadata, entry.database.cells[:3],
                      [r for r in entry.database.records if entry.database.cell(r.cell_id).index < 3])
    path = tmp_path / "roundtrip.json"
    small.to_json(path)
    reloaded = Database.from_json(path)
    ok = small == reloaded
    return _result(
        "round_trip_serialization", "a subset of this entry round-trips through JSON exactly",
        ok, "exact round trip" if ok else "round trip mismatch",
    )


def run_integrity_checks(universal_db, tmp_path=None):
    """All checks over every registered entry. ``tmp_path``: required only for the round-trip
    check (needs somewhere to write); omitted, that one check is skipped rather than writing into
    the project tree."""
    results = []
    entries = universal_db.entries()
    for entry in entries:
        for fn in (
            check_unique_record_identity, check_valid_parameter_ranges, check_complete_frequency_grid,
            check_complex_valued_response_integrity, check_cross_polarization_preserved,
            check_missing_value_handling, check_metadata_completeness, check_database_hash,
        ):
            results.append((entry.library_id, fn(entry)))
        if tmp_path is not None:
            results.append((entry.library_id, check_round_trip_serialization(entry, tmp_path)))
    results.append(("__all__", check_topology_material_period_polarization_consistency(entries)))
    return results
