"""Phase 13 Part 3: a small JSON-backed registry of unit-cell libraries the app knows about --
the read-only Phase 7 benchmark (auto-seeded if its checkpoint export is on disk) plus any library
the app itself has built (`routers/library.py`'s sweep endpoint). Deliberately lighter than
`phase8d_rev.universal_database.UniversalDatabase`: that class is an in-memory, per-script-run
registry with no persistence of its own, and its topology/family-aggregation machinery answers a
different question ("which topologies/periods exist across every registered library") than what
the app's reflectarray/HMSA pages need moment-to-moment ("give me phase candidates from library X
at frequency f", which `unit_cell_library.UnitCellLibrary` already answers directly from one
library's JSON). Wiring the app's libraries into a live `UniversalDatabase` too is a documented
follow-up (docs/ARCHITECTURE.md), not silently done here.
"""

import json
from pathlib import Path

from ..feasibility_schema import DEFAULT_LIBRARY_DIR
from .storage import output_root

BENCHMARK_LIBRARY_ID = "phase7_benchmark"
BENCHMARK_LIBRARY_PATH = DEFAULT_LIBRARY_DIR / "database" / "unit_cell_library.json"


def _registry_path():
    return output_root() / "library_registry.json"


def _load():
    path = _registry_path()
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _save(entries):
    output_root().mkdir(parents=True, exist_ok=True)
    _registry_path().write_text(json.dumps(entries, indent=1, sort_keys=True), encoding="utf-8", newline="\n")


def register(library_id, path, meta):
    """``meta``: free-form JSON-able summary (period_mm, backing, substrate, n_cells, ...) shown
    by the browse endpoint without re-reading the (potentially large) library file."""
    entries = _load()
    entries[library_id] = {"path": str(path), "meta": meta, "readonly": False}
    _save(entries)


def list_libraries():
    """The registered app-built libraries, plus the read-only Phase 7 benchmark if its export is
    present on disk (it is excluded from git -- CLAUDE.md -- so this is conditional, never
    assumed)."""
    out = dict(_load())
    if BENCHMARK_LIBRARY_PATH.exists() and BENCHMARK_LIBRARY_ID not in out:
        out[BENCHMARK_LIBRARY_ID] = {
            "path": str(BENCHMARK_LIBRARY_PATH),
            "meta": {
                "description": "Phase 7's own 324-cell, 4.8 mm ungrounded benchmark library.",
                "period_mm": 4.8,
                "backing": "air",
                "caveats": [
                    "one L/W pair (of 324) is mesh-duplicate-identical to another at this mesh "
                    "(Phase 12B's own duplicate-geometry audit) -- solves to bit-identical S-parameters "
                    "despite differing continuous geometry",
                    "ungrounded only",
                ],
            },
            "readonly": True,
        }
    return out


def resolve_path(library_id):
    entries = list_libraries()
    if library_id not in entries:
        raise KeyError(f"no registered library {library_id!r}")
    return Path(entries[library_id]["path"])
