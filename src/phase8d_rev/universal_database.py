"""Sec. 6, 17, 20: the registry itself. `UniversalDatabase` holds any number of `LibraryEntry`
objects (each one an independent solved `database.Database`) and answers "which families/periods/
materials/frequencies do we actually have data for" without ever assuming the answer is the
existing 4.8 mm Phase 7 build.
"""

from . import library_entry as _library_entry
from .family import family_from_entries


class UniversalDatabaseError(ValueError):
    pass


class UniversalDatabase:
    """A registry, not a solver. `register` never re-solves or copies the wrapped
    `database.Database`; it only records identity/provenance alongside a reference to it."""

    def __init__(self):
        self._entries = {}

    def register(self, entry: _library_entry.LibraryEntry):
        if entry.library_id in self._entries:
            raise UniversalDatabaseError(f"library_id {entry.library_id!r} already registered")
        self._entries[entry.library_id] = entry
        return entry

    def entries(self):
        return tuple(self._entries.values())

    def entry(self, library_id):
        try:
            return self._entries[library_id]
        except KeyError:
            raise KeyError(f"no registered library {library_id!r}") from None

    def topologies(self):
        return tuple(sorted({e.topology for e in self._entries.values()}))

    def entries_for_topology(self, topology):
        return tuple(e for e in self._entries.values() if e.topology == topology)

    def family(self, topology):
        """Sec. 9: the aggregate `UnitCellFamily` for ``topology``, derived from every registered
        entry that has it -- KeyError if none is registered (never a fabricated placeholder)."""
        entries = self.entries_for_topology(topology)
        if not entries:
            raise KeyError(f"no registered library entry has topology {topology!r}")
        return family_from_entries(topology, entries)

    def periods_for_topology(self, topology):
        """Sec. 10: which (p_x_mm, p_y_mm) pairs are actually validated for ``topology`` -- never
        assumes 4.8 mm or any other value."""
        return tuple(sorted({(e.p_x_mm, e.p_y_mm) for e in self.entries_for_topology(topology)}))


def with_benchmark(checkpoint_path=_library_entry.DEFAULT_CHECKPOINT_JSON,
                    library_dir=_library_entry.DEFAULT_LIBRARY_DIR):
    """A `UniversalDatabase` with exactly the existing, unmodified Phase 7 324-cell/4.8 mm library
    registered as its one `reference_validation_library` entry (Sec. 3, 20). This is the
    architecture's *benchmark* state, not its final state -- callers add further `LibraryEntry`
    objects via `register` as more libraries (other periods, topologies, substrates) are built."""
    udb = UniversalDatabase()
    udb.register(_library_entry.benchmark_library(checkpoint_path, library_dir))
    return udb
