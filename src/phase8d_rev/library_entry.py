"""Sec. 6, 17, 22: one registered `database.Database` plus the identity/provenance metadata a
`UniversalDatabase` needs to filter and report on it. A `LibraryEntry` never copies or mutates the
`database.Database` it wraps -- it is read-only, held by reference, exactly like `phase8c` already
treats the Phase 7 checkpoint.
"""

import json
from dataclasses import dataclass
from pathlib import Path

from ..phase8c.data_source import DEFAULT_CHECKPOINT_JSON, load_full_database
from ..feasibility_schema import DEFAULT_LIBRARY_DIR
from .validation_levels import library_validation_status
from .record_schema import record_from_database

BENCHMARK_LIBRARY_ID = "phase7_reference_validation_library_4p8mm"
BENCHMARK_TOPOLOGY = "rectangular_patch"

# Sec. 3: the existing Phase 7 324-cell/4.8 mm dataset's role in the universal architecture --
# never "the" database.
REFERENCE_VALIDATION_LIBRARY = "reference_validation_library"
PROJECT_LIBRARY = "project_library"
ROLES = (REFERENCE_VALIDATION_LIBRARY, PROJECT_LIBRARY)


@dataclass(frozen=True)
class LibraryEntry:
    library_id: str
    topology: str
    role: str
    database: object                 # database.Database, held read-only by reference
    p_x_mm: float
    p_y_mm: float
    substrate: dict                  # name, eps_r, tan_delta, thickness_mm
    frequency_min_hz: float
    frequency_max_hz: float
    n_frequency_points: int
    polarizations_supported: tuple
    incidence_supported: tuple
    solver_version: str | None
    configuration_hash: str | None
    source_type: str
    validation_status: str
    backing: str = "air"             # "air" | "pec" (Phase 12A stack/backing tag, Sec. 4)

    def __post_init__(self):
        if self.role not in ROLES:
            raise ValueError(f"role {self.role!r} not in {ROLES}")
        if self.backing not in ("air", "pec"):
            raise ValueError(f"backing must be 'air' or 'pec', got {self.backing!r}")

    @property
    def n_geometries(self):
        return len(self.database.cells)

    @property
    def n_ok_cells(self):
        return sum(1 for c in self.database.cells if c.status == "ok")

    def matches_period(self, p_x_mm, p_y_mm, rtol=1e-9):
        return (
            abs(self.p_x_mm - p_x_mm) <= rtol * max(abs(p_x_mm), 1.0)
            and abs(self.p_y_mm - p_y_mm) <= rtol * max(abs(p_y_mm), 1.0)
        )

    def has_frequency(self, frequency_hz, rtol=1e-9):
        freqs = self.database.frequencies
        return bool((abs(freqs - frequency_hz) <= rtol * max(abs(frequency_hz), 1.0)).any())

    def record(self, L_mm, W_mm, frequency_hz, extra=None):
        """One `record_schema.UnitCellRecord`, built from an exact (L, W, frequency) already in
        this entry's database -- KeyError if not present (Sec. 26: no interpolation)."""
        cell = next(
            (c for c in self.database.cells
             if abs(round(c.parameters["length_m"] * 1e3, 9) - L_mm) < 1e-6
             and abs(round(c.parameters["width_m"] * 1e3, 9) - W_mm) < 1e-6),
            None,
        )
        if cell is None or cell.status != "ok":
            raise KeyError(f"no OK cell at L={L_mm} mm, W={W_mm} mm in {self.library_id}")
        freqs = self.database.frequencies
        matches = [i for i, f in enumerate(freqs) if abs(f - frequency_hz) <= 1e-9 * max(abs(frequency_hz), 1.0)]
        if not matches:
            raise KeyError(f"frequency {frequency_hz} Hz is not exact grid point of {self.library_id}")
        db_record = self.database.record(cell.index, matches[0])
        merged_extra = {
            "solver_version": self.solver_version, "configuration_hash": self.configuration_hash,
            "source_type": self.source_type,
        }
        merged_extra.update(extra or {})
        return record_from_database(self.database, self.library_id, self.topology, cell, db_record, merged_extra)

    def all_ok_records(self, frequency_hz, extra=None):
        """One `UnitCellRecord` per OK cell at an exact stored frequency (Sec. 20's benchmark
        query: recovering all 324 states)."""
        freqs = self.database.frequencies
        matches = [i for i, f in enumerate(freqs) if abs(f - frequency_hz) <= 1e-9 * max(abs(frequency_hz), 1.0)]
        if not matches:
            raise KeyError(f"frequency {frequency_hz} Hz is not an exact grid point of {self.library_id}")
        idx = matches[0]
        merged_extra = {
            "solver_version": self.solver_version, "configuration_hash": self.configuration_hash,
            "source_type": self.source_type,
        }
        merged_extra.update(extra or {})
        out = []
        for cell in self.database.cells:
            if cell.status != "ok":
                continue
            db_record = self.database.record(cell.index, idx)
            out.append(record_from_database(self.database, self.library_id, self.topology, cell, db_record, merged_extra))
        return out


def from_database(database, library_id, topology, role, manifest=None):
    """Wrap an already-solved `database.Database` (Phase 5/7, unmodified) as a `LibraryEntry`.
    ``manifest``: the optional `manifest.json`-shaped dict (Phase 7's own build manifest) this
    entry's solver_version/configuration_hash come from -- left `None`/`None` when no manifest is
    available, never guessed."""
    ctx = database.metadata.get("context") or {}
    lattice = ctx.get("lattice", {})
    substrate_ctx = ctx.get("substrate", {})
    freqs = database.frequencies
    if freqs.size == 0:
        raise ValueError(f"{library_id}: database has no records")

    return LibraryEntry(
        library_id=library_id, topology=topology, role=role, database=database,
        p_x_mm=float(lattice["period_x_m"] * 1e3), p_y_mm=float(lattice["period_y_m"] * 1e3),
        substrate={
            "name": substrate_ctx.get("name", "unavailable"),
            "eps_r": substrate_ctx.get("eps_r", float("nan")),
            "tan_delta": substrate_ctx.get("tan_delta", float("nan")),
            "thickness_mm": float(substrate_ctx["thickness_m"] * 1e3) if "thickness_m" in substrate_ctx else float("nan"),
        },
        # Phase 12A: absent in a context written before Phase 12A -- that was always "air"
        # (sweep.py's own unit_cell_from_context default), never a mix.
        backing=ctx.get("backing", "air"),
        frequency_min_hz=float(freqs.min()), frequency_max_hz=float(freqs.max()),
        n_frequency_points=int(freqs.size),
        polarizations_supported=("x", "y"),          # database.CONVENTIONS: both linear pols, both ports
        incidence_supported=("normal",),               # phase8b.conventions.NORMAL_INCIDENCE_ONLY
        solver_version=(manifest or {}).get("solver_version"),
        configuration_hash=(manifest or {}).get("configuration_hash"),
        source_type="MoM",
        validation_status=library_validation_status(database.records),
    )


def benchmark_library(checkpoint_path=DEFAULT_CHECKPOINT_JSON, library_dir=DEFAULT_LIBRARY_DIR):
    """Sec. 3, 20: the existing, unmodified Phase 7 324-cell/4.8 mm build, read-only, as one
    `reference_validation_library` entry. Uses the full linear-basis checkpoint
    (`phase8c.data_source.load_full_database`, itself read-only and self-cross-checking) so cross-
    polarization data is preserved (Sec. 7's explicit requirement), not the `_xx`-only delivered
    library."""
    database = load_full_database(checkpoint_path)
    manifest_path = Path(library_dir) / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else None
    return from_database(database, BENCHMARK_LIBRARY_ID, BENCHMARK_TOPOLOGY, REFERENCE_VALIDATION_LIBRARY, manifest)
