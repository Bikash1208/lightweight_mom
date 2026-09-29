"""Phase 7: run the locked L/W x frequency sweep through the existing, unmodified solver.

This module does not solve anything itself -- it builds a `sweep.SweepDefinition` (18 L values x
18 W values, `axes=(length_m, width_m)`, Cartesian product = 324 geometries; 201-point frequency
grid) from a `LibraryConfig` and calls `runner.run_sweep` (Phase 5, unmodified), which already
provides everything req. 4/27/28 ask for: deterministic cell IDs and ordering, an append-only
`checkpoint.jsonl` (so a multi-hour build safely resumes after an interruption instead of
restarting), and `on_invalid: skip` recording of any a priori invalid geometry with its reason
(none are expected here -- every L, W in the locked sweep is well inside the 4.8 mm period; see
`tests/test_library_builder.py`). A genuine solve exception (as opposed to an a priori invalid
geometry) is *not* caught anywhere in `runner.run_sweep` and is not caught here either: per req.
4's "otherwise fail the build explicitly", this module relies on the existing, unmodified
architecture's behaviour (abort with the original exception, checkpoint left intact for a
resumed retry) rather than adding a new try/except around the solver to paper over a failure.
"""

from dataclasses import dataclass
from pathlib import Path

from .library_schema import LibraryConfig
from .runner import run_sweep
from .sweep import ParameterAxis, SweepDefinition, ValidationPolicy


@dataclass(frozen=True)
class LibraryBuildResult:
    config: LibraryConfig
    database: object            # database.Database
    computed: int
    reused: int
    checkpoint_path: object
    json_path: object
    csv_path: object
    n_ok_cells: int
    n_invalid_cells: int


def build_sweep_definition(config: LibraryConfig):
    """The `SweepDefinition` for ``config``'s locked L/W x frequency grid. `on_invalid='skip'`
    (req. 4: record, don't abort, an a priori invalid geometry); validation left non-strict
    (`ValidationPolicy(strict=False)`, its default) so a cell with high mesh staircase error is
    still solved and its `metal_area_error`/`area_ok` recorded (Phase 5, unmodified), not
    rejected outright -- rejecting it would be a second, silent way to "not solve" a geometry
    that req. 4 asks this module to avoid."""
    from .sweep import FrequencyGrid

    l_values = tuple(v * 1e-3 for v in config.l_values_mm)
    w_values = tuple(v * 1e-3 for v in config.w_values_mm)
    frequency = FrequencyGrid(config.f_min_hz, config.f_max_hz, config.f_step_hz)
    axes = (ParameterAxis("length_m", l_values), ParameterAxis("width_m", w_values))
    return SweepDefinition(frequency=frequency, axes=axes, on_invalid="skip", validation=ValidationPolicy())


def build_library_database(config: LibraryConfig, out_dir=None, *, resume=True, progress=None):
    """Run ``config``'s sweep through `runner.run_sweep` (Phase 5, unmodified) and return the
    result. ``out_dir`` (recommended for a build this size): checkpoint.jsonl/database.json/csv
    on disk, safely resumable; ``None`` runs the whole sweep in memory (tests only)."""
    definition = build_sweep_definition(config)
    sim_config = config.build_config()
    result = run_sweep(sim_config, out_dir, definition=definition, resume=resume, progress=progress)

    n_ok_cells = sum(1 for c in result.database.cells if c.status == "ok")
    n_invalid_cells = sum(1 for c in result.database.cells if c.status != "ok")
    if len(result.database.cells) != config.n_geometries:
        raise RuntimeError(
            f"expected exactly {config.n_geometries} geometries (18x18), "
            f"got {len(result.database.cells)} -- library build is inconsistent with its own config"
        )

    return LibraryBuildResult(
        config, result.database, result.computed, result.reused,
        result.checkpoint_path, result.json_path, result.csv_path,
        n_ok_cells, n_invalid_cells,
    )
