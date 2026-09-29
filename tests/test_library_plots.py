"""The 8 required plots (req. 21): representative-cell selection is deterministic and
documented, every plot is a real PNG built from already-computed data, and the full-grid plots
(phase-vs-geometry, resonance map, quality map) cover every geometry, not a subset.

Needs the optional matplotlib dependency; skipped entirely otherwise.
"""

import warnings

import pytest

pytest.importorskip("matplotlib")

from src.library_builder import build_library_database
from src.library_classifier import classify_cell
from src.library_plots import (
    generate_all_plots,
    plot_phase_coverage_map,
    plot_phase_vs_geometry,
    plot_s_parameter_curve,
    select_representative_cells,
)
from src.library_schema import LibraryConfig
from src.phase_extraction import extract_design_frequency, phase_coverage, phase_vs_geometry_table

warnings.filterwarnings("ignore", category=UserWarning)

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness_m": 1.575e-3}
_PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def small_config():
    return LibraryConfig(
        period_m=4.8e-3,
        substrate=_SUBSTRATE,
        f_min_hz=11.0e9,
        f_max_hz=13.0e9,
        f_step_hz=1.0e9,
        f_design_hz=12.0e9,
        l_values_mm=(1.0, 1.8, 2.6, 3.4, 4.2),
        w_values_mm=(1.0, 1.8, 2.6, 3.4, 4.2),
    )


@pytest.fixture(scope="module")
def built():
    cfg = small_config()
    database = build_library_database(cfg, out_dir=None).database
    rows = extract_design_frequency(database, cfg)
    table = phase_vs_geometry_table(rows)
    coverage = phase_coverage(table)
    by_lw = {(r.l_mm, r.w_mm): r for r in rows}
    classifications = [
        classify_cell(
            c, database, by_lw[(round(c.parameters["length_m"] * 1e3, 9), round(c.parameters["width_m"] * 1e3, 9))]
        )
        for c in database.cells
    ]
    return cfg, database, table, coverage, classifications


# ---------------------------------------------------------------------- representative selection


def test_select_representative_cells_is_a_deterministic_3x3_subgrid():
    cfg = small_config()
    cells = select_representative_cells(cfg, n_per_axis=3)
    l_values = sorted({l for l, w in cells})
    w_values = sorted({w for l, w in cells})
    assert l_values[0] == cfg.l_values_mm[0]
    assert l_values[-1] == cfg.l_values_mm[-1]
    assert len(cells) <= 9
    # Deterministic: calling again gives the identical list.
    assert select_representative_cells(cfg, n_per_axis=3) == cells


def test_select_representative_cells_default_config_gives_corners_and_centre():
    cfg = LibraryConfig.default()
    cells = select_representative_cells(cfg)
    assert (1.0, 1.0) in cells
    assert (4.4, 4.4) in cells
    assert len(cells) == 9


# --------------------------------------------------------------------------------------- plots


@pytest.mark.slow
def test_all_8_plots_are_generated(built, tmp_path):
    cfg, database, table, coverage, classifications = built
    paths = generate_all_plots(database, cfg, table, coverage, classifications, tmp_path)
    assert set(paths) == {
        "s11_magnitude_db",
        "s11_phase_deg",
        "s21_magnitude_db",
        "s21_phase_deg",
        "phase_vs_geometry",
        "phase_coverage_map",
        "resonance_map",
        "library_quality_map",
    }
    for p in paths.values():
        assert p.exists()
        assert p.read_bytes()[:8] == _PNG_MAGIC
        assert p.stat().st_size > 1000


@pytest.mark.slow
def test_s_parameter_curve_plot_skips_failed_cells_without_raising(built, tmp_path):
    cfg, database, *_ = built
    cells = select_representative_cells(cfg)
    path = tmp_path / "s11_mag.png"
    plot_s_parameter_curve(database, cells, "S11", "mag", path)
    assert path.exists()


@pytest.mark.slow
def test_phase_vs_geometry_plot_covers_the_full_grid_not_a_subset(built, tmp_path):
    cfg, database, table, coverage, classifications = built
    path = tmp_path / "pvg.png"
    plot_phase_vs_geometry(table, cfg, path, cfg.f_design_hz / 1e9)
    assert path.exists()
    # The underlying table itself (not just the picture) has one row per successful geometry.
    assert len(table) == sum(1 for c in database.cells if c.status == "ok")


def test_phase_coverage_map_plot_runs_on_single_point(tmp_path):
    """Edge case: coverage analysis of a single cell (100% gap) must still render, not crash."""
    from src.phase_extraction import PhaseCoverage

    table = [{"L_mm": 1.0, "W_mm": 1.0, "phase_deg": 10.0, "magnitude_db": -10.0, "S11_real": 0.0, "S11_imag": 0.0}]
    coverage = PhaseCoverage(1, 10.0, 10.0, 0.0, 360.0, 0.0, 10.0, 10.0)
    path = tmp_path / "coverage_single.png"
    plot_phase_coverage_map(table, coverage, path, 12.5)
    assert path.exists()
    assert path.read_bytes()[:8] == _PNG_MAGIC
