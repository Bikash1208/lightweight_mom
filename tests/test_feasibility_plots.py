"""Phase 8A Sec. 43: the 6 required plots render as valid PNGs from small deterministic data."""

import warnings

import pytest

pytest.importorskip("matplotlib")

from src.feasibility_coverage import library_table_at_frequency, phase_coverage_vs_frequency
from src.feasibility_matching import match_required_phase
from src.feasibility_plots import (
    plot_circular_phase_coverage,
    plot_coverage_vs_frequency,
    plot_matching_error,
    plot_phase_vs_lw,
    plot_required_phase_histogram,
    plot_required_phase_map,
)
from src.library_builder import build_library_database
from src.library_schema import LibraryConfig
from src.phase_extraction import phase_coverage
from src.unit_cell_library import write_library_json

warnings.filterwarnings("ignore", category=UserWarning)

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness_m": 1.575e-3}


def small_config():
    return LibraryConfig(
        period_m=4.8e-3,
        substrate=_SUBSTRATE,
        f_min_hz=11.0e9,
        f_max_hz=13.0e9,
        f_step_hz=1.0e9,
        f_design_hz=12.0e9,
        l_values_mm=(1.0, 2.4, 4.2),
        w_values_mm=(1.0, 2.4, 4.2),
    )


@pytest.fixture(scope="module")
def small_library_doc(tmp_path_factory):
    import json

    cfg = small_config()
    database = build_library_database(cfg, out_dir=None).database
    path = write_library_json(database, tmp_path_factory.mktemp("phase8a_plots") / "lib.json")
    return json.loads(path.read_text()), cfg


def _is_png(path):
    return path.exists() and path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_plot_phase_vs_lw(small_library_doc, tmp_path):
    doc, cfg = small_library_doc
    table = library_table_at_frequency(doc, 12.0)
    path = tmp_path / "phase_vs_lw.png"
    plot_phase_vs_lw(table, cfg.l_values_mm, cfg.w_values_mm, 12.0, path)
    assert _is_png(path)


def test_plot_circular_phase_coverage(small_library_doc, tmp_path):
    doc, _ = small_library_doc
    table = library_table_at_frequency(doc, 12.0)
    cov = phase_coverage(table)
    path = tmp_path / "coverage.png"
    plot_circular_phase_coverage(table, cov, 12.0, path)
    assert _is_png(path)


def test_plot_coverage_vs_frequency(small_library_doc, tmp_path):
    doc, _ = small_library_doc
    freq_cov = phase_coverage_vs_frequency(doc)
    path = tmp_path / "coverage_vs_freq.png"
    plot_coverage_vs_frequency(freq_cov, path)
    assert _is_png(path)


def test_plot_required_phase_map(tmp_path):
    path = tmp_path / "required_phase_map.png"
    plot_required_phase_map([0, 1, 2], [0, 1, 2], [10.0, -90.0, 170.0], path)
    assert _is_png(path)


def test_plot_required_phase_histogram(tmp_path):
    path = tmp_path / "hist.png"
    plot_required_phase_histogram([10.0, -90.0, 170.0, 0.0], path)
    assert _is_png(path)


def test_plot_matching_error(small_library_doc, tmp_path):
    doc, _ = small_library_doc
    table = library_table_at_frequency(doc, 12.0)
    matches = match_required_phase([0.0, 90.0, 200.0], table, tolerance_deg=10.0)
    path = tmp_path / "matching_error.png"
    plot_matching_error(matches, path)
    assert _is_png(path)
