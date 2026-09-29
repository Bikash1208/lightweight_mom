"""HFSS validation subset selection (req. 22) and JSON/CSV internal consistency (req. 20) --
deterministic, and honest about what has (nothing) and has not been HFSS-validated."""

import csv
import warnings

import pytest

from src.library_builder import build_library_database
from src.library_classifier import classify_cell
from src.library_schema import LibraryConfig
from src.library_validation import (
    check_json_csv_consistency,
    render_library_validation_md,
    select_hfss_validation_subset,
    write_hfss_validation_subset_csv,
)
from src.phase_extraction import extract_design_frequency, phase_vs_geometry_table
from src.unit_cell_library import write_library_csv, write_library_json

pytestmark = pytest.mark.slow

warnings.filterwarnings("ignore", category=UserWarning)

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness_m": 1.575e-3}


@pytest.fixture(scope="module")
def built():
    cfg = LibraryConfig(
        period_m=4.8e-3,
        substrate=_SUBSTRATE,
        f_min_hz=11.0e9,
        f_max_hz=13.0e9,
        f_step_hz=1.0e9,
        f_design_hz=12.0e9,
        l_values_mm=(1.0, 1.8, 2.6, 3.4, 4.2),
        w_values_mm=(1.0, 1.8, 2.6, 3.4, 4.2),
    )
    database = build_library_database(cfg, out_dir=None).database
    rows = extract_design_frequency(database, cfg)
    table = phase_vs_geometry_table(rows)
    by_lw = {(r.l_mm, r.w_mm): r for r in rows}
    classifications = [
        classify_cell(
            c, database, by_lw[(round(c.parameters["length_m"] * 1e3, 9), round(c.parameters["width_m"] * 1e3, 9))]
        )
        for c in database.cells
    ]
    return cfg, database, table, classifications


# --------------------------------------------------------------------------- subset selection


def test_subset_is_deterministic_across_repeated_calls(built):
    cfg, _, table, classifications = built
    a = select_hfss_validation_subset(cfg, table, classifications)
    b = select_hfss_validation_subset(cfg, table, classifications)
    assert [(e.cell_id, e.reasons) for e in a] == [(e.cell_id, e.reasons) for e in b]


def test_subset_includes_smallest_and_largest_patch(built):
    cfg, _, table, classifications = built
    subset = select_hfss_validation_subset(cfg, table, classifications)
    lw = {(e.l_mm, e.w_mm) for e in subset}
    assert (min(cfg.l_values_mm), min(cfg.w_values_mm)) in lw
    assert (max(cfg.l_values_mm), max(cfg.w_values_mm)) in lw


def test_subset_includes_extreme_aspect_ratios(built):
    cfg, _, table, classifications = built
    subset = select_hfss_validation_subset(cfg, table, classifications)
    lw = {(e.l_mm, e.w_mm) for e in subset}
    assert (min(cfg.l_values_mm), max(cfg.w_values_mm)) in lw
    assert (max(cfg.l_values_mm), min(cfg.w_values_mm)) in lw


def test_subset_is_much_smaller_than_the_full_library(built):
    cfg, _, table, classifications = built
    subset = select_hfss_validation_subset(cfg, table, classifications)
    assert 0 < len(subset) < cfg.n_geometries


def test_subset_every_entry_has_a_reason(built):
    cfg, _, table, classifications = built
    subset = select_hfss_validation_subset(cfg, table, classifications)
    for entry in subset:
        assert len(entry.reasons) >= 1
        assert all(isinstance(r, str) and r for r in entry.reasons)


def test_write_hfss_validation_subset_csv_has_required_columns(built, tmp_path):
    cfg, _, table, classifications = built
    subset = select_hfss_validation_subset(cfg, table, classifications)
    path = write_hfss_validation_subset_csv(subset, tmp_path / "subset.csv")
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert set(rows[0]) == {"cell_id", "L_mm", "W_mm", "reason_for_selection"}
    assert len(rows) == len(subset)


# --------------------------------------------------------------------------- consistency check


def test_consistency_passes_on_a_clean_export(built, tmp_path):
    _, database, _, _ = built
    json_path = write_library_json(database, tmp_path / "lib.json")
    csv_path = write_library_csv(database, tmp_path / "lib.csv")
    report = check_json_csv_consistency(json_path, csv_path)
    assert report.ok
    assert report.mismatches == ()
    assert report.n_json_records > 0
    assert report.n_csv_rows == report.n_json_records


def test_consistency_detects_a_corrupted_csv_value(built, tmp_path):
    """Mutation check: deliberately corrupt one CSV cell and confirm detection, so a passing
    result above is meaningful and not vacuous."""
    _, database, _, _ = built
    json_path = write_library_json(database, tmp_path / "lib.json")
    csv_path = write_library_csv(database, tmp_path / "lib.csv")

    with open(csv_path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    header = rows[0]
    mag_col = header.index("S11_mag")
    rows[1][mag_col] = repr(float(rows[1][mag_col]) + 1.0)  # corrupt one value
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh, lineterminator="\n").writerows(rows)

    report = check_json_csv_consistency(json_path, csv_path)
    assert not report.ok
    assert len(report.mismatches) >= 1
    assert "S11_mag" in report.mismatches[0]


# ------------------------------------------------------------------------------------ report


def test_render_library_validation_md_states_no_hfss_comparison(built):
    from types import SimpleNamespace

    from src.library_validation import ConsistencyReport

    cfg, database, table, classifications = built
    subset = select_hfss_validation_subset(cfg, table, classifications)
    n_ok = sum(1 for c in database.cells if c.status == "ok")
    build = SimpleNamespace(n_ok_cells=n_ok, n_invalid_cells=len(database.cells) - n_ok, computed=0, reused=0)
    consistency = ConsistencyReport(True, 10, 10, ())
    text = render_library_validation_md(cfg, subset, consistency, build)
    assert "No HFSS comparison has been performed" in text
    assert "HFSS_VALIDATION_SUBSET.csv" in text
    for entry in subset:
        assert entry.cell_id in text
