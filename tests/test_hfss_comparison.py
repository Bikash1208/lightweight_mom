"""Phase 6A: HFSS comparison export for a single already-exported validation case.

Checks that `src/hfss_comparison.py` reads only the already-written `mom_results.json` (never
re-solves), and that the CSV/plots/README content it derives from that file is exactly
`20*log10(|S|)` / the stored principal phase, with the exact requested column names, no
unwrapping, no normalization and no interpolation -- verified against an independent read of the
same `mom_results.json` via `Database.from_json` and plain `numpy`, not against this module's own
intermediate values.

Needs the optional `cadquery`/`matplotlib` dependencies (Phase 6's own case export uses
`cadquery`); skipped entirely otherwise.
"""

import csv
import math

import numpy as np
import pytest

pytest.importorskip("cadquery")
pytest.importorskip("matplotlib")

from src.database import Database
from src.hfss_comparison import (
    HFSS_CSV_COLUMNS,
    export_hfss_comparison,
    hfss_comparison_table,
    load_case_artifacts,
)
from src.validation_cases import CASES_BY_ID
from src.validation_export import _DB_FLOOR, export_case

CASE = CASES_BY_ID["case01_non_resonant_patch"]


@pytest.fixture(scope="module")
def case_dir(tmp_path_factory):
    """One already-exported Phase 6 case folder (export_case, unmodified), used as the read-only
    input every test in this file exercises -- the fixture solves the case exactly once."""
    out = tmp_path_factory.mktemp("case01")
    export_case(CASE, out)
    return out


@pytest.fixture(scope="module")
def result(case_dir):
    return export_hfss_comparison(case_dir)


@pytest.fixture(scope="module")
def reference_db(case_dir):
    """An independent read of the same mom_results.json this module reads, used as ground truth
    (never the module's own intermediate table)."""
    return Database.from_json(case_dir / "mom_results.json")


# --------------------------------------------------------------------------------- I/O contract

def test_load_case_artifacts_reads_only_existing_exported_files(case_dir):
    artifacts = load_case_artifacts(case_dir)
    assert artifacts.geometry["case_id"] == CASE.case_id
    assert artifacts.solver_config["mesh"]["cells_x"] == CASE.config.mesh.cells_x
    assert len(artifacts.database.records) == len(CASE.frequency.frequencies())


def test_load_case_artifacts_never_solves_missing_case(tmp_path):
    """No mom_results.json in this empty directory: must raise, not fall back to solving."""
    with pytest.raises(FileNotFoundError, match="mom_results.json"):
        load_case_artifacts(tmp_path)


def test_export_writes_exactly_the_requested_files(result):
    assert result.csv_path.name == "hfss_comparison_reference.csv"
    assert result.md_path.name == "HFSS_COMPARISON.md"
    assert result.csv_path.exists() and result.md_path.exists()
    expected_plot_names = {
        "hfss_s11_magnitude_db.png", "hfss_s21_magnitude_db.png",
        "hfss_s11_phase_deg.png", "hfss_s21_phase_deg.png", "hfss_cross_polarization_db.png",
    }
    assert {p.name for p in result.plot_paths.values()} == expected_plot_names
    for p in result.plot_paths.values():
        assert p.exists()
        assert p.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


# ------------------------------------------------------------------------------------- CSV content

def test_csv_header_matches_exact_requested_columns(result):
    with open(result.csv_path, newline="", encoding="utf-8") as fh:
        header = next(csv.reader(fh))
    assert tuple(header) == HFSS_CSV_COLUMNS == (
        "frequency_GHz",
        "S11xx_mag_dB", "S11yy_mag_dB", "S21xx_mag_dB", "S21yy_mag_dB",
        "S11xx_phase_deg", "S11yy_phase_deg", "S21xx_phase_deg", "S21yy_phase_deg",
        "S11xy_mag_dB", "S11yx_mag_dB", "S21xy_mag_dB", "S21yx_mag_dB",
    )


def test_csv_frequency_points_match_8_to_12_ghz_half_step(result):
    with open(result.csv_path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))[1:]
    freqs = [float(r[0]) for r in rows]
    assert freqs == pytest.approx([8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0])


def test_csv_values_match_independent_read_of_mom_results(result, reference_db):
    """Every CSV cell recomputed here from a fresh `Database.from_json` read (not this module's
    own table) must match, to full float precision -- proves the CSV is a faithful, unmutated
    transcription of the already-computed data, not a re-derivation."""
    with open(result.csv_path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    header, data_rows = rows[0], rows[1:]

    freqs_hz, s11xx = reference_db.series(0, "S11_xx")
    for i, row in enumerate(data_rows):
        values = dict(zip(header, (float(x) for x in row)))
        assert values["frequency_GHz"] == pytest.approx(freqs_hz[i] / 1e9)

        for comp, mag_col, phase_col in (
            ("S11_xx", "S11xx_mag_dB", "S11xx_phase_deg"),
            ("S11_yy", "S11yy_mag_dB", "S11yy_phase_deg"),
            ("S21_xx", "S21xx_mag_dB", "S21xx_phase_deg"),
            ("S21_yy", "S21yy_mag_dB", "S21yy_phase_deg"),
        ):
            _, s = reference_db.series(0, comp)
            _, phase_deg = reference_db.phase_series_deg(0, comp, unwrap=False)
            expected_mag_db = 20.0 * math.log10(max(abs(s[i]), 10 ** (_DB_FLOOR / 20.0)))
            assert values[mag_col] == pytest.approx(expected_mag_db, abs=1e-9)
            assert values[phase_col] == pytest.approx(phase_deg[i], abs=1e-9)

        for comp, mag_col in (
            ("S11_xy", "S11xy_mag_dB"), ("S11_yx", "S11yx_mag_dB"),
            ("S21_xy", "S21xy_mag_dB"), ("S21_yx", "S21yx_mag_dB"),
        ):
            _, s = reference_db.series(0, comp)
            expected_mag_db = 20.0 * math.log10(max(abs(s[i]), 10 ** (_DB_FLOOR / 20.0)))
            assert values[mag_col] == pytest.approx(expected_mag_db, abs=1e-9)


def test_unrotated_patch_cross_pol_is_at_the_db_floor(result):
    """case01 has rotation_deg = 0: every cross-pol column should sit at the documented floor
    (machine-precision zero magnitude), demonstrating -- not merely claiming -- negligible
    cross-polarization for this case."""
    with open(result.csv_path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    header, data_rows = rows[0], rows[1:]
    for row in data_rows:
        values = dict(zip(header, (float(x) for x in row)))
        for col in ("S11xy_mag_dB", "S11yx_mag_dB", "S21xy_mag_dB", "S21yx_mag_dB"):
            assert values[col] == pytest.approx(_DB_FLOOR, abs=1e-6)


# ---------------------------------------------------------------------------------- table helper

def test_hfss_comparison_table_matches_database_series_directly(reference_db):
    table, mags, phases = hfss_comparison_table(reference_db)
    freqs_hz, s = reference_db.series(0, "S11_xx")
    assert np.allclose(table["frequency_GHz"], freqs_hz / 1e9)
    assert np.allclose(mags["S11_xx"], np.abs(s))
    _, phase_deg = reference_db.phase_series_deg(0, "S11_xx", unwrap=False)
    assert np.allclose(phases["S11_xx"], phase_deg)
    # every phase value must be a principal value (-180, 180], never unwrapped beyond that range
    assert np.all(phases["S11_xx"] > -180.0 - 1e-9) and np.all(phases["S11_xx"] <= 180.0 + 1e-9)


# ------------------------------------------------------------------------------------ README content

def test_readme_states_every_required_fact(result):
    text = result.md_path.read_text()
    assert "exp(+j*omega*t)" in text
    assert "z_top" in text and "+5 mm" in text
    assert "z_bottom" in text and "-6.575 mm" in text
    assert "20*log10" in text
    assert "no unwrapping" in text.lower() or "no phase unwrapping" in text.lower()
    assert "no empirical correction" in text.lower()
    assert "no normalization" in text.lower()
    assert "conjugate" in text.lower()  # time-convention mismatch guidance
    assert "Floquet" in text
    for comp in ("S11xx", "S11yy", "S21xx", "S21yy", "S11xy", "S11yx", "S21xy", "S21yx"):
        assert comp in text


def test_readme_lists_every_plot_and_the_csv(result):
    text = result.md_path.read_text()
    assert result.csv_path.name in text
    for p in result.plot_paths.values():
        assert p.name in text


def test_readme_frequency_list_matches_the_csv(result):
    text = result.md_path.read_text()
    assert "8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12" in text


# ------------------------------------------------------------------------- determinism / no invention

def test_repeated_export_from_the_same_mom_results_is_numerically_identical(case_dir, tmp_path):
    """Re-running the Phase 6A export against the same, already-written mom_results.json (no
    re-solve) must reproduce identical CSV values -- this module performs no computation whose
    result could vary between calls."""
    r1 = export_hfss_comparison(case_dir)
    with open(r1.csv_path, newline="", encoding="utf-8") as fh:
        rows1 = list(csv.reader(fh))
    r2 = export_hfss_comparison(case_dir)  # re-export into the same folder
    with open(r2.csv_path, newline="", encoding="utf-8") as fh:
        rows2 = list(csv.reader(fh))
    assert rows1 == rows2
