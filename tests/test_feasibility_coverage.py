"""Phase 8A Sec. 8-16: library-only phase coverage analysis, read straight from
`unit_cell_library.json` and built on `phase_extraction`'s circular arithmetic (Phase 7,
unmodified)."""

import json
import re
import warnings

import pytest

from src.feasibility_coverage import (
    describe_coverage_vs_resonance,
    library_table_at_frequency,
    phase_coverage_vs_frequency,
    phase_density,
    resonance_summary,
    usable_geometry_count,
)
from src.feasibility_schema import DEFAULT_LIBRARY_DIR
from src.library_builder import build_library_database
from src.library_schema import LibraryConfig
from src.unit_cell_library import write_library_json

warnings.filterwarnings("ignore", category=UserWarning)

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness_m": 1.575e-3}


def small_config():
    return LibraryConfig(
        period_m=4.8e-3, substrate=_SUBSTRATE,
        f_min_hz=11.0e9, f_max_hz=13.0e9, f_step_hz=1.0e9, f_design_hz=12.0e9,
        l_values_mm=(1.0, 2.4, 4.2), w_values_mm=(1.0, 2.4, 4.2),
    )


@pytest.fixture(scope="module")
def small_library_doc(tmp_path_factory):
    cfg = small_config()
    database = build_library_database(cfg, out_dir=None).database
    path = write_library_json(database, tmp_path_factory.mktemp("phase8a_cov") / "lib.json")
    return json.loads(path.read_text()), cfg


def test_library_table_at_frequency_matches_number_of_ok_cells(small_library_doc):
    doc, cfg = small_library_doc
    table = library_table_at_frequency(doc, 12.0)
    assert len(table) == 9


def test_library_table_at_frequency_rejects_off_grid_frequency(small_library_doc):
    doc, _ = small_library_doc
    with pytest.raises(ValueError):
        library_table_at_frequency(doc, 12.3456)   # not one of 11/12/13 GHz


def test_phase_density_reports_consistent_geometry_count(small_library_doc):
    doc, _ = small_library_doc
    table = library_table_at_frequency(doc, 12.0)
    density = phase_density(table)
    assert density.n_valid_geometries == 9
    assert density.n_distinct_phase_states <= 9
    assert density.largest_gap_deg >= density.smallest_gap_deg
    assert density.phase_density_per_deg == pytest.approx(9 / 360.0)


def test_phase_coverage_vs_frequency_covers_every_stored_point(small_library_doc):
    doc, _ = small_library_doc
    result = phase_coverage_vs_frequency(doc)
    assert [r["frequency_GHz"] for r in result] == [11.0, 12.0, 13.0]
    for r in result:
        assert 0.0 <= r["coverage_deg"] <= 360.0


def test_resonance_summary_counts_are_consistent(small_library_doc):
    doc, _ = small_library_doc
    n_resonant, n_non_resonant, freqs = resonance_summary(doc)
    assert n_resonant + n_non_resonant == 9
    assert len(freqs) == n_resonant


def test_describe_coverage_vs_resonance_is_honest_when_no_resonance_found():
    freq_cov = [{"frequency_GHz": 11.0, "coverage_deg": 10.0}, {"frequency_GHz": 12.0, "coverage_deg": 20.0}]
    result = describe_coverage_vs_resonance(freq_cov, n_resonant_cells=0)
    assert result["resonant_cells_in_band"] == 0
    assert "no resonance-vs-coverage correlation can be computed" in result["note"]
    assert result["frequency_of_peak_coverage_GHz"] == 12.0


def test_describe_coverage_vs_resonance_does_not_claim_causality_when_resonance_found():
    freq_cov = [{"frequency_GHz": 11.0, "coverage_deg": 10.0}, {"frequency_GHz": 12.0, "coverage_deg": 40.0}]
    result = describe_coverage_vs_resonance(freq_cov, n_resonant_cells=2, resonant_freqs_hz=(11.5e9, 12.5e9))
    assert result["resonant_cells_in_band"] == 2
    assert "does not claim a causal relationship" in result["note"]


def test_usable_geometry_count_matches_phase_usability_threshold(small_library_doc):
    doc, _ = small_library_doc
    table = library_table_at_frequency(doc, 12.0)
    n_usable = usable_geometry_count(table)
    from src.library_classifier import phase_usability

    expected = sum(1 for row in table if phase_usability(row["magnitude_db"]) == "PHASE_USABLE")
    assert n_usable == expected


# --------------------------------------------------- Tests 4-5: regression against the real library

_PRODUCTION_JSON = DEFAULT_LIBRARY_DIR / "database" / "unit_cell_library.json"
_LIBRARY_REPORT = DEFAULT_LIBRARY_DIR / "reports" / "LIBRARY_REPORT.md"


@pytest.mark.skipif(not _PRODUCTION_JSON.exists(), reason="Phase 7 production library not built")
def test_recomputed_coverage_matches_the_stored_library_report_value():
    doc = json.loads(_PRODUCTION_JSON.read_text())
    table = library_table_at_frequency(doc, 12.5)
    from src.phase_extraction import phase_coverage

    cov = phase_coverage(table)

    report_text = _LIBRARY_REPORT.read_text()
    m = re.search(r"coverage = \*\*([0-9.]+) deg of 360 deg\*\*", report_text)
    assert m, "LIBRARY_REPORT.md coverage line not found"
    stored_coverage_deg = float(m.group(1))

    # Documented tolerance: LIBRARY_REPORT.md prints :.2f, so up to 0.005 deg of rounding.
    assert cov.coverage_deg == pytest.approx(stored_coverage_deg, abs=0.01)
    assert cov.coverage_deg == pytest.approx(24.84, abs=0.5)   # brief Sec. 11's expected ballpark


@pytest.mark.skipif(not _PRODUCTION_JSON.exists(), reason="Phase 7 production library not built")
def test_recomputed_largest_gap_matches_the_stored_library_report_value():
    doc = json.loads(_PRODUCTION_JSON.read_text())
    table = library_table_at_frequency(doc, 12.5)
    from src.phase_extraction import phase_coverage

    cov = phase_coverage(table)

    report_text = _LIBRARY_REPORT.read_text()
    m = re.search(r"\*\*([0-9.]+) deg\*\*, between", report_text)
    assert m, "LIBRARY_REPORT.md largest-gap line not found"
    stored_gap_deg = float(m.group(1))

    assert cov.largest_gap_deg == pytest.approx(stored_gap_deg, abs=0.01)
    assert cov.largest_gap_deg == pytest.approx(335.16, abs=0.5)   # brief Sec. 11's expected ballpark
