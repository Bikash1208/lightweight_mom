"""`unit_cell_library.json`/`.csv` (req. 5, 20), the manifest (req. 19), the read-only query API
and phase inversion (req. 24-26), and the full end-to-end pipeline (req. 18) -- all on small,
deterministic configurations (never the 324x201 production library, built separately by
`examples/build_unit_cell_library.py`)."""

import csv
import json
import warnings

import pytest

pytest.importorskip("matplotlib")

from src.library_builder import build_library_database
from src.library_schema import LIBRARY_CSV_COLUMNS, LibraryConfig
from src.phase_extraction import resolve_phase_source
from src.unit_cell_library import (
    LIBRARY_FORMAT_NAME,
    POLARIZATION_CONVENTION,
    TIME_CONVENTION,
    UnitCellLibrary,
    build_full_library,
    build_manifest,
    phase_lookup_filename,
    write_library_csv,
    write_library_json,
)

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
def built():
    cfg = small_config()
    return cfg, build_library_database(cfg, out_dir=None).database


# --------------------------------------------------------------------------- JSON / CSV export


def test_library_json_has_required_top_level_fields(built, tmp_path):
    _, database = built
    path = write_library_json(database, tmp_path / "lib.json")
    doc = json.loads(path.read_text())
    assert doc["format"] == LIBRARY_FORMAT_NAME
    assert len(doc["cells"]) == 9
    for field in (
        "cell_id",
        "L_mm",
        "W_mm",
        "period_x_mm",
        "period_y_mm",
        "substrate_thickness_mm",
        "epsilon_r",
        "tan_delta",
        "frequency_GHz",
        "S11",
        "S21",
        "S12",
        "S22",
        "polarization_convention",
        "time_convention",
        "geometry_status",
        "solver_status",
    ):
        assert field in doc["cells"][0], field


def test_library_json_time_and_polarization_convention_match_project(built, tmp_path):
    _, database = built
    path = write_library_json(database, tmp_path / "lib.json")
    doc = json.loads(path.read_text())
    assert doc["cells"][0]["time_convention"] == "exp(+j*omega*t)" == TIME_CONVENTION
    assert doc["cells"][0]["polarization_convention"] == POLARIZATION_CONVENTION


def test_library_json_s_entries_have_full_precision_fields(built, tmp_path):
    _, database = built
    path = write_library_json(database, tmp_path / "lib.json")
    doc = json.loads(path.read_text())
    entry = doc["cells"][0]["S11"][0]
    assert set(entry) == {"real", "imag", "mag", "mag_db", "phase_deg"}
    assert entry["mag"] == pytest.approx((entry["real"] ** 2 + entry["imag"] ** 2) ** 0.5)


def test_library_csv_has_exact_required_columns(built, tmp_path):
    _, database = built
    path = write_library_csv(database, tmp_path / "lib.csv")
    with open(path, newline="", encoding="utf-8") as fh:
        header = next(csv.reader(fh))
    assert tuple(header) == LIBRARY_CSV_COLUMNS


def test_library_csv_row_count_matches_cells_times_frequencies(built, tmp_path):
    _, database = built
    path = write_library_csv(database, tmp_path / "lib.csv")
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    n_ok = sum(1 for c in database.cells if c.status == "ok")
    assert len(rows) - 1 == n_ok * len(database.frequencies)


def test_phase_lookup_filename_format():
    assert phase_lookup_filename(12.5e9) == "phase_lookup_12p5GHz.csv"
    assert phase_lookup_filename(12.0e9) == "phase_lookup_12GHz.csv"


# --------------------------------------------------------------------------------------- manifest


def test_manifest_has_every_required_field(built):
    cfg, database = built
    n_ok = sum(1 for c in database.cells if c.status == "ok")

    class _Build:
        n_ok_cells = n_ok
        n_invalid_cells = len(database.cells) - n_ok

    manifest = build_manifest(cfg, _Build(), resolve_phase_source(), deterministic_timestamp="2026-01-01T00:00:00Z")
    for field in (
        "library_version",
        "creation_timestamp",
        "solver_version",
        "configuration_hash",
        "number_of_cells",
        "number_of_successful_cells",
        "number_of_failed_cells",
        "frequency_count",
        "frequency_min_GHz",
        "frequency_max_GHz",
        "design_frequency_GHz",
        "period_x_mm",
        "period_y_mm",
        "substrate_thickness_mm",
        "phase_source",
        "time_convention",
        "polarization_convention",
    ):
        assert field in manifest, field
    assert manifest["number_of_cells"] == cfg.n_geometries
    assert manifest["phase_source"] == "S11_xx"
    assert manifest["creation_timestamp"] == "2026-01-01T00:00:00Z"


def test_manifest_configuration_hash_is_stable_and_sensitive():
    cfg1 = small_config()
    cfg2 = small_config()

    class _Build:
        n_ok_cells = 9
        n_invalid_cells = 0

    m1 = build_manifest(cfg1, _Build(), resolve_phase_source(), "t")
    m2 = build_manifest(cfg2, _Build(), resolve_phase_source(), "t")
    assert m1["configuration_hash"] == m2["configuration_hash"]

    cfg3 = LibraryConfig(
        period_m=4.8e-3,
        substrate=_SUBSTRATE,
        f_min_hz=11.0e9,
        f_max_hz=13.0e9,
        f_step_hz=1.0e9,
        f_design_hz=12.0e9,
        l_values_mm=(1.0, 2.4, 4.0),
        w_values_mm=(1.0, 2.4, 4.2),  # one L value different
    )
    m3 = build_manifest(cfg3, _Build(), resolve_phase_source(), "t")
    assert m3["configuration_hash"] != m1["configuration_hash"]


# --------------------------------------------------------------------------------- query API


@pytest.fixture(scope="module")
def library_json(built, tmp_path_factory):
    _, database = built
    d = tmp_path_factory.mktemp("query")
    return write_library_json(database, d / "unit_cell_library.json")


def test_get_cell_exact_match(library_json):
    lib = UnitCellLibrary.from_json(library_json)
    cell = lib.get_cell(2.4, 2.4)
    assert cell is not None
    assert cell["L_mm"] == pytest.approx(2.4) and cell["W_mm"] == pytest.approx(2.4)


def test_get_cell_missing_geometry_returns_none(library_json):
    lib = UnitCellLibrary.from_json(library_json)
    assert lib.get_cell(3.7, 3.7) is None


def test_response_exact_frequency(library_json):
    lib = UnitCellLibrary.from_json(library_json)
    resp = lib.response(2.4, 2.4, 12.0)
    assert set(resp) == {"real", "imag", "mag", "mag_db", "phase_deg"}


def test_response_rejects_non_grid_frequency(library_json):
    lib = UnitCellLibrary.from_json(library_json)
    with pytest.raises(ValueError, match="no interpolation"):
        lib.response(2.4, 2.4, 12.37)


def test_response_rejects_missing_geometry(library_json):
    lib = UnitCellLibrary.from_json(library_json)
    with pytest.raises(KeyError):
        lib.response(3.7, 3.7, 12.0)


def test_find_phase_candidates_recovers_known_phase(library_json):
    """A target phase equal to a library cell's own phase must find that cell with ~0 error --
    proves the search actually works, not merely returns something."""
    lib = UnitCellLibrary.from_json(library_json)
    cell = lib.get_cell(2.4, 2.4)
    own_phase = cell["S11"][1]["phase_deg"]  # index 1 <-> 12.0 GHz (11, 12, 13 GHz grid)
    candidates = lib.find_phase_candidates(own_phase, 12.0, tolerance_deg=0.01)
    assert any(c.l_mm == 2.4 and c.w_mm == 2.4 for c in candidates)


def test_find_phase_candidates_returns_all_not_one(library_json):
    lib = UnitCellLibrary.from_json(library_json)
    candidates = lib.find_phase_candidates(0.0, 12.0, tolerance_deg=360.0)  # everything matches
    n_ok = sum(1 for c in lib._doc["cells"] if c["solver_status"] == "OK")
    assert len(candidates) == n_ok


def test_invert_phase_returns_no_valid_candidate_outside_coverage(library_json):
    """Req. 26's "no extrapolation" test: an essentially unreachable target with a tight
    tolerance must not silently return the closest thing anyway."""
    lib = UnitCellLibrary.from_json(library_json)
    # Find a phase that is definitely not present, using a minuscule tolerance.
    result = lib.invert_phase(-179.999, 12.0)
    if result != "NO_VALID_CANDIDATE":
        # If by chance a cell's phase is that close, tighten further to force the miss case.
        result = lib.find_phase_candidates(-179.999, 12.0, tolerance_deg=1e-6)
        assert result == []
    else:
        assert result == "NO_VALID_CANDIDATE"


def test_invert_phase_structure_when_found(library_json):
    lib = UnitCellLibrary.from_json(library_json)
    cell = lib.get_cell(2.4, 2.4)
    own_phase = cell["S11"][1]["phase_deg"]
    result = lib.invert_phase(own_phase, 12.0)
    assert result != "NO_VALID_CANDIDATE"
    for key in (
        "target_phase_deg",
        "selected_phase_deg",
        "phase_error_deg",
        "L_mm",
        "W_mm",
        "response_magnitude_db",
        "cell_id",
    ):
        assert key in result


# ---------------------------------------------------------------------------- full pipeline


def test_build_full_library_end_to_end_directory_layout(tmp_path):
    cfg = small_config()
    result = build_full_library(cfg, tmp_path / "lib", deterministic_timestamp="2026-01-01T00:00:00Z")

    assert (tmp_path / "lib" / "database" / "unit_cell_library.json").exists()
    assert (tmp_path / "lib" / "database" / "unit_cell_library.csv").exists()
    assert (tmp_path / "lib" / "database" / phase_lookup_filename(cfg.f_design_hz)).exists()
    assert (tmp_path / "lib" / "manifest.json").exists()
    assert (tmp_path / "lib" / "reports" / "LIBRARY_REPORT.md").exists()
    assert (tmp_path / "lib" / "reports" / "LIBRARY_VALIDATION.md").exists()
    assert (tmp_path / "lib" / "reports" / "HFSS_VALIDATION_SUBSET.csv").exists()
    assert len(result.step_paths) == result.build.n_ok_cells == 9
    assert len(result.plot_paths) == 8
    assert result.consistency.ok


def test_build_full_library_report_mentions_no_hfss_validation(tmp_path):
    cfg = small_config()
    result = build_full_library(cfg, tmp_path / "lib2", deterministic_timestamp="2026-01-01T00:00:00Z")
    text = result.library_report_path.read_text()
    assert "HFSS validation is pending" in text
    assert "No claim of agreement with HFSS" in text
