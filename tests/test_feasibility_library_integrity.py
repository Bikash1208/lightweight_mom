"""Phase 8A Sec. 32-34: library integrity, cross-checked from the actual delivered files (never
the manifest alone). Uses a small on-disk library built with `unit_cell_library.build_full_library`
for the general checks, and (when present) the real Phase 7 production library for the exact
counts brief Sec. 44 Tests 1-3 ask for."""

import json
import warnings

import pytest

pytest.importorskip("cadquery")
pytest.importorskip("matplotlib")

from src.feasibility_library_integrity import check_library_integrity, load_library_bundle
from src.feasibility_schema import DEFAULT_LIBRARY_DIR
from src.library_schema import LibraryConfig
from src.unit_cell_library import build_full_library

warnings.filterwarnings("ignore", category=UserWarning)

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness_m": 1.575e-3}


def small_config():
    return LibraryConfig(
        period_m=4.8e-3, substrate=_SUBSTRATE,
        f_min_hz=11.0e9, f_max_hz=13.0e9, f_step_hz=1.0e9, f_design_hz=12.0e9,
        l_values_mm=(1.0, 3.0), w_values_mm=(1.0, 3.0),
    )


@pytest.fixture(scope="module")
def built_library_dir(tmp_path_factory):
    base = tmp_path_factory.mktemp("phase8a_integrity_lib")
    cfg = small_config()
    build_full_library(cfg, base, deterministic_timestamp="2026-01-01T00:00:00Z")
    return base, cfg


def test_load_library_bundle_reads_all_four_files(built_library_dir):
    base, _ = built_library_dir
    bundle = load_library_bundle(base)
    assert bundle.manifest["number_of_cells"] == 4
    assert len(bundle.library_doc["cells"]) == 4
    assert len(bundle.step_files) == 4


def test_check_library_integrity_passes_on_a_clean_build(built_library_dir):
    base, cfg = built_library_dir
    bundle = load_library_bundle(base)
    result = check_library_integrity(bundle, expected_n_geometries=4, expected_n_frequency_points=3)
    assert result.ok, result.errors
    assert result.n_cells == 4
    assert result.n_ok_cells == 4
    assert result.n_step_files == 4


def test_check_library_integrity_detects_wrong_expected_count(built_library_dir):
    base, _ = built_library_dir
    bundle = load_library_bundle(base)
    result = check_library_integrity(bundle, expected_n_geometries=999, expected_n_frequency_points=3)
    assert not result.ok
    assert any("999" in e for e in result.errors)


# ------------------------------------------------------------------- Test 13: library corruption

def test_corrupted_json_record_is_detected_as_library_invalid(built_library_dir, tmp_path):
    base, cfg = built_library_dir
    src_json = base / "database" / "unit_cell_library.json"
    doc = json.loads(src_json.read_text())
    # Remove one frequency point (and its S11 entry) from the first OK cell: a genuine data
    # corruption, not a config mismatch.
    for cell in doc["cells"]:
        if cell["solver_status"] == "OK":
            cell["frequency_GHz"].pop(0)
            break

    corrupt_dir = tmp_path / "corrupt"
    (corrupt_dir / "database").mkdir(parents=True)
    (corrupt_dir / "geometry").mkdir()
    (corrupt_dir / "database" / "unit_cell_library.json").write_text(json.dumps(doc))
    (corrupt_dir / "database" / "unit_cell_library.csv").write_text(
        (base / "database" / "unit_cell_library.csv").read_text()
    )
    (corrupt_dir / "manifest.json").write_text((base / "manifest.json").read_text())
    for step in (base / "geometry").glob("*.step"):
        (corrupt_dir / "geometry" / step.name).write_text(step.read_text())

    bundle = load_library_bundle(corrupt_dir)
    result = check_library_integrity(bundle, expected_n_geometries=4, expected_n_frequency_points=3)
    assert not result.ok
    assert any("frequency points" in e for e in result.errors)


# --------------------------------------------------------- Tests 1-3: the real production library

_PRODUCTION_JSON = DEFAULT_LIBRARY_DIR / "database" / "unit_cell_library.json"


@pytest.mark.skipif(not _PRODUCTION_JSON.exists(), reason="Phase 7 production library not built")
def test_production_library_has_324_geometries():
    bundle = load_library_bundle(DEFAULT_LIBRARY_DIR)
    assert len(bundle.library_doc["cells"]) == 324


@pytest.mark.skipif(not _PRODUCTION_JSON.exists(), reason="Phase 7 production library not built")
def test_production_library_has_201_frequency_points_per_cell():
    bundle = load_library_bundle(DEFAULT_LIBRARY_DIR)
    for cell in bundle.library_doc["cells"]:
        if cell["solver_status"] == "OK":
            assert len(cell["frequency_GHz"]) == 201


@pytest.mark.skipif(not _PRODUCTION_JSON.exists(), reason="Phase 7 production library not built")
def test_production_library_has_125_ghz_as_an_exact_grid_point():
    bundle = load_library_bundle(DEFAULT_LIBRARY_DIR)
    cell = bundle.library_doc["cells"][0]
    assert any(abs(f - 12.5) < 1e-9 for f in cell["frequency_GHz"])


@pytest.mark.skipif(not _PRODUCTION_JSON.exists(), reason="Phase 7 production library not built")
def test_production_library_passes_integrity_check():
    bundle = load_library_bundle(DEFAULT_LIBRARY_DIR)
    result = check_library_integrity(bundle)
    assert result.ok, result.errors
