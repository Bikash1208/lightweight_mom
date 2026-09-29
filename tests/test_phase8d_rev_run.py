"""Phase 8D-REV end-to-end pipeline and CLI, against the real Phase 7 checkpoint (UDB-01, UDB-11,
UDB-12)."""

import subprocess
import sys
from pathlib import Path

import pytest

from src.phase8c import data_source
from src.phase8d_rev.run import capture_phase7_hashes, print_final_report, run_phase8d_rev

ROOT = Path(__file__).resolve().parent.parent
_HAS_CHECKPOINT = data_source.DEFAULT_CHECKPOINT_JSON.exists()
_LIBRARY_DIR = ROOT / "results" / "unit_cell_library"


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_produces_all_required_output_files(tmp_path, phase8c_full_database):
    result = run_phase8d_rev(tmp_path / "out", database=phase8c_full_database, root_dir=tmp_path)
    out = result["output_dir"]
    for name in (
        "database_schema.md", "database_schema.json", "universal_database_spec.md", "query_spec.md",
        "library_build_request_schema.json", "validation_report.md", "benchmark_query_results.csv",
        "database_integrity.json",
    ):
        assert (out / name).exists(), name
    assert (tmp_path / "UNIVERSAL_UNIT_CELL_DATABASE.md").exists()


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_recovers_the_324_cell_benchmark_and_rejects_the_wrong_period(tmp_path, phase8c_full_database):
    result = run_phase8d_rev(tmp_path / "out", database=phase8c_full_database, root_dir=tmp_path)
    assert result["positive_result"].status == "MATCH"
    assert sum(e.n_ok_cells for e in result["positive_result"].matching_entries) == 324
    assert result["negative_result"].status == "DATABASE_COVERAGE_MISSING"


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_does_not_modify_phase7_checkpoint_or_delivered_library(tmp_path, phase8c_full_database):
    before = capture_phase7_hashes(data_source.DEFAULT_CHECKPOINT_JSON, _LIBRARY_DIR)
    run_phase8d_rev(tmp_path / "out", database=phase8c_full_database, root_dir=tmp_path)
    after = capture_phase7_hashes(data_source.DEFAULT_CHECKPOINT_JSON, _LIBRARY_DIR)
    assert before == after


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_leaves_no_scratch_files_in_the_output_directory(tmp_path, phase8c_full_database):
    result = run_phase8d_rev(tmp_path / "out", database=phase8c_full_database, root_dir=tmp_path)
    children = {p.name for p in result["output_dir"].iterdir()}
    assert not any(name.startswith("_integrity_tmp") for name in children)


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_reports_all_integrity_checks_passed(tmp_path, phase8c_full_database):
    result = run_phase8d_rev(tmp_path / "out", database=phase8c_full_database, root_dir=tmp_path)
    assert all(r.passed for _, r in result["integrity_results"])


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_print_final_report_never_claims_hfss_performed(tmp_path, capsys, phase8c_full_database):
    result = run_phase8d_rev(tmp_path / "out", database=phase8c_full_database, root_dir=tmp_path)
    print_final_report(result)
    out = capsys.readouterr().out
    assert "HFSS validation:          NOT PERFORMED" in out
    assert "Phase 7 modified:         NO" in out


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_regression_suite_criterion_is_false_unless_explicitly_confirmed(tmp_path, phase8c_full_database):
    result = run_phase8d_rev(tmp_path / "out", database=phase8c_full_database, root_dir=tmp_path)
    by_name = dict(result["criteria_results"])
    assert by_name["Full regression suite passes"] is False


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_cli_runs_against_the_real_checkpoint(tmp_path):
    out_dir = tmp_path / "cli_out"
    proc = subprocess.run(
        [sys.executable, str(ROOT / "examples" / "run_phase8d_rev_universal_database.py"),
         "--output-dir", str(out_dir)],
        cwd=ROOT, capture_output=True, text=True, timeout=120,
    )
    assert "PHASE 8D-REV" in proc.stdout
    assert "MATCH (324 geometries)" in proc.stdout
    assert (out_dir / "database_integrity.json").exists()
