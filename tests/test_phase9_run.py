"""Phase 9 end-to-end pipeline and CLI, against the real Phase 7 checkpoint."""

import subprocess
import sys
from pathlib import Path

import pytest

from src.phase8c import data_source
from src.phase9.run import print_final_report, run_phase9

ROOT = Path(__file__).resolve().parent.parent
_HAS_CHECKPOINT = data_source.DEFAULT_CHECKPOINT_JSON.exists()
_LIBRARY_DIR = ROOT / "results" / "unit_cell_library"


def _hash(path):
    import hashlib
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_produces_all_required_output_files(tmp_path, phase8c_full_database):
    result = run_phase9(tmp_path / "out", database=phase8c_full_database)
    out = result["output_dir"]
    for name in ("phase9_report.md", "phase9_summary.json", "phase9_gate_status.json"):
        assert (out / name).exists(), name
    dir_9a, dir_9b = out / "phase9a_control_variable", out / "phase9b_local_periodicity"
    for name in ("hfss_reference_cells.csv", "s11_comparison.csv", "s21_comparison.csv", "architecture_audit.json"):
        assert (dir_9a / name).exists(), name
    for name in ("supercell_definitions.csv", "supercell_results.csv", "lpa_error_metrics.csv", "locality_analysis.json"):
        assert (dir_9b / name).exists(), name
    assert len(list((dir_9a / "plots").glob("*.png"))) >= 10
    assert len(list((dir_9b / "plots").glob("*.png"))) >= 8


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_honestly_reports_no_hfss_data(tmp_path, phase8c_full_database):
    result = run_phase9(tmp_path / "out", database=phase8c_full_database)
    assert result["hfss_stats_s11"] is None
    assert result["hfss_unavailable_reason"] is not None
    assert result["decision"].s11_status == "UNVALIDATED"
    assert result["decision"].s21_status == "UNVALIDATED"


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_finds_local_periodicity_invalidated_under_strong_contrast(tmp_path, phase8c_full_database):
    """The real, measured Phase 9 finding this whole phase was built to check."""
    result = run_phase9(tmp_path / "out", database=phase8c_full_database)
    strong = next(r for r in result["lpa_results"] if r.name == "3x3_strong")
    assert strong.relative_complex_error > 0.5
    assert result["decision"].local_periodicity_status == "INVALIDATED"


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_never_unlocks_geometry_inversion(tmp_path, phase8c_full_database):
    result = run_phase9(tmp_path / "out", database=phase8c_full_database)
    assert result["decision"].geometry_inversion_unlocked is False
    assert result["decision"].final_status in (
        "LOCAL_PERIODICITY_INVALIDATED", "NO_VALIDATED_FORMULATION",
    )


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_does_not_modify_phase7_checkpoint_or_delivered_library(tmp_path, phase8c_full_database):
    delivered = _LIBRARY_DIR / "database" / "unit_cell_library.json"
    checkpoint = data_source.DEFAULT_CHECKPOINT_JSON
    before = {"delivered": _hash(delivered), "checkpoint": _hash(checkpoint)}
    run_phase9(tmp_path / "out", database=phase8c_full_database)
    after = {"delivered": _hash(delivered), "checkpoint": _hash(checkpoint)}
    assert before == after


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_print_final_report_never_claims_geometry_unlocked(tmp_path, capsys, phase8c_full_database):
    result = run_phase9(tmp_path / "out", database=phase8c_full_database)
    print_final_report(result)
    out = capsys.readouterr().out
    assert "Geometry inversion:                   LOCKED" in out
    assert "Final MTS STEP generation:            LOCKED" in out
    assert "Solver / Phase 1-8 modified:          NO" in out


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_cli_runs_against_the_real_checkpoint(tmp_path):
    out_dir = tmp_path / "cli_out"
    proc = subprocess.run(
        [sys.executable, str(ROOT / "examples" / "run_phase9_validation.py"), "--output-dir", str(out_dir)],
        cwd=ROOT, capture_output=True, text=True, timeout=180,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "PHASE 9" in proc.stdout
    assert (out_dir / "phase9_report.md").exists()
