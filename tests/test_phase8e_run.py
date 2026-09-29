"""Phase 8E end-to-end pipeline and CLI, against the real Phase 7 checkpoint."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from src.phase8c import data_source
from src.phase8e.run import run_phase8e, print_final_report

ROOT = Path(__file__).resolve().parent.parent
_HAS_CHECKPOINT = data_source.DEFAULT_CHECKPOINT_JSON.exists()
_LIBRARY_DIR = ROOT / "results" / "unit_cell_library"
_UNIVERSAL_DB_DIR = ROOT / "results" / "universal_database"


def _hash(path):
    import hashlib
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_produces_all_required_output_files(tmp_path, phase8c_full_database):
    result = run_phase8e(tmp_path / "out", database=phase8c_full_database)
    out = result["output_dir"]
    for name in (
        "report.md", "formulation_candidates.json", "formulation_comparison.csv",
        "control_variable_analysis.json", "reference_spectrum.csv", "object_spectrum.csv",
        "required_response_map.csv", "library_response_map.csv", "reconstruction_metrics.json",
        "required_library_specification.json",
    ):
        assert (out / name).exists(), name
    plot_files = list((out / "plots").glob("*.png"))
    assert len(plot_files) >= 17


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_never_forces_a_positive_result(tmp_path, phase8c_full_database):
    result = run_phase8e(tmp_path / "out", database=phase8c_full_database)
    assert result["overall_status"] == "NO_VALIDATED_FORMULATION"
    assert all(not d.geometry_inversion_permitted for d in result["decisions"])


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_finds_formulation_a_and_d_architecture_incompatible(tmp_path, phase8c_full_database):
    result = run_phase8e(tmp_path / "out", database=phase8c_full_database)
    assert result["architecture_results"]["A"].status == "ARCHITECTURE_INCOMPATIBLE"
    assert result["architecture_results"]["D"].status == "ARCHITECTURE_INCOMPATIBLE"


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_beam_reconstruction_identifies_all_three_demo_beams(tmp_path, phase8c_full_database):
    """A real correctness check on the far-field math (Sec. 18): the continuous phase-difference
    hologram must actually place a detectable peak near each of the 3 demo beam directions."""
    result = run_phase8e(tmp_path / "out", database=phase8c_full_database)
    metrics = result["metrics_cont_b"]
    assert all(e < 5.0 for e in metrics.pointing_error_deg)   # within the far-field grid's own 2 deg resolution


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_run_does_not_modify_phase7_or_universal_database_outputs(tmp_path, phase8c_full_database):
    delivered = _LIBRARY_DIR / "database" / "unit_cell_library.json"
    checkpoint = data_source.DEFAULT_CHECKPOINT_JSON
    before = {"delivered": _hash(delivered), "checkpoint": _hash(checkpoint)}
    if _UNIVERSAL_DB_DIR.exists():
        before["udb_schema"] = _hash(_UNIVERSAL_DB_DIR / "database_schema.json")

    run_phase8e(tmp_path / "out", database=phase8c_full_database)

    after = {"delivered": _hash(delivered), "checkpoint": _hash(checkpoint)}
    if _UNIVERSAL_DB_DIR.exists():
        after["udb_schema"] = _hash(_UNIVERSAL_DB_DIR / "database_schema.json")
    assert before == after


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_print_final_report_never_claims_hfss_performed(tmp_path, capsys, phase8c_full_database):
    result = run_phase8e(tmp_path / "out", database=phase8c_full_database)
    print_final_report(result)
    out = capsys.readouterr().out
    assert "HFSS validation:          NOT PERFORMED" in out
    assert "Solver/Phase 7 modified:  NO" in out
    assert "Phase 8D-REV modified:    NO" in out


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 checkpoint not present")
def test_cli_runs_against_the_real_checkpoint(tmp_path):
    out_dir = tmp_path / "cli_out"
    proc = subprocess.run(
        [sys.executable, str(ROOT / "examples" / "run_phase8e_control_bridge.py"), "--output-dir", str(out_dir)],
        cwd=ROOT, capture_output=True, text=True, timeout=180,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "PHASE 8E" in proc.stdout
    assert (out_dir / "report.md").exists()
