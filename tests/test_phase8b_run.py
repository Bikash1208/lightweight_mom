"""Phase 8B end-to-end pipeline and CLI: Sec. 50's output directory layout, Sec. 55's acceptance
checklist items that are only checkable against a real run, and Sec. 56's terminal report."""

import subprocess
import sys
import warnings
from pathlib import Path

import pytest

pytest.importorskip("matplotlib")

from src.feasibility_schema import DEFAULT_LIBRARY_DIR
from src.phase8b.run import print_final_report, run_phase8b

warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path(__file__).resolve().parent.parent
_HAS_PRODUCTION_LIBRARY = (DEFAULT_LIBRARY_DIR / "database" / "unit_cell_library.json").exists()


@pytest.mark.skipif(not _HAS_PRODUCTION_LIBRARY, reason="Phase 7 production library not built")
def test_run_phase8b_produces_all_required_files(tmp_path):
    out_dir = tmp_path / "phase8b_out"
    result = run_phase8b(out_dir)
    for name in (
        "repository_scientific_audit.md",
        "Holographic_Formulation.md",
        "holographic_specification.json",
        "control_variable_analysis.md",
        "source_registry.json",
        "phase_control_quantity.csv",
        "phase8b_distinct_state_analysis.csv",
    ):
        assert (out_dir / name).exists(), name
    for name in (
        "object_wave_single_beam.png",
        "object_wave_multi_beam.png",
        "phase_gradient_test.png",
        "phase_control_quantity.png",
        "distinct_phase_states.png",
    ):
        assert (out_dir / "plots" / name).exists(), name


@pytest.mark.skipif(not _HAS_PRODUCTION_LIBRARY, reason="Phase 7 production library not built")
def test_run_phase8b_status_is_phase_source_unvalidated(tmp_path):
    result = run_phase8b(tmp_path / "out")
    assert result["status"] == "PHASE_SOURCE_UNVALIDATED"
    assert result["distinct_states"].n_geometries == 324


@pytest.mark.skipif(not _HAS_PRODUCTION_LIBRARY, reason="Phase 7 production library not built")
def test_run_phase8b_does_not_modify_phase7_files(tmp_path):
    json_path = DEFAULT_LIBRARY_DIR / "database" / "unit_cell_library.json"
    before = json_path.read_bytes()
    run_phase8b(tmp_path / "out")
    after = json_path.read_bytes()
    assert before == after


@pytest.mark.skipif(not _HAS_PRODUCTION_LIBRARY, reason="Phase 7 production library not built")
def test_print_final_report_never_claims_hfss_performed(tmp_path, capsys):
    result = run_phase8b(tmp_path / "out")
    print_final_report(result)
    out = capsys.readouterr().out
    assert "NOT PERFORMED" in out
    assert "Final scientific status:" in out
    final_line = out.split("Final scientific status:")[1].splitlines()[0].strip()
    assert final_line != "PASS"
    assert "Solver modified:" in out and "NO" in out


@pytest.mark.skipif(not _HAS_PRODUCTION_LIBRARY, reason="Phase 7 production library not built")
def test_cli_runs_against_the_real_production_library(tmp_path):
    out_dir = tmp_path / "cli_out"
    proc = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "research" / "run_phase8b_formulation.py"),
            "--output-dir",
            str(out_dir),
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "PHASE 8B" in proc.stdout
    assert "HFSS validation:" in proc.stdout and "NOT PERFORMED" in proc.stdout
    assert (out_dir / "holographic_specification.json").exists()
