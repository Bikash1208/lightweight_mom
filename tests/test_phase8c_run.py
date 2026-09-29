"""Phase 8C end-to-end pipeline and CLI, against the real Phase 7 checkpoint data.

Every test except the CLI subprocess test passes in the session-scoped `phase8c_full_database`
fixture (`conftest.py`) so the ~20 s checkpoint parse happens once for the whole test session, not
once per test.
"""

import subprocess
import sys
import warnings
from pathlib import Path

import pytest

pytest.importorskip("matplotlib")

from src.phase8c import data_source
from src.phase8c.run import print_final_report, run_phase8c

warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path(__file__).resolve().parent.parent
_HAS_CHECKPOINT = data_source.DEFAULT_CHECKPOINT_JSON.exists()


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 internal solve checkpoint not present")
def test_run_phase8c_produces_all_required_files(tmp_path, phase8c_full_database):
    result = run_phase8c(tmp_path / "out", database=phase8c_full_database)
    out = result["output_dir"]
    for name in (
        "unit_cell_topology_audit.md",
        "s_parameter_definition.md",
        "control_quantity_analysis.md",
        "control_quantity_validation.json",
        "unit_cell_scattering_12p5GHz.csv",
        "control_quantity_candidates.csv",
        "phase_state_analysis.csv",
    ):
        assert (out / name).exists(), name
    for i in range(1, 19):
        matches = list((out / "plots").glob(f"{i:02d}_*.png"))
        assert matches, f"plot {i:02d} missing"


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 internal solve checkpoint not present")
def test_run_phase8c_status_is_control_variable_unvalidated(tmp_path, phase8c_full_database):
    result = run_phase8c(tmp_path / "out", database=phase8c_full_database)
    assert result["status"] == "CONTROL_VARIABLE_UNVALIDATED"
    assert result["audit"].ground_plane == "ABSENT"
    assert result["reciprocity"].reciprocal


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 internal solve checkpoint not present")
def test_run_phase8c_does_not_modify_phase7_checkpoint_or_delivered_library(tmp_path, phase8c_full_database):
    checkpoint_before = data_source.DEFAULT_CHECKPOINT_JSON.stat().st_size
    delivered = ROOT / "results" / "unit_cell_library" / "database" / "unit_cell_library.json"
    delivered_before = delivered.read_bytes()
    run_phase8c(tmp_path / "out", database=phase8c_full_database)
    assert data_source.DEFAULT_CHECKPOINT_JSON.stat().st_size == checkpoint_before
    assert delivered.read_bytes() == delivered_before


@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 internal solve checkpoint not present")
def test_print_final_report_never_claims_hfss_performed_or_bare_pass(tmp_path, capsys, phase8c_full_database):
    result = run_phase8c(tmp_path / "out", database=phase8c_full_database)
    print_final_report(result)
    out = capsys.readouterr().out
    assert "HFSS validation: NOT PERFORMED" in out
    assert "Solver modified: NO" in out
    assert "Phase 7 modified: NO" in out
    final_line = out.split("FINAL SCIENTIFIC STATUS")[1].strip().splitlines()
    assert final_line[1].strip() != "PASS"


@pytest.mark.slow
@pytest.mark.skipif(not _HAS_CHECKPOINT, reason="Phase 7 internal solve checkpoint not present")
def test_cli_runs_against_the_real_checkpoint(tmp_path):
    out_dir = tmp_path / "cli_out"
    proc = subprocess.run(
        [sys.executable, str(ROOT / "examples" / "run_phase8c_validation.py"), "--output-dir", str(out_dir)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "PHASE 8C" in proc.stdout
    assert "CONTROL_VARIABLE_UNVALIDATED" in proc.stdout
    assert (out_dir / "control_quantity_validation.json").exists()
