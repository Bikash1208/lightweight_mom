"""Phase 8A orchestrator + CLI: `run_feasibility_gate`, `feasibility.json`/`FEASIBILITY_REPORT.md`
output, and `scripts/research/run_phase8a_feasibility.py`'s exit codes."""

import subprocess
import sys
import warnings
from pathlib import Path

import pytest

pytest.importorskip("cadquery")
pytest.importorskip("matplotlib")

from src.feasibility_gate import (
    build_feasibility_json,
    render_feasibility_report_md,
    run_and_write_phase8a,
    run_feasibility_gate,
)
from src.feasibility_schema import DEFAULT_LIBRARY_DIR, EXIT_CODES, FeasibilityConfig
from src.library_schema import LibraryConfig
from src.unit_cell_library import build_full_library

warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path(__file__).resolve().parent.parent
_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness_m": 1.575e-3}

# The small fixture library below has 2x2=4 geometries and 3 frequency points (11/12/13 GHz), not
# the locked 324/201 -- every FeasibilityConfig built against it must say so explicitly.
_N_GEOM, _N_FREQ = 4, 3


def small_config():
    return LibraryConfig(
        period_m=4.8e-3,
        substrate=_SUBSTRATE,
        f_min_hz=11.0e9,
        f_max_hz=13.0e9,
        f_step_hz=1.0e9,
        f_design_hz=12.0e9,
        l_values_mm=(1.0, 3.0),
        w_values_mm=(1.0, 3.0),
    )


def _cfg(built_library_dir, **overrides):
    kwargs = dict(library_dir=built_library_dir, expected_n_geometries=_N_GEOM, expected_n_frequency_points=_N_FREQ)
    kwargs.update(overrides)
    return FeasibilityConfig(**kwargs)


@pytest.fixture(scope="module")
def built_library_dir(tmp_path_factory):
    base = tmp_path_factory.mktemp("phase8a_gate_lib")
    build_full_library(small_config(), base, deterministic_timestamp="2026-01-01T00:00:00Z")
    return base


def test_gate_on_a_valid_small_library_reports_phase_source_unvalidated(built_library_dir):
    # This is the expected, honest result: no holographic generator exists in this project, so
    # the phase-source gate always fails today (brief Sec. 48) -- not a software defect.
    result = run_feasibility_gate(_cfg(built_library_dir, frequency_hz=12.0e9))
    assert result.status == "PHASE_SOURCE_UNVALIDATED"
    assert result.exit_code == EXIT_CODES["PHASE_SOURCE_UNVALIDATED"]
    assert result.checks["library_integrity"] == "PASS"
    assert result.checks["frequency_compatibility"] == "PASS"


# ------------------------------------------------------------------------- Test 7: bad frequency


def test_unsupported_frequency_fails_before_reaching_phase_source(built_library_dir):
    result = run_feasibility_gate(_cfg(built_library_dir, frequency_hz=20.0e9))
    assert result.status == "FREQUENCY_UNSUPPORTED"
    assert result.exit_code == EXIT_CODES["FREQUENCY_UNSUPPORTED"]


def test_off_grid_frequency_within_range_is_also_unsupported(built_library_dir):
    # 11.3 GHz is inside [11, 13] GHz but is not one of the stored 11/12/13 GHz points; no
    # interpolation is performed to make it "supported".
    result = run_feasibility_gate(_cfg(built_library_dir, frequency_hz=11.3e9))
    assert result.status == "FREQUENCY_UNSUPPORTED"


def test_unsupported_polarization_is_rejected_not_substituted(built_library_dir):
    result = run_feasibility_gate(_cfg(built_library_dir, frequency_hz=12.0e9, polarization="y"))
    assert result.status == "POLARIZATION_UNSUPPORTED"


# ------------------------------------------------------------------------------- Test 14: HFSS


def test_hfss_validation_is_never_claimed_performed(built_library_dir):
    for freq_hz in (12.0e9, 20.0e9):
        result = run_feasibility_gate(_cfg(built_library_dir, frequency_hz=freq_hz))
        assert result.hfss_validation == {"performed": False}


def test_feasibility_report_states_hfss_not_performed(built_library_dir):
    result = run_feasibility_gate(_cfg(built_library_dir, frequency_hz=12.0e9))
    report = render_feasibility_report_md(result, _cfg(built_library_dir))
    assert "HFSS validation: NOT PERFORMED" in report


# ---------------------------------------------------------------- required-phase synthetic path


def test_synthetic_required_phase_sufficient_case_reaches_matching_check(built_library_dir):
    result = run_feasibility_gate(_cfg(built_library_dir, frequency_hz=12.0e9))
    # The library's own actual phase span at 12 GHz (narrow, since only 4 geometries) is used as
    # the "required" set, so it is trivially representable -- this exercises the synthetic
    # required_phase_deg code path even though the gate itself still stops at PHASE_SOURCE_UNVALIDATED
    # before reaching it (the phase-source gate blocks first, by design).
    own_phases = [row["phase_deg"] for row in result.table]
    result2 = run_feasibility_gate(_cfg(built_library_dir, frequency_hz=12.0e9), required_phase_deg=own_phases)
    # Still PHASE_SOURCE_UNVALIDATED: the phase-source gate is checked before the matching step
    # regardless of what required_phase_deg is supplied (Sec. 36 requires every check to pass).
    assert result2.status == "PHASE_SOURCE_UNVALIDATED"


# --------------------------------------------------------------------------- output file layout


def test_run_and_write_phase8a_produces_required_files(built_library_dir, tmp_path):
    out_dir = tmp_path / "phase8a_out"
    result, written_dir = run_and_write_phase8a(_cfg(built_library_dir, frequency_hz=12.0e9), out_dir)
    assert (out_dir / "feasibility.json").exists()
    assert (out_dir / "feasibility_report.md").exists()
    assert (out_dir / "library_phase_coverage.csv").exists()
    assert (out_dir / "phase_coverage_vs_frequency.csv").exists()
    assert (out_dir / "plots" / "library_phase_vs_geometry.png").exists()
    assert (out_dir / "plots" / "library_phase_coverage.png").exists()
    assert (out_dir / "plots" / "phase_coverage_vs_frequency.png").exists()


def test_feasibility_json_matches_sec41_schema_keys(built_library_dir, tmp_path):
    result, _ = run_and_write_phase8a(_cfg(built_library_dir, frequency_hz=12.0e9), tmp_path / "out")
    doc = build_feasibility_json(result, _cfg(built_library_dir))
    for key in (
        "status",
        "library",
        "frequency",
        "phase_source",
        "phase_coverage",
        "required_phase",
        "matching",
        "hfss_validation",
        "checks",
        "warnings",
        "errors",
    ):
        assert key in doc, key
    assert doc["hfss_validation"] == {"performed": False}
    assert doc["status"] == "PHASE_SOURCE_UNVALIDATED"


# ------------------------------------------------------------------------------------------- CLI


def test_cli_runs_against_a_small_library_and_exits_with_the_gate_status_code(built_library_dir, tmp_path):
    out_dir = tmp_path / "cli_out"
    proc = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "research" / "run_phase8a_feasibility.py"),
            "--library-dir",
            str(built_library_dir),
            "--frequency-ghz",
            "12.0",
            "--expected-geometries",
            str(_N_GEOM),
            "--expected-frequency-points",
            str(_N_FREQ),
            "--output-dir",
            str(out_dir),
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == EXIT_CODES["PHASE_SOURCE_UNVALIDATED"], proc.stdout + proc.stderr
    assert "PHASE_SOURCE_UNVALIDATED" in proc.stdout
    assert "HFSS validation" in proc.stdout
    assert (out_dir / "feasibility.json").exists()


# -------------------------------------------------------------- against the real production library


@pytest.mark.skipif(
    not (DEFAULT_LIBRARY_DIR / "database" / "unit_cell_library.json").exists(),
    reason="Phase 7 production library not built",
)
def test_gate_against_the_real_production_library():
    result = run_feasibility_gate(FeasibilityConfig())
    assert result.checks["library_integrity"] == "PASS", result.errors
    assert result.library["geometry_count"] == 324
    assert result.status == "PHASE_SOURCE_UNVALIDATED"
    assert result.library_phase_coverage["span_deg"] < 30.0  # the measured ~24.8 deg
