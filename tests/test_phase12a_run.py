"""Phase 12A end-to-end pipeline, at REDUCED scale (fewer L values/frequencies, a smaller
convergence-size list) so the routine test suite stays bounded -- every stage is still a REAL MoM
solve, nothing mocked. The full, production-scale run is `scripts/research/run_phase12a_ground_plane.py`,
and its real, measured result is what `PHASE12A_GROUND_PLANE.md` reports. Module-scoped so the
(still nontrivial) reduced run happens once per test session, matching
``tests/test_phase11_run.py``'s own established pattern.
"""

import json

import pytest

from src.phase12a.run import run_phase12a

pytestmark = pytest.mark.slow

REDUCED_L_MM = (1.0, 2.0, 3.0, 4.4)
REDUCED_FREQUENCIES_HZ = (12.5e9,)  # must include DESIGN_FREQUENCY_HZ (12.5 GHz)


@pytest.fixture(scope="module")
def reduced_run(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("phase12a_reduced_run") / "out"
    return run_phase12a(
        out_dir,
        l_values_mm=REDUCED_L_MM,
        frequencies_hz=REDUCED_FREQUENCIES_HZ,
        convergence_sizes=(8, 16),
        skip_ungrounded_comparison=True,
    )


def test_run_produces_all_required_top_level_outputs(reduced_run):
    out = reduced_run.output_dir
    required = (
        "phase12a_report.md",
        "phase12a_gate_status.json",
        "sweep.csv",
        "sweep.json",
        "ungrounded_comparison.json",
        "analytic_checks.json",
        "surface_wave_pole.json",
        "convergence.json",
    )
    for name in required:
        assert (out / name).exists(), name


def test_run_produces_the_required_plots(reduced_run):
    out = reduced_run.output_dir
    plots = {p.name for p in (out / "plots").glob("*.png")}
    assert plots == {
        "phase_vs_length.png",
        "magnitude_vs_length.png",
        "grounded_vs_ungrounded.png",
        "convergence.png",
    }


def test_gate_status_json_is_well_formed_and_matches_decision(reduced_run):
    doc = json.loads((reduced_run.output_dir / "phase12a_gate_status.json").read_text())
    assert doc["status"] == reduced_run.decision.status
    assert doc["geometry_inversion_status"] == "GEOMETRY_INVERSION_LOCKED"
    assert doc["step_export_status"] == "STEP_EXPORT_LOCKED"


def test_run_never_unlocks_geometry_inversion_or_step_export(reduced_run):
    """Sec. 28-style mandatory negative check at the full pipeline level (every other phase's
    own convention): this phase never attempts to unlock either, regardless of the measured
    physics outcome."""
    d = reduced_run.decision
    assert d.geometry_inversion_status == "GEOMETRY_INVERSION_LOCKED"
    assert d.step_export_status == "STEP_EXPORT_LOCKED"


def test_run_status_is_one_of_the_fixed_set(reduced_run):
    assert reduced_run.decision.status in ("GROUNDED_STACK_VALIDATED", "GROUNDED_STACK_PARTIAL")


def test_sweep_json_has_one_point_per_l_value(reduced_run):
    doc = json.loads((reduced_run.output_dir / "sweep.json").read_text())
    assert len(doc["points"]) == len(REDUCED_L_MM) * len(REDUCED_FREQUENCIES_HZ)


def test_ungrounded_comparison_reports_not_available_when_skipped(reduced_run):
    doc = json.loads((reduced_run.output_dir / "ungrounded_comparison.json").read_text())
    assert doc["status"] == "NOT_AVAILABLE_FROM_CURRENT_DATA"


def test_report_contains_no_fabricated_hfss_claim(reduced_run):
    text = (reduced_run.output_dir / "phase12a_report.md").read_text()
    assert "HFSS_VALIDATED" not in text or "HFSS_NOT_AVAILABLE" in text
    assert "HFSS_NOT_AVAILABLE" in text
