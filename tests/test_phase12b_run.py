"""Phase 12B end-to-end pipeline, at REDUCED scale, so the routine test suite stays bounded --
every stage is still a REAL MoM solve, nothing mocked. The full, production-scale run is
`scripts/research/run_phase12b.py`, and its real, measured result is what `PHASE12B_MESH_PERIOD_STUDY.md`
reports. Module-scoped so the reduced run happens once per test session
(``tests/test_phase11_run.py``/``test_phase12a_run.py``'s own established pattern).
"""

import json

import pytest

from src.phase12b.run import run_phase12b

pytestmark = pytest.mark.slow

REDUCED_L_MM = (1.0, 2.0, 3.0, 4.4)
REDUCED_FREQUENCIES_HZ = (12.5e9,)
REDUCED_PERIOD_L_MM = {8.0: (2.0, 4.0, 6.0)}


@pytest.fixture(scope="module")
def reduced_run(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("phase12b_reduced_run") / "out"
    return run_phase12b(
        out_dir,
        l_values_mm=REDUCED_L_MM,
        frequencies_hz=REDUCED_FREQUENCIES_HZ,
        new_periods_mm=(8.0,),
        period_l_values_mm=REDUCED_PERIOD_L_MM,
        aligned_refinement_cells=(24, 48),
    )


def test_run_produces_all_required_top_level_outputs(reduced_run):
    out = reduced_run.output_dir
    required = (
        "phase12b_report.md",
        "phase12b_gate_status.json",
        "alignment_report.json",
        "reassessment.json",
        "duplicate_audit.json",
        "period_study.json",
        "period_sweep.csv",
    )
    for name in required:
        assert (out / name).exists(), name


def test_run_produces_the_required_plots(reduced_run):
    out = reduced_run.output_dir
    plots = {p.name for p in (out / "plots").glob("*.png")}
    assert plots == {
        "aligned_vs_unaligned_convergence.png",
        "phase_vs_length_per_period.png",
        "coverage_vs_frequency_per_period.png",
        "magnitude_vs_length_per_period.png",
        "duplicate_audit_table.png",
    }


def test_gate_status_json_matches_the_decision(reduced_run):
    doc = json.loads((reduced_run.output_dir / "phase12b_gate_status.json").read_text())
    assert doc["mesh_alignment"]["status"] == reduced_run.mesh_gate.status
    assert doc["audit"]["status"] == reduced_run.audit_gate.status
    assert doc["period_study"]["status"] == reduced_run.period_gate.status


def test_run_never_unlocks_geometry_inversion_or_step_export(reduced_run):
    text = (reduced_run.output_dir / "phase12b_report.md").read_text()
    assert "GEOMETRY_INVERSION_LOCKED" in text
    assert "STEP_EXPORT_LOCKED" in text


def test_audit_finds_the_known_phase10_c2_c3_collision(reduced_run):
    assert reduced_run.audit_gate.status == "CONCLUSIONS_AT_RISK"
    assert any("C2" in f and "C3" in f for f in reduced_run.audit_gate.at_risk_findings)


def test_period_results_cover_4p8_and_8_mm(reduced_run):
    assert set(reduced_run.period_results) == {4.8, 8.0}
