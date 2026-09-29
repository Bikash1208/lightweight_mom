"""Phase 12C end-to-end orchestration, at a deliberately reduced scale (real solves, tiny sizes)
so this test stays fast while still exercising the full A-E pipeline and every output file."""

import json

import pytest

from src.phase12c.run import run_phase12c


@pytest.fixture(scope="module")
def reduced_result(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("phase12c_reduced") / "phase12c_grounded_locality"
    return run_phase12c(
        out_dir,
        correctness_mesh_sequence=(2, 4),
        feasibility_sizes=(1, 3),
        feasibility_mesh_cells_per_subcell=4,
        periodic_n_cells=4,
        coupling_sizes=(1, 3),
        spot_check_mesh_cells_per_subcell=4,
        export_mesh_cells_per_side=8,   # 4 is too coarse to resolve L=1.6 mm (SINGLE_CELL_L_MM's
                                         # smallest value) at P=8 mm -- see test_phase12c_exports.py
    )


def test_reduced_run_produces_every_output_file(reduced_result):
    out = reduced_result.output_dir
    for name in ("task_a.json", "task_b.json", "task_c.json", "task_d.json",
                 "task_d_spot_check.csv", "phase12c_gate_status.json", "phase12c_report.md"):
        assert (out / name).exists(), name


def test_reduced_run_gate_status_json_has_every_required_status(reduced_result):
    data = json.loads((reduced_result.output_dir / "phase12c_gate_status.json").read_text())
    assert data["hfss_status"] == "HFSS_NOT_AVAILABLE"
    assert data["geometry_inversion"] == "GEOMETRY_INVERSION_LOCKED"
    assert data["step_export"] == "STEP_EXPORT_LOCKED"
    assert "task_a_status" in data and "task_b_local_periodicity" in data


def test_reduced_run_task_a_correctness_is_near_machine_precision(reduced_result):
    assert reduced_result.correctness.relative_error < 1e-6


def test_reduced_run_task_b_has_three_centers_and_a_periodic_supercell(reduced_result):
    assert reduced_result.centers.off_resonance != reduced_result.centers.at_resonance
    assert len(reduced_result.periodic_results) == 4


def test_reduced_run_task_d_gate_has_all_four_findings(reduced_result):
    names = [f[0] for f in reduced_result.task_d_gate.findings]
    assert set(names) == {"phase10_c2_c3", "phase9_strong_contrast", "phase11_case_d", "phase11_coupling_radius"}


def test_reduced_run_report_mentions_every_task(reduced_result):
    text = (reduced_result.output_dir / "phase12c_report.md").read_text()
    for marker in ("Task A", "Task B", "Task C", "Task D", "Task E"):
        assert marker in text
