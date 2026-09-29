"""Phase 12D end-to-end orchestration, at a deliberately reduced scale (real solves, tiny sizes)
so this test stays fast while still exercising the full A-C pipeline and every output file."""

import json

import pytest

from src.phase12d.run import run_phase12d


def _small_step_family(mesh_cells_per_subcell):
    unit_mm = 2.0 * 8.0 / mesh_cells_per_subcell
    max_k = int((8.0 * 0.9) // unit_mm)  # stay safely inside the 8mm period
    return tuple(round(unit_mm * k, 6) for k in range(2, max_k + 1))


@pytest.fixture(scope="module")
def reduced_result(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("phase12d_reduced") / "phase12d_local_periodicity"
    return run_phase12d(
        out_dir,
        max_k=1,
        step_curve_mesh_cells_per_subcell=12,
        step_curve_l_family_mm=_small_step_family(12),
        aperture_sizes=(3,),
        aperture_mesh_cells_per_subcell=12,
        aperture_l_family_mm=(2.0, 4.0, 6.0),
    )


def test_reduced_run_produces_every_output_file(reduced_result):
    out = reduced_result.output_dir
    for name in ("task_a.json", "task_b.json", "task_c.json", "phase12d_gate_status.json", "phase12d_report.md"):
        assert (out / name).exists(), name


def test_reduced_run_gate_status_json_has_every_required_status(reduced_result):
    data = json.loads((reduced_result.output_dir / "phase12d_gate_status.json").read_text())
    assert data["geometry_inversion"] == "GEOMETRY_INVERSION_LOCKED"
    assert data["step_export"] == "STEP_EXPORT_LOCKED"
    assert "task_a_step_limit" in data and "task_b_aperture" in data and "task_c_hfss" in data


def test_reduced_run_task_a_has_step_cases_and_centers(reduced_result):
    assert reduced_result.centers.off_resonance != reduced_result.centers.at_resonance
    assert len(reduced_result.step_cases) > 0


def test_reduced_run_task_b_has_aperture_result(reduced_result):
    assert len(reduced_result.aperture_results) == 1
    assert reduced_result.aperture_results[0].n_cells == 3


def test_reduced_run_report_mentions_every_task(reduced_result):
    text = (reduced_result.output_dir / "phase12d_report.md").read_text()
    for marker in ("Task A", "Task B", "Task C"):
        assert marker in text
