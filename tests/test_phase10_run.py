"""Phase 10 end-to-end pipeline, at REDUCED scale (2 centers, a handful of random samples) so the
routine test suite stays fast -- every stage is still a REAL MoM solve, nothing mocked. The full,
production-scale run (all 5 centers, 350 random samples) is performed separately via
``examples/run_phase10_environment_model.py`` and its real, measured result is what
``PHASE10_ENVIRONMENT_MODEL.md`` reports.
"""

import json

import pytest

from src.phase10 import run as run_mod


@pytest.fixture(scope="module")
def reduced_run(tmp_path_factory):
    """Module-scoped: `run_phase10` is a real, multi-minute-if-run-repeatedly MoM pipeline, so it
    runs ONCE and every test below shares the result (mirrors `phase8c_full_database`'s own
    session-scoped-fixture rationale in `tests/conftest.py`)."""
    import pytest as _pytest

    mp = _pytest.MonkeyPatch()
    try:
        mp.setattr(run_mod, "N_RANDOM_SAMPLES", 10)
        reduced_centers = {"C3": (2.8, 2.8), "C5": (4.2, 4.2)}
        # a single module-object patch: every Phase 10 submodule holds the SAME
        # `environment_geometry` module instance under its own `eg` alias, so one attribute patch
        # is visible everywhere.
        mp.setattr("src.phase10.environment_geometry.CENTER_CELLS_MM", reduced_centers)
        result = run_mod.run_phase10(tmp_path_factory.mktemp("phase10_reduced_run") / "out")
    finally:
        mp.undo()
    return result


def test_run_produces_all_required_top_level_outputs(reduced_run):
    out = reduced_run["output_dir"]
    for name in ("phase10_report.md", "phase10_summary.json", "phase10_gate_status.json"):
        assert (out / name).exists(), name


def test_run_produces_all_required_subdirectory_outputs(reduced_run):
    out = reduced_run["output_dir"]
    for name in ("five_cell_cases.csv", "diagonal_cases.csv", "long_range_cases.csv", "neighborhood_decision.json"):
        assert (out / "neighborhood_study" / name).exists(), name
    for name in ("train.csv", "validation.csv", "blind_test.csv", "dataset_metadata.json"):
        assert (out / "environment_dataset" / name).exists(), name
    for name in ("model_metadata.json", "predictions.csv", "error_metrics.csv"):
        assert (out / "surrogate" / name).exists(), name
    assert (out / "surrogate" / "model" / "surrogate.json").exists()
    for name in ("target_responses.csv", "inverse_results.csv", "mom_verification.csv"):
        assert (out / "inverse_validation" / name).exists(), name
    assert (out / "beam_validation" / "beam_validation.json").exists()


def test_run_produces_the_required_plot_counts(reduced_run):
    out = reduced_run["output_dir"]
    assert len(list((out / "neighborhood_study" / "plots").glob("*.png"))) >= 5
    assert len(list((out / "surrogate" / "plots").glob("*.png"))) >= 8
    assert len(list((out / "beam_validation" / "plots").glob("*.png"))) >= 7


def test_summary_json_is_well_formed_and_matches_decision(reduced_run):
    out = reduced_run["output_dir"]
    doc = json.loads((out / "phase10_summary.json").read_text())
    assert doc["final_status"] == reduced_run["decision"].final_status
    assert doc["gates"]["10A_control_variable"] == "CLOSED"   # no HFSS reference in this repo


def test_run_never_unlocks_geometry_inversion_without_hfss_data(reduced_run):
    """Sec. 47's mandatory negative test at the full pipeline level: with no HFSS/measurement
    reference anywhere in this repository, geometry inversion must stay LOCKED regardless of how
    every other stage performs."""
    decision = reduced_run["decision"]
    assert decision.gate_10a.status == "CLOSED"
    assert decision.geometry_inversion_unlocked is False


def test_run_honestly_reports_no_hfss_reference(reduced_run):
    assert reduced_run["hfss_exists"] is False


def test_solver_and_earlier_phases_are_not_imported_as_modified_copies():
    """Sanity: Phase 10 imports Phase 1-9 modules directly (no vendored/forked copies)."""
    import src.phase9.supercell_solver as p9_solver
    import src.solver as solver_mod
    from src.phase10 import grid_solver as p10_grid_solver

    assert p10_grid_solver.solver_mod is solver_mod
    assert p9_solver.solver_mod is solver_mod
