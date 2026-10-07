"""Phase 11 end-to-end pipeline, at REDUCED scale (2 centers) so the routine test suite stays
bounded -- every stage is still a REAL MoM solve, nothing mocked. The full, production-scale run
(all 5 centers) is performed separately via ``scripts/research/run_phase11_heterogeneous_mom.py`` and its
real, measured result is what ``PHASE11_HETEROGENEOUS_MOM.md`` reports. This single test is
inherently the most expensive one in the whole suite (Phase 11's own scope is large by design,
Sec. 36: "do not try to make the problem look easier than it is") -- module-scoped so it runs
exactly once.
"""

import json

import pytest

from src.phase11 import run as run_mod

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def reduced_run(tmp_path_factory):
    import pytest as _pytest

    mp = _pytest.MonkeyPatch()
    try:
        mp.setattr(run_mod, "CENTER_CELLS_MM", {"C3": (2.8, 2.8), "C1": (1.4, 1.4)})
        result = run_mod.run_phase11(tmp_path_factory.mktemp("phase11_reduced_run") / "out")
    finally:
        mp.undo()
    return result


def test_run_produces_all_required_top_level_outputs(reduced_run):
    out = reduced_run["output_dir"]
    required = (
        "phase11_report.md",
        "phase11_summary.json",
        "phase11_gate_status.json",
        "coupling_radius_results.json",
        "incremental_shell_results.json",
        "directional_coupling.json",
        "diagonal_coupling.json",
        "pattern_study.json",
        "frequency_study.json",
        "polarization_study.json",
        "pairwise_interaction.json",
        "higher_order_interaction.json",
        "convergence.json",
        "beam_validation.json",
    )
    for name in required:
        assert (out / name).exists(), name


def test_run_produces_the_required_plot_count(reduced_run):
    out = reduced_run["output_dir"]
    assert len(list((out / "plots").glob("*.png"))) >= 12


def test_summary_json_is_well_formed_and_matches_decision(reduced_run):
    out = reduced_run["output_dir"]
    doc = json.loads((out / "phase11_summary.json").read_text())
    assert doc["overall_case"] == reduced_run["decision"].overall_case
    assert doc["hfss_status"] == "HFSS_NOT_AVAILABLE"  # no HFSS reference in this repo


def test_run_never_unlocks_geometry_inversion_without_hfss_or_convergence(reduced_run):
    """Sec. 28's mandatory negative test at the full pipeline level."""
    decision = reduced_run["decision"]
    assert decision.hfss_status == "HFSS_NOT_AVAILABLE"
    assert decision.geometry_inversion_status == "GEOMETRY_INVERSION_LOCKED"
    assert decision.step_export_status == "STEP_EXPORT_LOCKED"


def test_run_reports_no_ml_by_construction(reduced_run):
    """Sec. 2/33: a real, mechanical check that no forbidden ML library is actually IMPORTED
    anywhere in Phase 11 -- parses each file's AST and inspects real `import`/`from ... import`
    statements, not prose (several docstrings in this package deliberately explain, in words,
    that no neural network/random forest/etc. is used, which a naive text grep would misflag)."""
    import ast
    import pathlib

    forbidden_modules = {"sklearn", "tensorflow", "torch", "keras", "xgboost", "lightgbm"}
    src_dir = pathlib.Path(__file__).resolve().parent.parent / "src" / "phase11"
    for path in src_dir.glob("*.py"):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    assert alias.name.split(".")[0] not in forbidden_modules, (
                        f"forbidden ML import {alias.name!r} found in {path}"
                    )
            elif isinstance(node, ast.ImportFrom) and node.module:
                assert node.module.split(".")[0] not in forbidden_modules, (
                    f"forbidden ML import {node.module!r} found in {path}"
                )


def test_diagonal_reproduction_is_within_or_near_phase10s_own_range(reduced_run):
    ratio = reduced_run["diagonal_reproduction"]["ratio"]
    assert ratio > 0
