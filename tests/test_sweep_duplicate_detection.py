"""Phase 12B: SweepPlan.mesh_duplicate_groups / n_mesh_distinct and plan_sweep's warn_duplicates."""

import warnings

import pytest

from src.config import config_from_dict
from src.runner import run_sweep
from src.sweep import ParameterAxis, SweepDefinition, plan_sweep

BASE = {
    "period_x": 4.8e-3,
    "period_y": 4.8e-3,
    "substrate": {"name": "RT5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3},
    "geometry": {"type": "rectangle", "length": 2.0e-3, "width": 2.0e-3},
    "excitation": {"frequency": 12.5e9},
    "mesh": {"cells_x": 32, "cells_y": 32, "floquet_order": 32},
    "reference_planes": {"top_air_thickness": 5.0e-3, "bottom_air_thickness": 5.0e-3},
}


def _config(**mesh_override):
    import copy

    data = copy.deepcopy(BASE)
    data["mesh"].update(mesh_override)
    return config_from_dict(data)


def test_plan_detects_the_known_32x32_duplicate_pairs():
    cfg = _config()
    defn = SweepDefinition(axes=(ParameterAxis("length_m", tuple(l * 1e-3 for l in (1.4, 1.6, 2.0, 2.2, 2.8))),))
    with pytest.warns(UserWarning, match="mesh-identical"):
        plan = plan_sweep(cfg, defn)
    assert len(plan.mesh_duplicate_groups) == 2
    assert plan.n_mesh_distinct == 3  # 5 cells, 2 duplicate pairs -> 3 distinct


def test_no_duplicates_and_no_warning_on_the_aligned_48_cell_mesh():
    cfg = _config(cells_x=48, cells_y=48, floquet_order=48)
    defn = SweepDefinition(axes=(ParameterAxis("length_m", tuple(l * 1e-3 for l in (1.4, 1.6, 2.0, 2.2, 2.8))),))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        plan = plan_sweep(cfg, defn)
    assert plan.mesh_duplicate_groups == ()
    assert plan.n_mesh_distinct == 5


def test_warn_duplicates_false_suppresses_the_warning_but_keeps_the_data():
    cfg = _config()
    defn = SweepDefinition(axes=(ParameterAxis("length_m", tuple(l * 1e-3 for l in (1.4, 1.6))),))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        plan = plan_sweep(cfg, defn, warn_duplicates=False)
    assert plan.mesh_duplicate_groups == ((0, 1),)


def test_run_sweep_never_warns_about_duplicates_and_records_them_in_metadata(tmp_path):
    cfg = _config()
    defn = SweepDefinition(axes=(ParameterAxis("length_m", tuple(l * 1e-3 for l in (1.4, 1.6))),))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = run_sweep(cfg, out_dir=tmp_path, definition=defn)
    assert result.database.metadata["mesh_duplicate_groups"] == [[0, 1]]
    assert result.database.metadata["n_mesh_distinct"] == 1


def test_default_construction_has_no_duplicate_groups():
    """SweepPlan built without the new keyword defaults to empty (backward compatible field)."""
    from src.sweep import SweepPlan

    plan = SweepPlan(config=None, definition=None, frequencies=None, context={}, cells=())
    assert plan.mesh_duplicate_groups == ()
