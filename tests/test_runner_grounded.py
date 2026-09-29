"""Phase 12A: ``run_sweep`` end-to-end for a grounded (one-port) config -- checkpoint resume,
refinement study, and the CSV-skip (JSON-only) behaviour for a one-port database.
"""

import copy

from src.config import config_from_dict
from src.database import Database
from src.sweep import ParameterAxis, RefinementPolicy, SweepDefinition, ValidationPolicy
from src.runner import run_sweep

BASE = {
    "period_x": 4.8e-3, "period_y": 4.8e-3,
    "substrate": {"name": "RT5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3, "backing": "pec"},
    "geometry": {"type": "rectangle", "length": 2.0e-3, "width": 2.0e-3},
    "excitation": {"frequency": 12.5e9},
    "mesh": {"cells_x": 8, "cells_y": 8, "floquet_order": 16},
    "reference_planes": {"top_air_thickness": 5.0e-3},
}


def _config():
    return config_from_dict(copy.deepcopy(BASE))


def test_grounded_sweep_produces_one_port_records_only(tmp_path):
    cfg = _config()
    defn = SweepDefinition(axes=(ParameterAxis("length_m", (1.5e-3, 2.0e-3, 2.5e-3)),))
    result = run_sweep(cfg, out_dir=tmp_path, definition=defn)
    assert len(result.database.records) == 3
    assert all(r.port_topology == "one_port" for r in result.database.records)
    assert all(r.s.shape == (2, 2) for r in result.database.records)


def test_grounded_sweep_writes_json_but_not_csv(tmp_path):
    cfg = _config()
    defn = SweepDefinition(axes=(ParameterAxis("length_m", (1.5e-3, 2.0e-3)),))
    result = run_sweep(cfg, out_dir=tmp_path, definition=defn)
    assert result.json_path is not None and result.json_path.exists()
    assert result.csv_path is None
    assert not (tmp_path / "database.csv").exists()


def test_grounded_sweep_checkpoint_resume_reuses_every_result(tmp_path):
    cfg = _config()
    defn = SweepDefinition(axes=(ParameterAxis("length_m", (1.5e-3, 2.0e-3, 2.5e-3)),))
    first = run_sweep(cfg, out_dir=tmp_path, definition=defn)
    assert first.computed == 3 and first.reused == 0
    second = run_sweep(cfg, out_dir=tmp_path, definition=defn)
    assert second.computed == 0 and second.reused == 3
    assert second.database == first.database


def test_grounded_sweep_json_round_trip_equals_in_memory_database(tmp_path):
    cfg = _config()
    defn = SweepDefinition(axes=(ParameterAxis("length_m", (1.5e-3, 2.0e-3)),))
    result = run_sweep(cfg, out_dir=tmp_path, definition=defn)
    reloaded = Database.from_json(result.json_path)
    assert reloaded == result.database


def test_grounded_sweep_with_refinement_study_enabled(tmp_path):
    cfg = _config()
    defn = SweepDefinition(
        axes=(ParameterAxis("length_m", (2.0e-3,)),),
        validation=ValidationPolicy(refinement=RefinementPolicy(enabled=True, factor=2, tolerance=0.5)),
    )
    result = run_sweep(cfg, out_dir=tmp_path, definition=defn)
    validation = result.database.cell(0).validation
    assert validation["refinement"] is not None
    assert "floquet_error" in validation["refinement"]
    assert "mesh_error" in validation["refinement"]


def test_grounded_sweep_records_per_cell_surface_wave_proximity(tmp_path):
    cfg = _config()
    defn = SweepDefinition(axes=(ParameterAxis("length_m", (2.0e-3,)),))
    result = run_sweep(cfg, out_dir=tmp_path, definition=defn)
    sw = result.database.cell(0).validation["surface_wave"]
    assert sw is not None
    assert sw["any_in_band"] is False
    assert "frequency_hz" in sw


def test_ungrounded_sweep_has_no_surface_wave_entry(tmp_path):
    import copy

    data = copy.deepcopy(BASE)
    del data["substrate"]["backing"]
    data["reference_planes"] = {"top_air_thickness": 5.0e-3, "bottom_air_thickness": 5.0e-3}
    cfg = config_from_dict(data)
    defn = SweepDefinition(axes=(ParameterAxis("length_m", (2.0e-3,)),))
    result = run_sweep(cfg, out_dir=tmp_path, definition=defn)
    assert result.database.cell(0).validation["surface_wave"] is None


def test_diagnostics_reciprocity_and_passivity_for_grounded_cell(tmp_path):
    cfg = _config()
    defn = SweepDefinition(axes=(ParameterAxis("length_m", (2.0e-3,)),))
    result = run_sweep(cfg, out_dir=tmp_path, definition=defn)
    diag = result.database.records[0].diagnostics
    assert diag["reciprocity_error"] < 1e-8
    assert diag["passivity_min_eig"] > -1e-8   # lossy stack: strictly passive
