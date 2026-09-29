"""Running sweeps: frequency consistency, determinism, resume/restart, phase unwrapping on real
data, and physical checks of the stored linear/circular S-parameters."""

import dataclasses
import json
import os
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
import yaml

import src.runner as runner
from src.config import config_from_dict
from src.constants import C0
from src.database import Database
from src.green_function import LayeredGreenFunction
from src.phase import max_wrapped_step_deg, unwrap_deg, wrap_deg
from src.polarization import E_R, block_to_circular, matrix4_to_circular
from src.runner import CHECKPOINT_NAME, CheckpointError, run_sweep
from src.sparameters import compute_s_parameters, stack_s_parameters
from sweep_helpers import config_dict, sweep_config

ROOT = Path(__file__).resolve().parent.parent

SWEEP = {
    "frequency": {"start": "8 GHz", "stop": "12 GHz", "step": "2 GHz"},
    "parameters": {"length": {"values": ["3 mm", "5 mm"]}},
}
N_CELLS, N_FREQ = 2, 3


def rot(psi_deg):
    c, s = np.cos(np.radians(psi_deg)), np.sin(np.radians(psi_deg))
    return np.array([[c, -s], [s, c]])


def with_length(cfg, length):
    g = dataclasses.replace(cfg.unit_cell.geometry, length=length)
    return dataclasses.replace(cfg.unit_cell, geometry=g)


class SolveCounter:
    """Counts calls to the solver as seen by the runner."""

    def __init__(self, monkeypatch):
        self.n = 0
        original = runner.compute_s_parameters

        def spy(*args, **kwargs):
            self.n += 1
            return original(*args, **kwargs)

        monkeypatch.setattr(runner, "compute_s_parameters", spy)


# ==================================================== frequency sweep consistency

def test_every_record_equals_an_independent_single_frequency_solve():
    cfg = sweep_config(SWEEP)
    db = run_sweep(cfg).database
    for cell, length in zip(db.cells, (3e-3, 5e-3)):
        uc = with_length(cfg, length)                     # built without any sweep code
        for rec in db.records_of(cell.index):
            fresh = compute_s_parameters(uc, rec.frequency, cfg.mesh, cfg.reference_planes, warn_grating=False)
            assert np.max(np.abs(fresh.matrix() - rec.s)) < 1e-13
            assert rec.diagnostics["incoming_residual"] == pytest.approx(fresh.incoming_residual, abs=1e-14)


def test_frequency_grid_is_stored_as_requested_and_records_follow_it():
    db = run_sweep(sweep_config(SWEEP)).database
    assert db.frequencies.tolist() == [8e9, 10e9, 12e9]
    assert [r.frequency_index for r in db.records_of(0)] == [0, 1, 2]
    assert db.metadata["frequency_grid"] == {
        "start_hz": 8e9, "stop_hz": 12e9, "step_hz": 2e9, "count": 3, "frequencies_hz": [8e9, 10e9, 12e9]}


def test_results_do_not_depend_on_the_other_frequencies_in_the_sweep():
    """Sub-grid, different step and reversed order (single-frequency sweeps) give bit-identical
    values at the shared frequencies: no state leaks between frequencies."""
    fine = run_sweep(sweep_config({**SWEEP, "frequency": {"start": "8 GHz", "stop": "12 GHz", "step": "1 GHz"}})).database
    coarse = run_sweep(sweep_config(SWEEP)).database
    for cell in (0, 1):
        for i in range(3):
            assert np.array_equal(coarse.record(cell, i).s, fine.record(cell, 2 * i).s)
    for f_index, f in reversed(list(enumerate(fine.frequencies))):
        one = run_sweep(sweep_config({**SWEEP, "frequency": {"start": f"{f / 1e9} GHz", "stop": f"{f / 1e9} GHz", "step": "1 GHz"}})).database
        assert np.array_equal(one.record(0, 0).s, fine.record(0, f_index).s)


def test_lossless_cell_conserves_energy_and_is_reciprocal_at_every_swept_frequency():
    cfg = sweep_config(
        {"frequency": {"start": "6 GHz", "stop": "14 GHz", "step": "1 GHz"},
         "parameters": {"length": {"values": ["3 mm", "5 mm"]}, "rotation_deg": [0, 20]}},
        **{"substrate.tan_delta": 0.0},
    )
    db = run_sweep(cfg).database
    assert len(db.records) == 4 * 9
    strongest_cross = 0.0
    for r in db.records:
        assert np.max(np.abs(r.s.conj().T @ r.s - np.eye(4))) < 1e-11
        assert np.max(np.abs(r.s - r.s.T)) < 1e-11
        assert r.diagnostics["passivity_min_eig"] == pytest.approx(0.0, abs=1e-11)
        assert r.diagnostics["incoming_residual"] < 1e-12
        strongest_cross = max(strongest_cross, r.magnitude("S21_yx"))
        assert 0.05 < r.magnitude("S11_xx") + r.magnitude("S21_xx") < 2.0
    assert strongest_cross > 0.01                          # the rotated cells really cross-polarize


def test_lossy_cell_is_strictly_passive_at_every_frequency():
    db = run_sweep(sweep_config(SWEEP, **{"substrate.tan_delta": 0.02})).database
    for r in db.records:
        assert r.diagnostics["passivity_min_eig"] > 1e-4


def test_neighbouring_frequencies_vary_continuously_and_the_step_controls_the_change():
    def max_step(step):
        db = run_sweep(sweep_config({"frequency": {"start": "9 GHz", "stop": "11 GHz", "step": step}})).database
        s = np.array([r.component("S21_xx") for r in db.records_of(0)])
        return np.max(np.abs(np.diff(s)))

    fine, coarse = max_step("0.1 GHz"), max_step("0.4 GHz")
    assert fine < 0.05
    assert coarse / fine > 2.5                             # smooth: change per step ~ proportional to step


def test_grating_onset_is_flagged_per_record_without_warnings():
    f_c = C0 / 8e-3
    cfg = sweep_config({"frequency": {"start": "30 GHz", "stop": "40 GHz", "step": "5 GHz"}})
    assert 35e9 < f_c < 40e9
    with warnings.catch_warnings():
        warnings.simplefilter("error")                      # any UserWarning would fail the test
        db = run_sweep(cfg).database
    flags = [r.diagnostics["specular_only"] for r in db.records]
    assert flags == [True, True, False]
    assert [r.diagnostics["checks_passed"] for r in db.records] == [True, True, False]
    assert len(db.failed_validation()) == 1


# ==================================================================== phase unwrapping

PLANES_100MM = {"reference_planes.top_air_thickness": 0.1, "reference_planes.bottom_air_thickness": 0.1}


def s21_phase_sweep(step):
    sw = {"frequency": {"start": "8 GHz", "stop": "10 GHz", "step": step}}
    db = run_sweep(sweep_config(sw, **PLANES_100MM)).database
    return db, db.phase_series_deg(0, "S21_xx"), db.phase_series_deg(0, "S21_xx", unwrap=True)


def test_unwrapped_phase_of_a_real_sweep_is_independent_of_the_sampling_when_it_resolves_the_phase():
    _, (f_fine, p_fine), (_, u_fine) = s21_phase_sweep("0.125 GHz")
    _, (f_coarse, _), (_, u_coarse) = s21_phase_sweep("0.25 GHz")
    assert max_wrapped_step_deg(p_fine) < 60.0
    assert np.allclose(u_coarse, u_fine[::2], atol=1e-8)          # shared frequencies, same continuous phase
    assert np.allclose(f_coarse, f_fine[::2])


def test_unwrapped_phase_follows_the_reference_plane_delay_under_exp_plus_jwt():
    """S21 carries exp(-j k0 (z_top - z_bottom)): the unwrapped phase must fall by
    360 deg * df * d / c (d = plane separation); the cell itself only adds a small dispersion."""
    db, _, (f, u) = s21_phase_sweep("0.125 GHz")
    d = 0.1 + 1.575e-3 + 0.1
    expected = -360.0 * (f - f[0]) * d / C0
    assert u[-1] - u[0] == pytest.approx(expected[-1], abs=8.0)
    assert np.all(np.diff(u) < 0.0)                                # a delay: phase decreases monotonically
    assert u[0] == wrap_deg(u[0])                                  # anchored on the first principal value


def test_undersampled_phase_unwraps_to_the_wrong_answer_and_is_flagged():
    _, _, (_, u_fine) = s21_phase_sweep("0.125 GHz")
    with pytest.warns(UserWarning, match="aliased"):
        db, (_, p_coarse), (_, u_coarse) = s21_phase_sweep("1 GHz")
    assert max_wrapped_step_deg(p_coarse) > 90.0
    assert abs(u_coarse[1] - u_fine[8]) > 300.0                    # off by a turn (or worse) at 9 GHz
    assert np.allclose(np.mod(u_coarse - p_coarse + 180, 360) - 180, 0.0, atol=1e-9)   # still a valid unwrapping of the samples


def test_unwrapping_the_stored_phase_equals_unwrapping_the_complex_values():
    db, _, (_, u) = s21_phase_sweep("0.25 GHz")
    _, s = db.series(0, "S21_xx")
    assert np.allclose(u, unwrap_deg(np.degrees(np.angle(s))), atol=1e-8)


# ============================================================================ determinism

def read_outputs(path):
    return {n: (Path(path) / n).read_bytes() for n in ("database.json", "database.csv", CHECKPOINT_NAME)}


def test_repeated_runs_give_byte_identical_files_and_equal_databases(tmp_path):
    cfg = sweep_config({**SWEEP, "validation": {"refinement": {"enabled": True, "tolerance": 0.05}}})
    runs = [run_sweep(cfg, tmp_path / f"run{i}") for i in range(3)]
    first = read_outputs(tmp_path / "run0")
    for i in (1, 2):
        assert read_outputs(tmp_path / f"run{i}") == first
    assert runs[0].database == runs[1].database == runs[2].database
    assert run_sweep(cfg).database == runs[0].database            # in memory equals on disk


def test_files_contain_no_timestamps_paths_or_host_information(tmp_path):
    run_sweep(sweep_config(SWEEP), tmp_path)
    text = (tmp_path / "database.json").read_text()
    assert str(tmp_path) not in text and str(ROOT) not in text
    assert set(json.loads(text)) == {"format", "version", "metadata", "environment", "cells"}
    assert set(json.loads(text)["environment"]) == {"python", "numpy"}


def run_in_subprocess(cfg_path, out_dir, hash_seed, *, kill_after=None):
    code = (
        "import os, sys\n"
        "from src.config import load_config\n"
        "from src.runner import run_sweep\n"
        f"kill = {kill_after!r}\n"
        "def progress(e):\n"
        "    if kill is not None and e['done'] >= kill:\n"
        "        os._exit(3)\n"
        f"run_sweep(load_config({str(cfg_path)!r}), {str(out_dir)!r}, progress=progress)\n"
    )
    env = {**os.environ, "PYTHONHASHSEED": str(hash_seed), "PYTHONPATH": str(ROOT)}
    return subprocess.run([sys.executable, "-c", code], cwd=ROOT, env=env, capture_output=True, text=True)


def write_config(tmp_path, sweep=SWEEP):
    p = tmp_path / "cfg.yaml"
    p.write_text(yaml.safe_dump(config_dict(sweep)))
    return p


def test_output_does_not_depend_on_the_python_hash_seed(tmp_path):
    cfg = write_config(tmp_path)
    for seed in (1, 12345):
        proc = run_in_subprocess(cfg, tmp_path / f"seed{seed}", seed)
        assert proc.returncode == 0, proc.stderr
    assert read_outputs(tmp_path / "seed1") == read_outputs(tmp_path / "seed12345")
    run_sweep(sweep_config(SWEEP), tmp_path / "inproc")
    assert (tmp_path / "inproc" / "database.json").read_bytes() == (tmp_path / "seed1" / "database.json").read_bytes()


# ================================================================== resume and restart

class Stop(BaseException):
    """Simulated interruption (like KeyboardInterrupt, it is not an Exception)."""


def stop_after(n):
    def progress(event):
        if event["done"] >= n:
            raise Stop
    return progress


def checkpoint_lines(path):
    return (Path(path) / CHECKPOINT_NAME).read_text().splitlines()


def test_interrupted_sweep_resumes_and_matches_an_uninterrupted_run_exactly(tmp_path, monkeypatch):
    cfg = sweep_config(SWEEP)
    reference = run_sweep(cfg, tmp_path / "ref")
    total = N_CELLS * N_FREQ

    with pytest.raises(Stop):
        run_sweep(cfg, tmp_path / "run", progress=stop_after(4))
    assert len(checkpoint_lines(tmp_path / "run")) == 1 + 4      # header + the 4 finished results
    assert not (tmp_path / "run" / "database.json").exists()      # nothing exported from a partial run

    counter = SolveCounter(monkeypatch)
    resumed = run_sweep(cfg, tmp_path / "run")
    assert (resumed.computed, resumed.reused, resumed.dropped) == (total - 4, 4, 0)
    assert counter.n == total - 4                                  # only the missing results were solved
    assert resumed.database == reference.database
    assert read_outputs(tmp_path / "run")["database.json"] == read_outputs(tmp_path / "ref")["database.json"]
    assert read_outputs(tmp_path / "run")["database.csv"] == read_outputs(tmp_path / "ref")["database.csv"]


def test_finished_sweep_reruns_without_any_solve(tmp_path, monkeypatch):
    cfg = sweep_config({**SWEEP, "validation": {"refinement": {"enabled": True, "tolerance": 0.05}}})
    first = run_sweep(cfg, tmp_path)
    assert first.computed == N_CELLS * N_FREQ + N_CELLS and first.reused == 0
    counter = SolveCounter(monkeypatch)
    again = run_sweep(cfg, tmp_path)
    assert (again.computed, again.reused) == (0, N_CELLS * N_FREQ + N_CELLS) and counter.n == 0
    assert again.database == first.database


def test_interruption_during_a_second_cell_keeps_the_first_cells_refinement(tmp_path):
    cfg = sweep_config({**SWEEP, "validation": {"refinement": {"enabled": True, "tolerance": 0.05}}})
    with pytest.raises(Stop):
        run_sweep(cfg, tmp_path, progress=stop_after(N_FREQ + 1))      # into the second cell
    kinds = [json.loads(l)["type"] for l in checkpoint_lines(tmp_path)[1:]]
    assert kinds == ["record"] * N_FREQ + ["refinement"] + ["record"]
    done = run_sweep(cfg, tmp_path)
    assert done.reused == N_FREQ + 1 + 1 and done.computed == (N_FREQ - 1) + 1


def test_hard_kill_of_the_process_loses_nothing_that_was_reported_done(tmp_path):
    cfg_path = write_config(tmp_path)
    killed = run_in_subprocess(cfg_path, tmp_path / "out", 0, kill_after=4)
    assert killed.returncode == 3                                    # os._exit: no cleanup ran
    assert len(checkpoint_lines(tmp_path / "out")) >= 1 + 4
    resumed = run_sweep(sweep_config(SWEEP), tmp_path / "out")
    assert resumed.reused >= 4 and resumed.reused + resumed.computed == N_CELLS * N_FREQ
    assert resumed.database == run_sweep(sweep_config(SWEEP)).database


def test_torn_last_line_is_dropped_and_recomputed(tmp_path):
    cfg = sweep_config(SWEEP)
    reference = run_sweep(cfg, tmp_path / "ref")
    with pytest.raises(Stop):
        run_sweep(cfg, tmp_path / "run", progress=stop_after(4))
    path = tmp_path / "run" / CHECKPOINT_NAME
    text = path.read_text()
    path.write_text(text[: len(text) - len(text.splitlines()[-1]) // 2 - 1])     # cut the last line in half
    resumed = run_sweep(cfg, tmp_path / "run")
    assert resumed.dropped == 1 and resumed.reused == 3 and resumed.computed == N_CELLS * N_FREQ - 3
    assert resumed.database == reference.database
    assert all(json.loads(l) for l in checkpoint_lines(tmp_path / "run"))          # compacted: every line valid


def test_corrupt_and_structurally_invalid_lines_are_dropped_not_trusted(tmp_path):
    cfg = sweep_config(SWEEP)
    reference = run_sweep(cfg, tmp_path / "ref")
    run_sweep(cfg, tmp_path / "run")
    path = tmp_path / "run" / CHECKPOINT_NAME
    lines = path.read_text().splitlines()
    bad = json.loads(lines[3])
    bad["data"]["s"]["re"] = bad["data"]["s"]["re"][:3]                # wrong shape
    lines[2] = "this is not json"
    lines[3] = json.dumps(bad)
    lines[4] = json.dumps({"type": "record", "key": lines[4][:5]})       # missing data
    path.write_text("\n".join(lines) + "\n")
    resumed = run_sweep(cfg, tmp_path / "run")
    assert resumed.dropped == 3 and resumed.computed == 3 and resumed.reused == N_CELLS * N_FREQ - 3
    assert resumed.database == reference.database


def test_a_file_that_is_not_a_checkpoint_is_never_overwritten(tmp_path):
    cfg = sweep_config(SWEEP)
    tmp_path.joinpath(CHECKPOINT_NAME).write_text("precious user data\n")
    for resume in (True, False):
        with pytest.raises(CheckpointError, match="not a version-1 checkpoint"):
            run_sweep(cfg, tmp_path, resume=resume)
    assert (tmp_path / CHECKPOINT_NAME).read_text() == "precious user data\n"
    tmp_path.joinpath(CHECKPOINT_NAME).write_text(json.dumps({"format": "lightweight_mom.checkpoint", "version": 99}) + "\n")
    with pytest.raises(CheckpointError):
        run_sweep(cfg, tmp_path)


def test_empty_checkpoint_file_is_treated_as_new(tmp_path):
    tmp_path.joinpath(CHECKPOINT_NAME).write_text("")
    assert run_sweep(sweep_config(SWEEP), tmp_path).computed == N_CELLS * N_FREQ


def test_restart_discards_the_checkpoint(tmp_path, monkeypatch):
    cfg = sweep_config(SWEEP)
    first = run_sweep(cfg, tmp_path)
    counter = SolveCounter(monkeypatch)
    again = run_sweep(cfg, tmp_path, resume=False)
    assert again.reused == 0 and again.computed == N_CELLS * N_FREQ and counter.n == N_CELLS * N_FREQ
    assert again.database == first.database
    assert len(checkpoint_lines(tmp_path)) == 1 + N_CELLS * N_FREQ


def test_extending_a_sweep_computes_only_the_new_results(tmp_path, monkeypatch):
    small = {"frequency": {"start": "8 GHz", "stop": "10 GHz", "step": "2 GHz"}, "parameters": {"length": ["3 mm"]}}
    big = {"frequency": {"start": "8 GHz", "stop": "12 GHz", "step": "2 GHz"}, "parameters": {"length": ["3 mm", "5 mm"]}}
    run_sweep(sweep_config(small), tmp_path)
    counter = SolveCounter(monkeypatch)
    extended = run_sweep(sweep_config(big), tmp_path)
    assert extended.reused == 2 and extended.computed == 6 - 2 == counter.n
    assert extended.database == run_sweep(sweep_config(big)).database


def test_anything_that_changes_a_result_invalidates_the_cached_result(tmp_path):
    run_sweep(sweep_config(SWEEP), tmp_path)
    for change in ({"mesh.cells_x": 8, "mesh.cells_y": 8, "mesh.floquet_order": 16}, {"mesh.floquet_order": 33},
                   {"substrate.eps_r": 2.3}, {"substrate.thickness": 1.6e-3},
                   {"reference_planes.top_air_thickness": 1e-3}, {"geometry.width": 1e-3}):
        r = run_sweep(sweep_config(SWEEP, **change), tmp_path)
        assert r.reused == 0 and r.computed == N_CELLS * N_FREQ, change


def test_a_priori_thresholds_are_part_of_the_key_but_the_refinement_tolerance_is_not(tmp_path):
    ref = {"refinement": {"enabled": True, "tolerance": 0.05}}
    run_sweep(sweep_config({**SWEEP, "validation": ref}), tmp_path)
    same_solves = run_sweep(sweep_config({**SWEEP, "validation": {"refinement": {"enabled": True, "tolerance": 1e-9}}}), tmp_path)
    assert same_solves.computed == 0                              # only the pass/fail verdict changes
    assert all(c.validation["refinement"]["passed"] is False for c in same_solves.database.cells)
    stricter = run_sweep(sweep_config({**SWEEP, "validation": {"min_cells_per_wavelength": 60, **ref}}), tmp_path)
    assert stricter.reused == N_CELLS and stricter.computed == N_CELLS * N_FREQ     # refinements reused, records not
    assert all(r.diagnostics["checks_passed"] is False for r in stricter.database.records if r.frequency_index == 0)


def test_progress_callback_reports_every_result_in_order(tmp_path):
    events = []
    run_sweep(sweep_config(SWEEP), tmp_path, progress=events.append)
    assert [e["done"] for e in events] == list(range(1, N_CELLS * N_FREQ + 1))
    assert {e["total"] for e in events} == {N_CELLS * N_FREQ} and not any(e["reused"] for e in events)
    events.clear()
    run_sweep(sweep_config(SWEEP), tmp_path, progress=events.append)
    assert all(e["reused"] for e in events)
    assert [(e["cell_index"], e["frequency_index"]) for e in events] == [(c, f) for c in range(2) for f in range(3)]


def test_skipped_cells_are_neither_solved_nor_checkpointed(tmp_path, monkeypatch):
    sw = {**SWEEP, "parameters": {"length": ["3 mm", "9 mm"]}, "on_invalid": "skip"}
    counter = SolveCounter(monkeypatch)
    r = run_sweep(sweep_config(sw), tmp_path)
    assert counter.n == N_FREQ and r.computed == N_FREQ
    assert [c.status for c in r.database.cells] == ["ok", "invalid"] and len(checkpoint_lines(tmp_path)) == 1 + N_FREQ
    assert Database.from_json(r.json_path) == r.database


def test_out_dir_is_created_and_in_memory_runs_write_nothing(tmp_path):
    target = tmp_path / "a" / "b"
    r = run_sweep(sweep_config(SWEEP), target)
    assert r.checkpoint_path == target / CHECKPOINT_NAME and r.json_path.exists() and r.csv_path.exists()
    mem = run_sweep(sweep_config(SWEEP))
    assert mem.checkpoint_path is None and mem.json_path is None
    assert not list(target.glob("*.tmp"))                        # atomic writes leave no temp files


# ============================================================ refinement recorded per cell

def test_refinement_result_is_recorded_per_cell_with_the_documented_fields():
    sw = {**SWEEP, "validation": {"refinement": {"enabled": True, "factor": 2, "tolerance": 0.05, "frequencies": "all"}}}
    db = run_sweep(sweep_config(sw)).database
    for c in db.cells:
        ref = c.validation["refinement"]
        assert ref["frequencies_hz"] == [8e9, 10e9, 12e9] and ref["factor"] == 2 and ref["tolerance"] == 0.05
        assert 0.0 < ref["floquet_error"] < 0.05 and 0.0 < ref["mesh_error"] < 0.05 and ref["passed"] is True
        assert c.validation["area_ok"] is True and c.validation["metal_area_error"] < 1e-12
    only_max = run_sweep(sweep_config({**SWEEP, "validation": {"refinement": {"enabled": True, "tolerance": 0.05}}})).database
    assert only_max.cells[0].validation["refinement"]["frequencies_hz"] == [12e9]
    assert only_max.cells[0].validation["refinement"]["mesh_error"] <= db.cells[0].validation["refinement"]["mesh_error"] + 1e-15
    off = run_sweep(sweep_config(SWEEP)).database
    assert off.cells[0].validation["refinement"] is None


def test_a_deliberately_poor_discretization_fails_the_refinement_tolerance():
    cfg = sweep_config({**SWEEP, "validation": {"refinement": {"enabled": True, "tolerance": 0.002}}},
                       **{"mesh.cells_x": 8, "mesh.cells_y": 8, "mesh.floquet_order": 8, "geometry.length": 4e-3})
    db = run_sweep(cfg, definition=dataclasses.replace(cfg.sweep, axes=())).database
    ref = db.cells[0].validation["refinement"]
    assert ref["passed"] is False and max(ref["floquet_error"], ref["mesh_error"]) > 0.002
    assert [c.cell_id for c in db.failed_validation()] == [db.cells[0].cell_id]


# ============================================================ physical checks of stored data

def rect_rotation_db(**extra):
    sw = {"frequency": {"start": "10 GHz", "stop": "10 GHz", "step": "1 GHz"},
          "parameters": {"rotation_deg": [0, 90]}}
    return run_sweep(sweep_config(sw, **extra)).database


def test_rotating_the_cell_by_90_degrees_swaps_x_and_y_exactly():
    db = rect_rotation_db()
    s0, s90 = db.record(0, 0).s, db.record(1, 0).s
    r = rot(90.0)
    for i in (0, 1):
        for j in (0, 1):
            blk0 = s0[2 * i : 2 * i + 2, 2 * j : 2 * j + 2]
            blk90 = s90[2 * i : 2 * i + 2, 2 * j : 2 * j + 2]
            assert np.allclose(blk90, r @ blk0 @ r.T, atol=1e-10)
    assert abs(db.record(0, 0).component("S21_xx") - db.record(0, 0).component("S21_yy")) > 0.05     # non-vacuous


def test_geometric_phase_of_a_90_degree_rotation_in_the_circular_basis():
    """RR, LL unchanged; LR' = exp(+2j psi) LR and RL' = exp(-2j psi) RL with psi = 90 deg (factor -1)."""
    db = rect_rotation_db()
    a, b = db.record(0, 0), db.record(1, 0)
    for fam in ("S11", "S12", "S21", "S22"):
        assert b.component(f"{fam}_RR") == pytest.approx(a.component(f"{fam}_RR"), abs=1e-10)
        assert b.component(f"{fam}_LL") == pytest.approx(a.component(f"{fam}_LL"), abs=1e-10)
        assert b.component(f"{fam}_LR") == pytest.approx(-a.component(f"{fam}_LR"), abs=1e-10)
        assert b.component(f"{fam}_RL") == pytest.approx(-a.component(f"{fam}_RL"), abs=1e-10)
    assert abs(a.component("S21_LR")) > 0.05 and abs(a.component("S21_RL")) > 0.05


def rotated_db():
    sw = {"frequency": {"start": "9 GHz", "stop": "11 GHz", "step": "1 GHz"},
          "parameters": {"rotation_deg": [30, 60]}}
    return run_sweep(sweep_config(sw, **{"substrate.tan_delta": 0.0, "geometry.length": 4e-3})).database


def test_circular_basis_reciprocity_relations_on_cross_polarizing_cells():
    """Linear S = S^T becomes S12_RR = S21_LL, S12_LL = S21_RR, S12_LR = S21_LR, S12_RL = S21_RL
    (conjugation swaps R and L, derivation: DATABASE.md)."""
    db = rotated_db()
    for r in db.records:
        assert r.component("S12_RR") == pytest.approx(r.component("S21_LL"), abs=1e-10)
        assert r.component("S12_LL") == pytest.approx(r.component("S21_RR"), abs=1e-10)
        assert r.component("S12_LR") == pytest.approx(r.component("S21_LR"), abs=1e-10)
        assert r.component("S12_RL") == pytest.approx(r.component("S21_RL"), abs=1e-10)
        assert abs(r.component("S21_LR")) > 1e-3 and abs(r.component("S21_yx")) > 1e-3


def test_circular_matrix_of_a_lossless_cell_is_unitary_and_power_is_conserved_per_polarization():
    for r in rotated_db().records:
        c = matrix4_to_circular(r.s)
        assert np.max(np.abs(c.conj().T @ c - np.eye(4))) < 1e-10
        for inc in ("R", "L"):                               # every circularly polarized input, port 1
            p = sum(abs(r.component(f"{fam}_{rec}{inc}")) ** 2 for fam in ("S11", "S21") for rec in "RL")
            assert p == pytest.approx(1.0, abs=1e-10)


def helicity_of_field(e_xy, direction):
    """RHCP / LHCP by the right-hand rule on the physical field Re{e exp(j w t)} (independent of src)."""
    t = np.linspace(0.0, 1.0, 2001)[:-1]
    field = np.real(np.outer(np.exp(2j * np.pi * t), e_xy))
    torque = np.mean(field[:, 0] * np.gradient(field[:, 1], t) - field[:, 1] * np.gradient(field[:, 0], t))
    about_propagation = torque if direction == "+z" else -torque
    return "RHCP" if about_propagation > 0 else "LHCP"


def test_pec_sheet_reflects_the_R_label_as_R_label_but_reverses_the_physical_handedness():
    sw = {"frequency": {"start": "8 GHz", "stop": "12 GHz", "step": "2 GHz"}}
    db = run_sweep(sweep_config(sw, **{"geometry.length": 8e-3, "geometry.width": 8e-3})).database
    for r in db.records:
        assert r.component("S11_RR") == pytest.approx(-1.0, abs=1e-8)
        assert r.component("S11_LL") == pytest.approx(-1.0, abs=1e-8)
        assert abs(r.component("S11_LR")) < 1e-8 and abs(r.component("S11_RL")) < 1e-8
        assert abs(r.component("S21_RR")) < 1e-8 and abs(r.component("S21_LR")) < 1e-8
        incident = E_R                                          # R-labelled wave, travelling -z
        reflected = r.s[:2, :2] @ incident                       # what the solver says comes back
        assert helicity_of_field(incident, "-z") == "RHCP"
        assert helicity_of_field(reflected, "+z") == "LHCP"     # same label, opposite handedness
        assert np.allclose(reflected, -incident, atol=1e-8)


def test_transmission_through_the_cell_keeps_the_handedness_of_an_incident_circular_wave():
    r = run_sweep(sweep_config({})).database.record(0, 0)
    incident = E_R
    transmitted = r.s[2:, :2] @ incident                        # port 2 output, travelling -z
    assert helicity_of_field(incident, "-z") == helicity_of_field(transmitted, "-z") == "RHCP"
    assert abs(r.component("S21_RR")) > 10 * abs(r.component("S21_LR"))         # dominantly the same label


def test_thinner_strips_approach_the_bare_dielectric_slab():
    """Vanishing metal limit on stored data: the distance to the Fresnel/Airy slab result grows
    monotonically with the strip width."""
    cfg = sweep_config({"parameters": {"width": ["1 mm", "2 mm", "3 mm", "4 mm"]}})
    db = run_sweep(cfg).database
    green = LayeredGreenFunction(cfg.unit_cell.substrate, 10e9)
    slab = stack_s_parameters(green, cfg.reference_planes).matrix()
    dist = [np.max(np.abs(db.record(i, 0).s - slab)) for i in range(4)]
    assert dist == sorted(dist) and dist[0] < 0.6 * dist[-1]


def test_swept_geometry_parameters_are_stored_exactly_and_deterministically():
    sw = {"parameters": {"length": {"start": "2 mm", "stop": "6 mm", "step": "1 mm"}, "center_x": ["0 mm", "1 mm"]}}
    a, b = run_sweep(sweep_config(sw)).database, run_sweep(sweep_config(sw)).database
    assert a == b and len(a.cells) == 10
    assert a.parameter_array("length_m").tolist() == pytest.approx([2e-3] * 2 + [3e-3] * 2 + [4e-3] * 2 + [5e-3] * 2 + [6e-3] * 2)
    assert a.parameter_array("center_x_m").tolist() == [0.0, 1e-3] * 5
    assert len({c.cell_id for c in a.cells}) == 10


# ============================================================ rectangular (non-square) lattice

def test_non_square_lattice_is_stored_and_rebuilt_with_the_periods_in_the_right_order(tmp_path):
    """Every other test uses a square cell, which cannot notice x/y period mix-ups."""
    overrides = {"period_x": 8e-3, "period_y": 6e-3, "mesh.cells_x": 16, "mesh.cells_y": 12}
    sw = {"frequency": {"start": "9 GHz", "stop": "11 GHz", "step": "2 GHz"}, "parameters": {"length": ["3 mm", "4 mm"]}}
    cfg = sweep_config(sw, **overrides)
    r = run_sweep(cfg, tmp_path)
    db = Database.from_json(r.json_path)
    assert db.metadata["context"]["lattice"] == {"period_x_m": 8e-3, "period_y_m": 6e-3}
    for cell in db.cells:
        uc = db.unit_cell(cell.index)
        assert (uc.lattice.period_x, uc.lattice.period_y) == (8e-3, 6e-3)
        length = cell.parameters["length_m"]
        for rec in db.records_of(cell.index):
            fresh = compute_s_parameters(with_length(cfg, length), rec.frequency, cfg.mesh, cfg.reference_planes, warn_grating=False)
            assert np.max(np.abs(fresh.matrix() - rec.s)) < 1e-13
    swapped = sweep_config(sw, **{"period_x": 6e-3, "period_y": 8e-3, "mesh.cells_x": 12, "mesh.cells_y": 16})
    other = run_sweep(swapped).database
    assert np.max(np.abs(other.record(0, 0).s - db.record(0, 0).s)) > 1e-3        # x and y really differ
    assert other.cells[0].cell_id != db.cells[0].cell_id


def test_grating_onset_of_a_non_square_cell_is_set_by_the_larger_period():
    sw = {"frequency": {"start": "30 GHz", "stop": "45 GHz", "step": "5 GHz"}}
    db = run_sweep(sweep_config(sw, **{"period_x": 8e-3, "period_y": 6e-3, "mesh.cells_x": 16, "mesh.cells_y": 12})).database
    assert C0 / 8e-3 < 40e9 and 35e9 < C0 / 8e-3                                 # onset 37.5 GHz (x), 50 GHz (y)
    assert [r.diagnostics["specular_only"] for r in db.records] == [True, True, False, False]


def test_refinement_factor_is_part_of_the_refinement_key_but_not_of_the_record_key(tmp_path):
    ref = lambda factor: {"validation": {"refinement": {"enabled": True, "factor": factor, "tolerance": 0.05}}}
    first = run_sweep(sweep_config({**SWEEP, **ref(2)}), tmp_path)
    again = run_sweep(sweep_config({**SWEEP, **ref(3)}), tmp_path)
    assert again.reused == N_CELLS * N_FREQ and again.computed == N_CELLS                 # only the studies redone
    a, b = first.database.cells[0].validation["refinement"], again.database.cells[0].validation["refinement"]
    assert (a["factor"], b["factor"]) == (2, 3) and a["mesh_error"] != b["mesh_error"]
