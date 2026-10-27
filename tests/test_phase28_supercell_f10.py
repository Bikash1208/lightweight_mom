"""Phase 28: the 10.0 GHz supercell-table driver -- frequency override, task list and cost model, assembly into the Phase 26 table format (no MoM is run here)."""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))
sys.path.insert(0, str(ROOT))

import phase28_supercell_f10 as drv  # noqa: E402

from src.frequency import free_space_wavenumber  # noqa: E402
from src.hologram.tensor_leaky_cp import SupercellTable  # noqa: E402


def test_set_frequency_points_both_phase26_scripts_at_the_new_frequency_and_back():
    st = drv.set_frequency(10.0e9)
    import phase26_ramp_locality as rl

    assert rl.F == st.F == 10.0e9 and rl.K0 == st.K0 == free_space_wavenumber(10.0e9)
    assert abs(st.LAMBDA0 - 299792458.0 / 10.0e9) < 1e-12
    drv.set_frequency(10.54e9)  # restore: other tests import the same modules
    assert rl.F == st.F == 10.54e9


def test_periods_cover_the_committed_set_and_extend_both_ways():
    committed = {3, 4, 5, 6, 8}  # K / k0 = lambda0 / (N P) of results/phase26/supercell_table_*.json
    assert committed == set(drv.TIER1)
    assert min(drv.TIER2) < min(drv.TIER1) and max(drv.TIER2) > max(drv.TIER1)
    assert len(drv.tasks_for("all")) == len(drv.SLOTS) * (len(drv.TIER1) + len(drv.TIER2))
    assert {s for s, _ in drv.tasks_for("1")} == set(drv.SLOTS)


def test_cost_and_memory_models_grow_with_the_period_count():
    ns = sorted(drv.TIER1 + drv.TIER2)
    assert all(drv.est_seconds(a) < drv.est_seconds(b) for a, b in zip(ns[:-1], ns[1:], strict=True))
    assert all(drv.est_mem_gb(a) < drv.est_mem_gb(b) for a, b in zip(ns[:-1], ns[1:], strict=True))
    assert drv.est_seconds(3) == drv.COST_S_AT_N3


def test_assemble_writes_tables_the_phase26_reader_accepts(tmp_path):
    for slot in drv.SLOTS:
        for n in (3, 4, 8):
            rec = {
                "n_cells": n,
                "period_mm": n * 6.4,
                "K_over_k0": 4.685 / n,
                "slot_L_W_mm": list(slot),
                "diameter_mm": 6.0,
                "beta_real_root_over_k0": 1.05,
                "beta_over_k0": 1.05,
                "alpha_np_per_m": 1.0 + n,
                "beta_uniform_over_k0": 1.04,
                "alpha_uniform_np_per_m": 0.01,
                "alpha_leak_np_per_m": 0.5 * n,
                "task_seconds": 1.0,
                "peak_working_set_gb": 0.1,
                "frequency_ghz": 10.0,
            }
            (tmp_path / (drv.task_name(slot, n) + ".json")).write_text(json.dumps(rec), encoding="utf-8")
    drv.assemble(tmp_path)
    t = SupercellTable.from_files(str(tmp_path / "supercell_table_*.json"))
    assert t.alpha.shape == (4, 3) and np.all(np.isfinite(t.alpha)) and t.pitch_mm == 6.4 and t.diameter_mm == 6.0
    assert np.allclose(sorted(t.k_grid), sorted(4.685 / np.array([3, 4, 8])), atol=1e-3)
    assert (
        "task_seconds"
        not in json.loads(next(tmp_path.glob("supercell_table_*.json")).read_text(encoding="utf-8"))["points"][0]
    )
