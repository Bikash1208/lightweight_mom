"""The band-refinement job driver: stages, task lists and order, whole-job deadline, checkpoint names, consolidation of finished tasks; no solves."""

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import phase20_r4_band_refine_run as br  # noqa: E402


def _reset():
    br.lib.configure_grid(f_step=0.25, theta_step=7.5, phi_step=30.0, corr_f_list=[8.0 + 0.5 * i for i in range(9)])


def test_stage_task_counts_and_order():
    try:
        g, c90, c120 = (br.build_tasks(s) for s in br.STAGES)
        assert len(g) == 5 * 217 and {t[3] for t in g} == {60} and {t[0] for t in g} == set(br.BAND_F)
        assert len(c90) == 50 and len(c120) == 50 and {t[3] for t in c90} == {90} and {t[3] for t in c120} == {120}
        assert {t[0] for t in c90} == set(br.CORR_NODES) and len(br.CORR_NODES) == 10
        assert br.build_tasks() == g + c90 + c120  # job order: grid, all mesh 90, all mesh 120
        # the nodes are exactly the 0.05 GHz steps of 10.4-10.9 that carry no correction yet (10.50 has one from R4)
        steps = {round(10.4 + 0.05 * k, 2) for k in range(11)}
        assert set(br.CORR_NODES) == steps - {10.5} and set(br.BAND_F) == {
            f for f in steps if round(f * 10) != f * 10 and abs(f * 10 - round(f * 10)) > 1e-9
        }
    finally:
        _reset()


def test_job_deadline_is_for_the_whole_job(tmp_path):
    d1 = br.job_deadline(tmp_path, 8.0)
    assert (tmp_path / "job_start.txt").exists()
    t = float((tmp_path / "job_start.txt").read_text())
    (tmp_path / "job_start.txt").write_text(str(t - 3600.0))  # an hour has passed since the first run
    d2 = br.job_deadline(tmp_path, 8.0)
    assert abs((d1 - d2) - 3600.0) < 5.0 and abs(d2 - time.monotonic() - 7.0 * 3600) < 5


def test_consolidate_builds_library_layout_chunks_from_task_files(tmp_path):
    try:
        # fabricate the 5 correction-point task files of one correction node at mesh 90 and 120
        for mesh in (90, 120):
            for f, th, ph, m in [t for t in br.build_tasks({90: "corr90", 120: "corr120"}[mesh]) if t[0] == 10.4]:
                (tmp_path / "tasks").mkdir(exist_ok=True)
                np.save(tmp_path / "tasks" / br.fg.task_name(f, th, ph, m), np.full((4, 2, 2), mesh, dtype=complex))
        br.consolidate(tmp_path)
        z = np.load(tmp_path / "chunks" / "corr120_f10.40.npz")
        assert z["s"].shape == (5, 4, 2, 2) and int(z["mesh"]) == 120 and (z["s"] == 120).all()
        assert (tmp_path / "chunks" / "corr90_f10.40.npz").exists() and not (
            tmp_path / "chunks" / "corr90_f10.45.npz"
        ).exists()
        import json

        m = json.loads((tmp_path / "MANIFEST.json").read_text())
        assert m["chunks_complete"]["corr90"] == [10.4] and m["tasks_defined"] == 5 * 217 + 100
    finally:
        _reset()
