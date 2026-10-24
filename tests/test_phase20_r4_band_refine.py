"""The band-refinement job driver: task list, ordering, checkpoint names; no solves."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import phase20_r4_band_refine_run as br  # noqa: E402


def test_task_counts_and_order():
    try:
        t = br.build_tasks(mesh120=True)
        assert len(t) == 25 + 25 + 5 * 217
        assert [x[3] for x in t[:25]] == [120] * 25 and [x[3] for x in t[25:50]] == [90] * 25 and {x[3] for x in t[50:]} == {60}
        assert {x[0] for x in t} == set(br.BAND_F) and not any(abs(f * 10 - round(f * 10)) < 1e-9 for f in br.BAND_F)  # none on the 0.1 GHz grid
        assert len(br.build_tasks(mesh120=False)) == 25 + 5 * 217
        mesh60 = {(x[0], x[1], x[2]) for x in t if x[3] == 60}
        assert all((f, th, ph) in mesh60 for f, th, ph, m in t if m == 90)  # every correction point is also a mesh-60 grid point
    finally:
        br.lib.configure_grid(f_step=0.25, theta_step=7.5, phi_step=30.0, corr_f_list=[8.0 + 0.5 * i for i in range(9)])
