"""Fewer-diodes engine: D2 groups, D3 one-axis rule, D6 core mask, D5 board geometry, rotated-cell option bookkeeping."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import sidequest_arch_study as st  # noqa: E402
import sidequest_coded_k4 as K  # noqa: E402
import sidequest_diode_reduction as D  # noqa: E402


def _run(**spec):
    bd = D.board("K2", spec.pop("dphi", None), spec.pop("dtheta", None))
    opt = K.options(2, D.VS_AN, D.G2)
    kw = {}
    if "gid" in spec:
        kw["gid"] = D.group_ids(bd, spec.pop("gid"))
    amp = np.sqrt(2 * K.A0 / bd.rho) * np.exp(-0.5 * (K.A0 + K.AA0) * (bd.rho - st.RHO_KEEP))
    return bd, D.design_g(bd, amp, opt, [0.0, 0.0], spec.get("lam", 0.4), 0.05, 2.0, **kw, **{k: v for k, v in spec.items() if k not in ("lam",)}), opt


def test_baseline_k2_diode_count_and_groups_never_cost_more_objective():
    _, (m, nd, loss, ids, J, best), _ = _run(lam=0.266)
    assert 700 <= nd.sum() <= 790
    _, (_, nd2, _, _, J2, _), _ = _run(lam=0.266, gid="2x2")
    assert J2 <= J + 1e-9  # a group is shared only where it lowers the objective


def test_one_axis_rule_has_no_cell_with_two_switched_elements_and_core_mask_holds():
    bd, (m, nd, loss, ids, J, best), opt = _run(lam=0.4, no_both=True)
    assert not np.any((ids[0] >= 0) & (ids[1] >= 0))
    bd, (m, nd, loss, ids, J, best), opt = _run(lam=0.2, forbid_sw=D.board("K2").rho > 0.06)
    out = D.board("K2").rho > 0.06
    assert nd[out].sum() == 0 and nd[~out].sum() > 0


def test_pointing_offsets_move_the_beam_phase_and_group_ids_partition_the_board():
    b0, b1 = D.board("K2"), D.board("K2", [2.0, -2.0], [0.0, 0.0])
    assert not np.allclose(b0.arr[0], b1.arr[0]) and not np.allclose(b0.arr[1], b1.arr[1])
    assert np.allclose(b0.arr[0] / b0.arr[0], 1.0)
    gid = D.group_ids(b0, "2x2")
    sizes = np.bincount(gid)
    assert sizes.sum() == b0.n and sizes.max() <= 4
