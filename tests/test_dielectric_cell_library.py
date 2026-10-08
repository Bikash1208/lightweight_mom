"""Phase 19 Part B -- nylon cell library: geometry limits, ordering/window logic, recorded results."""

import json
import math
import pathlib

import numpy as np
import pytest

from src.dielectric_cell.library import (
    COS_FLAT,
    LibraryParams,
    best_window,
    cell_list,
    make_cell,
    normal_incidence_table,
    order_by_fill,
    sweep,
)

SMALL = LibraryParams(28e9, 5e-3, 25e-3, 1e-3, 2.45, 0.0185, 0.6e-3, 3, (0.0, 20.0), (0.0, 45.0), 4)
SUMMARY = pathlib.Path(__file__).resolve().parents[1] / "results" / "nylon_cells_28ghz" / "summary.json"


def test_every_cell_respects_the_minimum_feature():
    for kind, r in cell_list(SMALL):
        c = make_cell(SMALL, kind, r)  # constructor also rejects octagons that do not fit
        if r > 0:
            flat_to_flat = 2 * r * COS_FLAT
            assert flat_to_flat >= SMALL.min_feature_m * (1 - 1e-9)  # pillar width / hole width
            assert SMALL.period_m - flat_to_flat >= SMALL.min_feature_m * (1 - 1e-9)  # gap / wall
        assert c.total_thickness_m == pytest.approx(25e-3)


def test_cell_list_has_empty_and_solid_endpoints():
    cl = cell_list(SMALL)
    assert cl[0] == ("pillar", 0.0) and cl[-1] == ("hole", 0.0)
    assert len(cl) == 2 * SMALL.n_per_family + 2


def test_fill_fraction_of_the_two_families_is_complementary():
    a = make_cell(SMALL, "pillar", 2e-3).fill_fraction
    b = make_cell(SMALL, "hole", 2e-3).fill_fraction
    assert a + b == pytest.approx(1.0)
    assert a == pytest.approx(2 * math.sqrt(2) * (2e-3) ** 2 / (5e-3) ** 2)


def test_sweep_endpoints_empty_and_solid():
    rows = sweep(SMALL, workers=1)
    tab = {(k, round(r, 9)): (ph, ls) for k, r, _f, ph, ls in normal_incidence_table(rows)}
    # solid cell is the slowest-phase / lossiest, the empty cell is the fastest / least lossy
    assert tab[("hole", 0.0)][1] < tab[("pillar", 0.0)][1] < 0
    # oblique rows present for both azimuths, normal only once
    assert sum(1 for r in rows if r[0] == "pillar" and r[1] == 0.0 and r[3] == 20.0) == 2
    assert sum(1 for r in rows if r[0] == "pillar" and r[1] == 0.0 and r[3] == 0.0) == 1


def test_order_by_fill_unwraps_and_window_finds_span():
    # synthetic: phase falls 100 deg per unit step, wrapped; loss worst in the middle
    tab = [
        ("pillar", i, i / 10, ((-100.0 * i + 180) % 360) - 180, -1.0 - (0.5 if 5 <= i <= 6 else 0)) for i in range(10)
    ]
    ordered = order_by_fill(tab)
    assert ordered[-1][3] - ordered[0][3] == pytest.approx(-900.0)
    i0, i1, worst, span = best_window(ordered, 360.0)
    assert span >= 360.0 and i1 - i0 == 4
    assert best_window(ordered, 2000.0) is None


@pytest.mark.skipif(not SUMMARY.exists(), reason="library not built (examples/build_nylon_cell_library_28ghz.py)")
def test_recorded_library_meets_the_phase_span_requirement():
    s = json.loads(SUMMARY.read_text())
    assert s["full_library_phase_span_deg"] >= 360.0
    w = s["best_360_window"]
    assert w is not None and w["span_deg"] >= 360.0
    # loss is physical: bounded by the 25 mm solid-nylon absorption (~1.8 dB) plus Fresnel ripple
    assert -3.0 < w["worst_case_loss_db"] < -0.5
    assert s["hfss_comparison"] == "NOT_AVAILABLE_FROM_CURRENT_DATA"
    assert np.isclose(s["params"]["min_feature_m"], 0.6e-3)
