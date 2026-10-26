"""Phase 27 item 1: the tensor-CP design module reproduces the committed Phase 26 layouts and predictions exactly."""

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from src.hologram import tensor_cp_design as t
from src.hologram.tensor_leaky_cp import SupercellTable

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def table():
    return SupercellTable.from_files(t.TABLE_GLOB)


def test_single_hexagon_cells_match_committed_layout_bit_identically(table):
    cells = t.single_hexagon_cells(table)
    rows = list(csv.reader((ROOT / "hfss_exports/hex150_tensor_cp_single/layout.csv").read_text().splitlines()))[1:]
    assert len(cells) == len(rows)
    for c, r in zip(cells, rows, strict=True):
        assert (
            f"{c['x_mm']:.3f}",
            f"{c['y_mm']:.3f}",
            f"{c['slot_L_mm']:.2f}",
            f"{c['slot_W_mm']:.2f}",
            f"{c['slot_angle_deg']:.3f}",
        ) == (r[0], r[1], r[3], r[4], r[5])
    # and the full-precision cell list of the committed design JSON
    ref = json.loads((ROOT / "results/phase26/tensor_cp_hexagon_cells.json").read_text())
    assert [c["x_mm"] for c in cells] == ref["x_mm"] and [c["slot_angle_deg"] for c in cells] == ref["slot_angle_deg"]


def test_fourhex_boards_match_committed_layouts(table):
    for board, csv_path in (("A", ROOT / "hfss_exports/hex4_tensor_cp_static/layout.csv"),):
        bd = t.design_board(board, table)
        rows = list(csv.reader(csv_path.read_text().splitlines()))[1:]
        k = 0
        for j in t.FEEDS:
            d = bd[j]["design"]
            for i in range(bd[j]["x_mm"].size):
                r = rows[k]
                assert (
                    f"{bd[j]['x_mm'][i]:.3f}",
                    f"{bd[j]['y_mm'][i]:.3f}",
                    f"{d['slot_L_mm'][i]:.2f}",
                    f"{d['slot_W_mm'][i]:.2f}",
                    f"{np.degrees(d['slot_angle_rad'][i]):.3f}",
                    str(j),
                ) == (r[0], r[1], r[3], r[4], r[5], r[6])
                k += 1
        assert k == len(rows)


@pytest.mark.slow
def test_full_single_scan_and_fourhex_results_reproduce_committed_json():
    res, cells = t.design_single_hexagon(log=lambda *a, **k: None)
    ref = json.loads((ROOT / "results/phase26/tensor_cp_hexagon.json").read_text())
    res = json.loads(json.dumps(res, default=float))
    res.pop("seconds"), ref.pop("seconds")
    assert res == ref
    assert json.loads(json.dumps(cells)) == json.loads(
        (ROOT / "results/phase26/tensor_cp_hexagon_cells.json").read_text()
    )
    r4, c4, _ = t.design_fourhex(log=lambda *a, **k: None)
    ref4 = json.loads((ROOT / "results/phase26/fourhex/fourhex_design.json").read_text())
    r4 = json.loads(json.dumps(r4, default=float))
    r4.pop("seconds"), ref4.pop("seconds")
    assert r4 == ref4
    for b in ("A", "B"):
        assert json.loads(json.dumps({"board": b, "cells": c4[b]})) == json.loads(
            (ROOT / f"results/phase26/fourhex/fourhex_{b}_cells.json").read_text()
        )


def test_prediction_builders_reproduce_committed_predictions_exactly():
    res = json.loads((ROOT / "results/phase26/tensor_cp_hexagon.json").read_text())
    got = json.loads(json.dumps(t.single_hex_predictions(res), default=float))
    assert got == json.loads((ROOT / "hfss_exports/hex150_tensor_cp_single/predictions.json").read_text())
    design = json.loads((ROOT / "results/phase26/fourhex/fourhex_design.json").read_text())
    scope = json.loads((ROOT / "results/phase26/fourhex_scope.json").read_text())
    n = len(json.loads((ROOT / "results/phase26/fourhex/fourhex_A_cells.json").read_text())["cells"])
    got4 = json.loads(json.dumps(t.fourhex_a_predictions(design, scope, n), default=float))
    assert got4 == json.loads((ROOT / "hfss_exports/hex4_tensor_cp_static/predictions.json").read_text())
