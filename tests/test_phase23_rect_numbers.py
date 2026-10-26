"""`scripts/research/phase23_rect_numbers.py`: the numeric version of the Phase 23 rectangle gate. The tests exercise the flatten / compare logic on synthetic results
(identical, a last-digit change, a structural change); the gate itself (a full synthesis of six layouts, ~15 s) is run as documented in docs/DEVELOPING.md."""

import copy
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))

import phase23_rect_numbers as rn  # noqa: E402

RESULT = {
    "frequency_hz": 10.54e9,
    "ports": [{"n_cells": 49, "phase": np.array([[0.1, -0.2], [0.3, 0.4]]), "label": "a", "flag": True, "none": None}],
    "patch_layout_png_url": "/api/x/abc.png",  # excluded: carries the job id
    "achieved_phase_deg": [1.0, 2.0],  # excluded: key added later
    "z": 1 + 2j,
}


def test_flatten_skips_the_two_excluded_keys_and_keeps_every_other_leaf():
    flat = rn.flatten(RESULT)
    assert "/patch_layout_png_url" not in flat and "/achieved_phase_deg#len" not in flat
    assert flat["/frequency_hz"] == 10.54e9
    assert flat["/ports[0]/phase"] == {"nd": "float64", "shape": [2, 2], "data": [0.1, -0.2, 0.3, 0.4]}
    assert flat["/ports[0]/n_cells"] == 49 and flat["/ports[0]/flag"] is True and flat["/ports[0]/none"] == "None"
    assert flat["/z"] == {"re": 1.0, "im": 2.0}


def test_identical_results_compare_identical():
    a = {"c": rn.flatten(RESULT)}
    info = rn.compare_numbers(a, copy.deepcopy(a))
    assert info["identical"] and info["n_diff"] == 0 and info["n_float"] > 0


def test_a_last_digit_change_is_reported_with_its_size_and_place():
    other = copy.deepcopy(RESULT)
    other["frequency_hz"] = np.nextafter(10.54e9, 1e12)
    other["ports"][0]["phase"][1, 0] = np.nextafter(0.3, 1.0)
    info = rn.compare_numbers({"c": rn.flatten(RESULT)}, {"c": rn.flatten(other)})
    assert not info["identical"] and info["n_diff"] == 2 and not info["structure"]
    assert 0.0 < info["max_rel"] < 1e-15 and info["where_rel"] in ("c/frequency_hz", "c/ports[0]/phase")


def test_a_changed_shape_count_or_label_is_a_structure_difference():
    base = {"c": rn.flatten(RESULT)}
    for mutate in (
        lambda r: r["ports"].append({"n_cells": 1}),
        lambda r: r["ports"][0].__setitem__("label", "b"),
        lambda r: r["ports"][0].__setitem__("phase", np.zeros((4,))),
    ):
        other = copy.deepcopy(RESULT)
        mutate(other)
        assert rn.compare_numbers(base, {"c": rn.flatten(other)})["identical"] is False


def test_nan_equals_nan_and_missing_case_is_flagged():
    a = {"c": rn.flatten({"x": float("nan")})}
    assert rn.compare_numbers(a, copy.deepcopy(a))["identical"]
    assert rn.compare_numbers(a, {"d": a["c"]})["structure"]
