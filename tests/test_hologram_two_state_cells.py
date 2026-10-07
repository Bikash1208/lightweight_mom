"""Phase 15 Part E: two-state (1-bit) cell selection for the Coded design method."""

import math

import pytest

from src.hologram.two_state_cells import (
    TwoStateCell,
    best_separated_pair,
    candidates_from_annular_slot_library,
    candidates_from_square_library,
    separation_warning,
    two_state_from_values,
)


def test_best_separated_pair_picks_the_closest_to_180deg():
    candidates = [
        TwoStateCell("a", "a", complex(1.0, 0.0), {}),  # 0deg
        TwoStateCell("b", "b", 1.0 * (0.5 + 0.866j), {}),  # 60deg
        TwoStateCell("c", "c", 1.0 * (-1.0 + 0.0j), {}),  # 180deg -- best pair with "a"
    ]
    choice = best_separated_pair(candidates)
    assert {choice.state_0.label, choice.state_1.label} == {"a", "c"}
    assert choice.achieved_separation_deg == pytest.approx(180.0, abs=1e-6)


def test_best_separated_pair_reports_an_imperfect_separation_honestly():
    candidates = [
        TwoStateCell("a", "a", complex(1.0, 0.0), {}),  # 0deg
        TwoStateCell("b", "b", 1.0 * (0.5 + 0.866j), {}),  # 60deg -- best available separation
    ]
    choice = best_separated_pair(candidates)
    assert choice.achieved_separation_deg == pytest.approx(60.0, abs=1e-3)


def test_best_separated_pair_requires_at_least_two_candidates():
    with pytest.raises(ValueError):
        best_separated_pair([TwoStateCell("a", "a", complex(1.0, 0.0), {})])


def test_candidates_from_square_library_wraps_each_entry():
    raw = [("c1", 2.0, 2.0, "lib", complex(1.0, 0.0)), ("c2", 3.0, 3.0, "lib", complex(-1.0, 0.0))]
    wrapped = candidates_from_square_library(raw)
    assert len(wrapped) == 2
    assert wrapped[0].geometry == {"kind": "rectangle", "L_mm": 2.0, "W_mm": 2.0}
    assert wrapped[1].value == complex(-1.0, 0.0)


def test_candidates_from_annular_slot_library_skips_failed_cells_and_picks_nearest_frequency():
    docs = [
        {
            "solver_status": "OK",
            "cell_id": "ok1",
            "a_mm": 4.0,
            "b_mm": 4.0,
            "R_mm": 1.5,
            "r_mm": 1.0,
            "frequency_GHz": [12.0, 12.5, 13.0],
            "S11": [{"real": 0.1, "imag": 0.1}, {"real": -1.0, "imag": 0.0}, {"real": 0.5, "imag": 0.5}],
        },
        {"solver_status": "FAILED", "cell_id": "bad", "frequency_GHz": [], "S11": []},
    ]
    wrapped = candidates_from_annular_slot_library(docs, frequency_ghz=12.5)
    assert len(wrapped) == 1
    assert wrapped[0].value == complex(-1.0, 0.0)
    assert wrapped[0].geometry == {"kind": "annular_slot", "a_mm": 4.0, "b_mm": 4.0, "R_mm": 1.5, "r_mm": 1.0}


def test_two_state_from_values_computes_separation():
    choice = two_state_from_values("full patch", complex(1.0, 0.0), "slotted patch", complex(-0.98, 0.0))
    assert choice.achieved_separation_deg == pytest.approx(180.0, abs=2.0)


def test_separation_warning_none_when_close_to_180deg():
    choice = two_state_from_values(
        "a", complex(1.0, 0.0), "b", complex(math.cos(math.radians(175)), math.sin(math.radians(175)))
    )
    assert separation_warning(choice) is None


def test_separation_warning_flags_a_poor_pair():
    choice = two_state_from_values(
        "a", complex(1.0, 0.0), "b", complex(math.cos(math.radians(60)), math.sin(math.radians(60)))
    )
    warning = separation_warning(choice)
    assert warning is not None
    assert "60.0" in warning
