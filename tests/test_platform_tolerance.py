"""The platform-aware reference comparison (`tests/_platform_tol.py`) must accept last-digit noise and reject a real change (1e-9 relative)."""

import copy

import pytest

from tests._platform_tol import (
    CROSS_PLATFORM_RTOL,
    assert_matches_reference,
    assert_text_matches_reference,
)

REF = {"name": "cell", "n": 7, "ok": True, "x": [30.757828525351485, -1.25e-3, 0.0], "s": {"s11": [0.4, 1.6e-16]}}
TEXT = 'FEEDS = [(64.952, 112.500)]\nGAIN = 23.030000000000001\n"n_cells": 1430,\n'


def _perturbed(rel):
    other = copy.deepcopy(REF)
    other["x"][0] *= 1 + rel
    return other


def test_the_tolerance_is_the_documented_one():
    assert CROSS_PLATFORM_RTOL == 1e-12


def test_tolerance_path_accepts_last_digit_noise_and_noise_level_entries():
    other = _perturbed(1e-15)
    other["s"]["s11"][1] = 3.0e-17  # rounding-noise entry, relative to its own size it is a different number
    assert_matches_reference(other, REF, "ref", reference_platform=False)


@pytest.mark.parametrize("rel", [1e-9, -1e-9])
def test_tolerance_path_really_fails_on_a_1e_9_perturbation(rel):
    with pytest.raises(AssertionError):
        assert_matches_reference(_perturbed(rel), REF, "ref", reference_platform=False)


def test_reference_platform_path_is_exact():
    assert_matches_reference(copy.deepcopy(REF), REF, "ref", reference_platform=True)
    with pytest.raises(AssertionError):
        assert_matches_reference(_perturbed(1e-15), REF, "ref", reference_platform=True)


def test_non_numeric_content_must_still_match_exactly():
    for mutate in (
        lambda r: r.__setitem__("name", "cel1"),
        lambda r: r.__setitem__("n", 8),
        lambda r: r.__setitem__("ok", False),
        lambda r: r["x"].append(1.0),
        lambda r: r.__setitem__("extra", 1.0),
        lambda r: r.__setitem__("n", 7.0),  # int vs float
    ):
        other = copy.deepcopy(REF)
        mutate(other)
        with pytest.raises(AssertionError):
            assert_matches_reference(other, REF, "ref", reference_platform=False)


def test_text_tolerance_path_parses_numbers_and_still_fails_on_1e_9():
    assert_text_matches_reference(
        TEXT.replace("23.030000000000001", "23.030000000000004"), TEXT, "t", reference_platform=False
    )
    with pytest.raises(AssertionError):
        assert_text_matches_reference(
            TEXT.replace("23.030000000000001", "23.030000023"), TEXT, "t", reference_platform=False
        )


def test_text_tolerance_path_is_exact_for_text_integers_and_printed_digits():
    for bad in (TEXT.replace("FEEDS", "FEEDZ"), TEXT.replace("1430", "1431"), TEXT.replace("64.952", "64.953")):
        with pytest.raises(AssertionError):
            assert_text_matches_reference(bad, TEXT, "t", reference_platform=False)
    assert_text_matches_reference(TEXT.encode(), TEXT.encode(), "t", reference_platform=False)
