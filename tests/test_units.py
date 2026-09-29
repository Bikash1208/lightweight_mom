import math

import pytest

from src import units


def test_mm_and_ghz_round_trip():
    assert units.from_mm(1.575) == pytest.approx(1.575e-3, rel=1e-15)
    assert units.to_mm(1.575e-3) == pytest.approx(1.575, rel=1e-15)
    assert units.to_mm(units.from_mm(7.123)) == pytest.approx(7.123, rel=1e-14)
    assert units.from_ghz(10.0) == 1e10
    assert units.to_ghz(units.from_ghz(2.45)) == pytest.approx(2.45, rel=1e-14)


def test_degree_radian():
    assert units.from_deg(180.0) == pytest.approx(math.pi, rel=1e-15)
    assert units.from_deg(90.0) == pytest.approx(math.pi / 2, rel=1e-15)
    assert units.to_deg(units.from_deg(37.5)) == pytest.approx(37.5, rel=1e-14)


@pytest.mark.parametrize(
    "text, metres",
    [
        ("8 mm", 8e-3),
        ("8mm", 8e-3),
        ("1.5 cm", 1.5e-2),
        ("25 um", 25e-6),
        ("25 µm", 25e-6),
        ("0.5 m", 0.5),
        ("8e-3", 8e-3),        # what PyYAML gives for `8e-3`: bare number string => SI
        ("  8.0E-3  ", 8e-3),
        ("-1.5 mm", -1.5e-3),
        (8e-3, 8e-3),
        (2, 2.0),
    ],
)
def test_parse_length(text, metres):
    assert units.parse_length(text) == pytest.approx(metres, rel=1e-15)


@pytest.mark.parametrize(
    "text, hertz",
    [
        ("10 GHz", 10e9),
        ("2.4GHz", 2.4e9),
        ("500 MHz", 500e6),
        ("3 kHz", 3e3),
        ("1 THz", 1e12),
        ("60 Hz", 60.0),
        ("10e9", 10e9),
        (10e9, 10e9),
    ],
)
def test_parse_frequency(text, hertz):
    assert units.parse_frequency(text) == pytest.approx(hertz, rel=1e-15)


@pytest.mark.parametrize(
    "bad",
    ["8 Mm", "8 MM", "8 furlong", "abc", "", "8 mm mm", True, None, float("nan"),
     float("inf"), "nan", "inf", 1 + 2j, [8e-3]],
)
def test_parse_length_rejects(bad):
    with pytest.raises(ValueError):
        units.parse_length(bad)


@pytest.mark.parametrize("bad", ["5 mHz", "10 ghz", "10 GHZ", "10 m", "GHz", -float("inf")])
def test_parse_frequency_rejects_case_and_unknown_units(bad):
    # 'mHz' (milli) vs 'MHz' (mega) and 'mm' vs 'Mm' must never be confused.
    with pytest.raises(ValueError):
        units.parse_frequency(bad)
