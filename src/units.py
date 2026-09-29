"""Unit handling at the input/output boundary.

Internally every module uses SI: metres, hertz, seconds, radians. This module
is the *only* place where mm / GHz / degrees are converted.

Input rules for config values:

* A plain number (or a bare-number string such as "8e-3") is already SI.
* A string with a unit suffix ("8 mm", "10GHz") is converted to SI.
* Unit suffixes are case-sensitive on purpose: "mm" is millimetre but "Mm" is
  megametre, and "mHz" is millihertz but "MHz" is megahertz.
"""

import math
import numbers
import re

_LENGTH_UNITS = {"m": 1.0, "cm": 1e-2, "mm": 1e-3, "um": 1e-6, "µm": 1e-6}
_FREQUENCY_UNITS = {"Hz": 1.0, "kHz": 1e3, "MHz": 1e6, "GHz": 1e9, "THz": 1e12}

_QUANTITY_RE = re.compile(
    r"^\s*(?P<num>[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*(?P<unit>\S*)\s*$"
)


def parse_quantity(value, unit_table, what):
    """Convert a config value to SI using ``unit_table`` (unit -> SI factor).

    ``what`` names the quantity in error messages. Returns a finite float.
    """
    if isinstance(value, bool):
        raise ValueError(f"{what}: expected a number or string, got bool")
    if isinstance(value, numbers.Real):
        result = float(value)
    elif isinstance(value, str):
        match = _QUANTITY_RE.match(value)
        if match is None:
            raise ValueError(f"{what}: cannot parse {value!r}")
        unit = match.group("unit")
        if unit == "":
            factor = 1.0  # bare number => SI
        elif unit in unit_table:
            factor = unit_table[unit]
        else:
            raise ValueError(
                f"{what}: unknown unit {unit!r} in {value!r}; "
                f"accepted: {sorted(unit_table)} (case-sensitive)"
            )
        result = float(match.group("num")) * factor
    else:
        raise ValueError(f"{what}: expected a number or string, got {type(value).__name__}")
    if not math.isfinite(result):
        raise ValueError(f"{what}: value must be finite, got {result!r}")
    return result


def parse_length(value, what="length"):
    """Config length -> metres."""
    return parse_quantity(value, _LENGTH_UNITS, what)


def parse_frequency(value, what="frequency"):
    """Config frequency -> hertz."""
    return parse_quantity(value, _FREQUENCY_UNITS, what)


# ---- explicit output-boundary conversions (SI -> display units and back) ----

def from_mm(x_mm):
    return x_mm * 1e-3


def to_mm(x_m):
    return x_m * 1e3


def from_ghz(f_ghz):
    return f_ghz * 1e9


def to_ghz(f_hz):
    return f_hz * 1e-9


def from_deg(angle_deg):
    return angle_deg * (math.pi / 180.0)


def to_deg(angle_rad):
    return angle_rad * (180.0 / math.pi)
