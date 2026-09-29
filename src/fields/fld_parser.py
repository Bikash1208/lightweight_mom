"""Parses the HFSS field-calculator ``.fld`` export format (PHASE13_WEBAPP_PROMPT.md "Reference-
field import"), e.g. ``inputs/reference_fields/efl_phase_EzZ_example.fld`` (151 x 151 points,
2 mm pitch, +-150 mm, z = 1.7 mm in "RelativeCS2", 22803 lines):

```
Grid Output Min: [-150mm -150mm 1.7mm] Max: [150mm 150mm 1.7mm] Grid Size: [2mm 2mm 0mm] CS: "RelativeCS2" Unit: "mm"
X, Y, Z, Scalar data "CmplxPhase(ScalarZ(<Ex,Ey,Ez>))"
-1.5000000000000000e+02 -1.5000000000000000e+02 1.7000000000000000e+00  -3.1073796258277993e+00
...
```

Line 1: grid extent/step (with per-axis units, e.g. ``-150mm``), the named coordinate system, and
the display unit. Line 2: the field-calculator expression that produced the scalar column. Line
3+: whitespace-separated ``X Y Z value`` rows. Ordering (X- or Y-major) is DETECTED from the data,
never assumed (PHASE13_WEBAPP_PROMPT.md's explicit instruction).

Recognised expressions (``FIELD_EXPRESSION_UNSUPPORTED`` -- never a guess -- for anything else):
``CmplxPhase(...)`` -> radians on (-pi, pi]; ``CmplxMag(...)``/``Mag(...)`` -> magnitude;
``CmplxReal(...)``/``CmplxImag(...)`` -> real/imaginary part; the inner ``ScalarX/Y/Z(<Ex,Ey,Ez>)``
names which field component the scalar column is (recorded, not further interpreted here).
"""

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..errors import field_expression_unsupported

_LENGTH_TOKEN_RE = re.compile(r"([-+0-9.eE]+)\s*([a-zA-Z]+)")
_TRIPLE_RE = re.compile(r"\[\s*([^\]]+?)\s*\]")
_HEADER1_RE = re.compile(
    r"Grid Output Min:\s*\[(?P<min>[^\]]+)\]\s*Max:\s*\[(?P<max>[^\]]+)\]\s*"
    r"Grid Size:\s*\[(?P<step>[^\]]+)\]\s*CS:\s*\"(?P<cs>[^\"]*)\"\s*Unit:\s*\"(?P<unit>[^\"]*)\""
)
_EXPRESSION_RE = re.compile(r'"(?P<expr>[^"]+)"\s*$')
_OUTER_FN_RE = re.compile(r"^(?P<fn>CmplxPhase|CmplxMag|CmplxReal|CmplxImag|Mag)\((?P<inner>.*)\)$")
_SCALAR_RE = re.compile(r"^Scalar(?P<comp>[XYZ])\(")

_UNIT_TO_M = {"mm": 1e-3, "m": 1.0, "cm": 1e-2, "um": 1e-6}

_OUTER_KIND = {
    "CmplxPhase": "phase_rad", "CmplxMag": "magnitude", "Mag": "magnitude",
    "CmplxReal": "real", "CmplxImag": "imag",
}


def _length_to_m(token):
    m = _LENGTH_TOKEN_RE.match(token.strip())
    if not m:
        raise ValueError(f".fld header: unparseable length token {token!r}")
    value, unit = float(m.group(1)), m.group(2)
    if unit not in _UNIT_TO_M:
        raise ValueError(f".fld header: unknown length unit {unit!r} in {token!r}")
    return value * _UNIT_TO_M[unit]


def _parse_triple_m(text):
    parts = text.split()
    if len(parts) != 3:
        raise ValueError(f".fld header: expected 3 values, got {text!r}")
    return tuple(_length_to_m(p) for p in parts)


@dataclass(frozen=True)
class FldExpression:
    raw: str
    kind: str          # "phase_rad" | "magnitude" | "real" | "imag"
    component: "str | None"   # "X" | "Y" | "Z" | None


def parse_expression(expr):
    """Interprets a field-calculator expression string; raises `errors.LwmomError`
    (`FIELD_EXPRESSION_UNSUPPORTED`) for anything not in the recognised table -- never a guess."""
    m = _OUTER_FN_RE.match(expr.strip())
    if not m:
        raise field_expression_unsupported(expr)
    fn, inner = m.group("fn"), m.group("inner")
    kind = _OUTER_KIND[fn]
    comp_match = _SCALAR_RE.match(inner.strip())
    component = comp_match.group("comp") if comp_match else None
    return FldExpression(raw=expr, kind=kind, component=component)


@dataclass(frozen=True)
class FldHeader:
    min_m: tuple
    max_m: tuple
    step_m: tuple
    coordinate_system: str
    unit: str
    expression: FldExpression


def _parse_header(line1, line2):
    m1 = _HEADER1_RE.search(line1)
    if not m1:
        raise ValueError(f".fld file: unparseable header line 1: {line1!r}")
    m2 = _EXPRESSION_RE.search(line2)
    if not m2:
        raise ValueError(f".fld file: unparseable header line 2 (no quoted expression): {line2!r}")
    return FldHeader(
        min_m=_parse_triple_m(m1.group("min")), max_m=_parse_triple_m(m1.group("max")),
        step_m=_parse_triple_m(m1.group("step")), coordinate_system=m1.group("cs"), unit=m1.group("unit"),
        expression=parse_expression(m2.group("expr")),
    )


@dataclass(frozen=True)
class FldGrid:
    header: FldHeader
    x_m: np.ndarray    # 2-D, shape (ny, nx)
    y_m: np.ndarray    # 2-D, shape (ny, nx)
    value: np.ndarray  # 2-D, shape (ny, nx); NaN where the source file had NaN (e.g. the feed point)
    ordering: str       # "x_major" | "y_major" -- detected, not assumed
    path: str

    @property
    def nan_mask(self):
        return np.isnan(self.value)


def read_fld(path):
    """Reads one `.fld` file end to end. Raises `ValueError` for a malformed file (wrong row
    count, unparseable header, ...) and `errors.LwmomError` (`FIELD_EXPRESSION_UNSUPPORTED`) for
    an unrecognised field-calculator expression."""
    path = Path(path)
    with open(path, encoding="utf-8") as fh:
        line1 = fh.readline()
        line2 = fh.readline()
        header = _parse_header(line1, line2)
        rows = fh.readlines()

    xs, ys, vals = [], [], []
    for line in rows:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != 4:
            raise ValueError(f"{path}: malformed data row (expected 4 columns, got {len(parts)}): {line!r}")
        xs.append(float(parts[0]))
        ys.append(float(parts[1]))
        vals.append(float(parts[3]))   # parts[2] is Z (the fixed sampling-plane height)
    xs, ys, vals = np.array(xs), np.array(ys), np.array(vals)

    dx_m, dy_m, _ = header.step_m
    nx = round((header.max_m[0] - header.min_m[0]) / dx_m) + 1 if dx_m else 1
    ny = round((header.max_m[1] - header.min_m[1]) / dy_m) + 1 if dy_m else 1
    if nx * ny != xs.size:
        raise ValueError(
            f"{path}: row count {xs.size} does not match the header's declared grid "
            f"({nx} x {ny} = {nx * ny}); file may be truncated or the header malformed"
        )

    # Detect ordering from the data itself (brief: "detect ordering, don't assume it"): if the
    # first `ny` rows share the same X value, X is the outer (slower-varying) loop.
    ordering = "x_major" if ny > 1 and np.isclose(xs[0], xs[min(ny - 1, xs.size - 1)]) else "y_major"
    x_mm_unit = 1e-3 if header.unit == "mm" else _UNIT_TO_M.get(header.unit, 1.0)
    xs_m, ys_m = xs * x_mm_unit, ys * x_mm_unit

    if ordering == "x_major":
        x2d, y2d, v2d = (a.reshape(nx, ny).T for a in (xs_m, ys_m, vals))
    else:
        x2d, y2d, v2d = (a.reshape(ny, nx) for a in (xs_m, ys_m, vals))

    return FldGrid(header=header, x_m=x2d, y_m=y2d, value=v2d, ordering=ordering, path=str(path))
