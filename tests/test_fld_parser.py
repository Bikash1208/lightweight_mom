"""Part 3: the HFSS field-calculator `.fld` parser. Round-trips the user's real example file
(counts, grid, NaN handling, units, ordering), plus synthetic malformed-file and expression tests."""

from pathlib import Path

import numpy as np
import pytest

from src import errors
from src.fields import fld_parser as fp

EXAMPLE_PATH = (
    Path(__file__).resolve().parent.parent / "inputs" / "reference_fields" / "efl_phase_EzZ_example.fld"
)


# ------------------------------------------------------------------------- expression table

@pytest.mark.parametrize("expr, kind, component", [
    ("CmplxPhase(ScalarZ(<Ex,Ey,Ez>))", "phase_rad", "Z"),
    ("CmplxMag(ScalarX(<Ex,Ey,Ez>))", "magnitude", "X"),
    ("Mag(ScalarY(<Ex,Ey,Ez>))", "magnitude", "Y"),
    ("CmplxReal(ScalarZ(<Ex,Ey,Ez>))", "real", "Z"),
    ("CmplxImag(ScalarZ(<Ex,Ey,Ez>))", "imag", "Z"),
])
def test_expression_interpretation_table(expr, kind, component):
    parsed = fp.parse_expression(expr)
    assert parsed.kind == kind
    assert parsed.component == component


def test_unrecognised_expression_raises_typed_error():
    with pytest.raises(errors.LwmomError) as exc_info:
        fp.parse_expression("SomethingElse(ScalarZ(<Ex,Ey,Ez>))")
    assert exc_info.value.code == "FIELD_EXPRESSION_UNSUPPORTED"


# ------------------------------------------------------------------------- real example file

@pytest.mark.skipif(not EXAMPLE_PATH.exists(), reason="example .fld file not present in this checkout")
def test_real_example_file_header():
    grid = fp.read_fld(EXAMPLE_PATH)
    assert grid.header.min_m == pytest.approx((-0.150, -0.150, 1.7e-3))
    assert grid.header.max_m == pytest.approx((0.150, 0.150, 1.7e-3))
    assert grid.header.step_m == pytest.approx((2e-3, 2e-3, 0.0))
    assert grid.header.coordinate_system == "RelativeCS2"
    assert grid.header.unit == "mm"
    assert grid.header.expression.kind == "phase_rad"
    assert grid.header.expression.component == "Z"


@pytest.mark.skipif(not EXAMPLE_PATH.exists(), reason="example .fld file not present in this checkout")
def test_real_example_file_shape_and_ordering():
    grid = fp.read_fld(EXAMPLE_PATH)
    assert grid.value.shape == (151, 151)
    assert grid.x_m.shape == (151, 151)
    assert grid.ordering == "x_major"
    # x varies along axis 1 (columns), y along axis 0 (rows), for a (ny, nx)-shaped array
    assert np.allclose(grid.x_m[0, :], np.linspace(-0.150, 0.150, 151))
    assert np.allclose(grid.y_m[:, 0], np.linspace(-0.150, 0.150, 151))


@pytest.mark.skipif(not EXAMPLE_PATH.exists(), reason="example .fld file not present in this checkout")
def test_real_example_file_nan_at_feed_point():
    grid = fp.read_fld(EXAMPLE_PATH)
    assert grid.nan_mask.sum() == 1
    iy, ix = np.argwhere(grid.nan_mask)[0]
    assert grid.x_m[iy, ix] == pytest.approx(0.0, abs=1e-9)
    assert grid.y_m[iy, ix] == pytest.approx(0.0, abs=1e-9)


@pytest.mark.skipif(not EXAMPLE_PATH.exists(), reason="example .fld file not present in this checkout")
def test_real_example_file_phase_is_on_principal_branch():
    grid = fp.read_fld(EXAMPLE_PATH)
    finite = grid.value[~grid.nan_mask]
    assert np.all(finite >= -np.pi - 1e-6) and np.all(finite <= np.pi + 1e-6)


# ------------------------------------------------------------------------- synthetic / malformed

def _write_fld(tmp_path, header1, header2, rows, name="synthetic.fld"):
    path = tmp_path / name
    with open(path, "w") as fh:
        fh.write(header1 + "\n")
        fh.write(header2 + "\n")
        for row in rows:
            fh.write(" ".join(f"{v:.10e}" for v in row) + "\n")
    return path


def _synthetic_header(nx=3, ny=3, step=2.0):
    lo = -(nx // 2) * step
    hi = (nx // 2) * step
    header1 = f'Grid Output Min: [{lo}mm {lo}mm 1.0mm] Max: [{hi}mm {hi}mm 1.0mm] Grid Size: [{step}mm {step}mm 0mm] CS: "Test" Unit: "mm"'
    header2 = 'X, Y, Z, Scalar data "CmplxPhase(ScalarZ(<Ex,Ey,Ez>))"'
    return header1, header2, lo, hi


def test_synthetic_x_major_roundtrip(tmp_path):
    header1, header2, lo, hi = _synthetic_header(nx=3, ny=3, step=2.0)
    rows = []
    for x in (lo, 0.0, hi):
        for y in (lo, 0.0, hi):
            rows.append((x, y, 1.0, x + y))   # value = x+y (mm units, just a synthetic marker)
    path = _write_fld(tmp_path, header1, header2, rows)
    grid = fp.read_fld(path)
    assert grid.ordering == "x_major"
    assert grid.value.shape == (3, 3)
    # value[iy, ix] should equal x_mm + y_mm at that point
    x_mm, y_mm = grid.x_m * 1e3, grid.y_m * 1e3
    assert np.allclose(grid.value, x_mm + y_mm)


def test_synthetic_y_major_roundtrip(tmp_path):
    header1, header2, lo, hi = _synthetic_header(nx=3, ny=3, step=2.0)
    rows = []
    for y in (lo, 0.0, hi):
        for x in (lo, 0.0, hi):
            rows.append((x, y, 1.0, x + y))
    path = _write_fld(tmp_path, header1, header2, rows)
    grid = fp.read_fld(path)
    assert grid.ordering == "y_major"
    x_mm, y_mm = grid.x_m * 1e3, grid.y_m * 1e3
    assert np.allclose(grid.value, x_mm + y_mm)


def test_nan_value_parses_and_is_masked(tmp_path):
    header1, header2, lo, hi = _synthetic_header(nx=3, ny=3, step=2.0)
    rows = []
    for x in (lo, 0.0, hi):
        for y in (lo, 0.0, hi):
            v = float("nan") if (x, y) == (0.0, 0.0) else x + y
            rows.append((x, y, 1.0, v))
    path = _write_fld(tmp_path, header1, header2, rows)
    grid = fp.read_fld(path)
    assert grid.nan_mask.sum() == 1


def test_malformed_row_count_raises(tmp_path):
    header1, header2, lo, hi = _synthetic_header(nx=3, ny=3, step=2.0)
    rows = [(lo, lo, 1.0, 0.0), (lo, 0.0, 1.0, 0.0)]   # only 2 of the declared 9 rows
    path = _write_fld(tmp_path, header1, header2, rows)
    with pytest.raises(ValueError, match="row count"):
        fp.read_fld(path)


def test_malformed_header_raises(tmp_path):
    path = tmp_path / "bad.fld"
    path.write_text("not a valid header line at all\nX, Y, Z, Scalar data \"CmplxPhase(ScalarZ(<Ex,Ey,Ez>))\"\n")
    with pytest.raises(ValueError, match="header"):
        fp.read_fld(path)


def test_malformed_data_row_raises(tmp_path):
    header1, header2, lo, hi = _synthetic_header(nx=1, ny=1, step=2.0)
    path = tmp_path / "bad2.fld"
    with open(path, "w") as fh:
        fh.write(header1 + "\n" + header2 + "\n")
        fh.write("1.0 2.0 3.0\n")   # only 3 columns, not 4
    with pytest.raises(ValueError, match="malformed data row"):
        fp.read_fld(path)
