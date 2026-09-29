"""Part 3: HFSS eigenmode dispersion table import."""

import pytest

from src import errors
from src.dispersion import hfss_table as ht


def _write_csv(tmp_path, rows, columns=("geometry_param", "beta_rad_per_m", "x_s_ohm")):
    path = tmp_path / "dispersion.csv"
    with open(path, "w") as fh:
        fh.write(",".join(columns) + "\n")
        for row in rows:
            fh.write(",".join(str(v) for v in row) + "\n")
    return path


def test_reads_a_well_formed_table(tmp_path):
    path = _write_csv(tmp_path, [(1.0, 100.0, 50.0), (2.0, 150.0, 60.0), (3.0, 220.0, 75.0)])
    table = ht.read_dispersion_csv(path, frequency_hz=10e9)
    assert table.label == "DISPERSION_FROM_HFSS"
    assert table.geometry_range == (1.0, 3.0)
    assert len(table.samples) == 3


def test_beta_interpolation_is_linear(tmp_path):
    path = _write_csv(tmp_path, [(1.0, 100.0, 50.0), (3.0, 300.0, 50.0)])
    table = ht.read_dispersion_csv(path, frequency_hz=10e9)
    assert table.beta_at(2.0) == pytest.approx(200.0)
    assert table.beta_at(1.0) == pytest.approx(100.0)
    assert table.beta_at(3.0) == pytest.approx(300.0)


def test_beta_out_of_range_raises_impedance_unreachable(tmp_path):
    path = _write_csv(tmp_path, [(1.0, 100.0, 50.0), (3.0, 300.0, 50.0)])
    table = ht.read_dispersion_csv(path, frequency_hz=10e9)
    with pytest.raises(errors.LwmomError) as exc_info:
        table.beta_at(5.0)
    assert exc_info.value.code == "IMPEDANCE_UNREACHABLE"


def test_inverse_lookup_finds_nearest_achievable(tmp_path):
    path = _write_csv(tmp_path, [(1.0, 100.0, 50.0), (2.0, 150.0, 60.0), (3.0, 220.0, 75.0)])
    table = ht.read_dispersion_csv(path, frequency_hz=10e9)
    gp, achieved, within = table.geometry_for_beta(148.0)
    assert gp == 2.0 and within is True
    gp2, achieved2, within2 = table.geometry_for_beta(1000.0)
    assert within2 is False   # out of range, but still returns the nearest achievable


def test_missing_required_column_raises_config_invalid(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("foo,bar\n1,2\n")
    with pytest.raises(errors.LwmomError) as exc_info:
        ht.read_dispersion_csv(path, frequency_hz=10e9)
    assert exc_info.value.code == "CONFIG_INVALID"


def test_neither_beta_nor_xs_column_raises(tmp_path):
    path = tmp_path / "bad2.csv"
    path.write_text("geometry_param\n1\n2\n")
    with pytest.raises(errors.LwmomError):
        ht.read_dispersion_csv(path, frequency_hz=10e9)


def test_non_numeric_value_raises(tmp_path):
    path = _write_csv(tmp_path, [(1.0, "oops", 50.0)])
    with pytest.raises(errors.LwmomError) as exc_info:
        ht.read_dispersion_csv(path, frequency_hz=10e9)
    assert exc_info.value.code == "CONFIG_INVALID"


def test_duplicate_geometry_param_raises(tmp_path):
    path = _write_csv(tmp_path, [(1.0, 100.0, 50.0), (1.0, 999.0, 50.0)])
    with pytest.raises(errors.LwmomError) as exc_info:
        ht.read_dispersion_csv(path, frequency_hz=10e9)
    assert exc_info.value.code == "CONFIG_INVALID"


def test_missing_file_raises_config_invalid(tmp_path):
    with pytest.raises(errors.LwmomError) as exc_info:
        ht.read_dispersion_csv(tmp_path / "nope.csv", frequency_hz=10e9)
    assert exc_info.value.code == "CONFIG_INVALID"


def test_x_s_only_table_works_without_beta(tmp_path):
    path = _write_csv(tmp_path, [(1.0, 50.0), (2.0, 60.0)], columns=("geometry_param", "x_s_ohm"))
    table = ht.read_dispersion_csv(path, frequency_hz=10e9)
    assert table.x_s_at(1.5) == pytest.approx(55.0)
    with pytest.raises(ValueError, match="beta_rad_per_m"):
        table.beta_at(1.5)
