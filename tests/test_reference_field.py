"""Phase 8A Sec. 7: generic reference-field and beam-definition input validation."""

import pytest

from src.reference_field import validate_beam_definition, validate_reference_field

_HEADER = "x_mm,y_mm,frequency_GHz,Ex_real,Ex_imag,Ey_real,Ey_imag\n"


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_valid_reference_field_parses_cleanly(tmp_path):
    text = _HEADER + "0,0,12.5,1,0,0,0\n" + "2,0,12.5,0.9,0.1,0,0\n" + "0,2,12.5,0.8,0.2,0,0\n"
    path = _write(tmp_path, "ref.csv", text)
    result = validate_reference_field(path, expected_frequency_ghz=12.5)
    assert result.ok, result.errors
    assert result.n_samples == 3
    assert result.dx_mm == pytest.approx(2.0)
    assert result.dy_mm == pytest.approx(2.0)


def test_missing_file_is_reported_not_raised(tmp_path):
    result = validate_reference_field(tmp_path / "missing.csv")
    assert not result.ok
    assert "not found" in result.errors[0]


def test_missing_required_column_is_detected(tmp_path):
    path = _write(tmp_path, "ref.csv", "x_mm,y_mm,frequency_GHz\n0,0,12.5\n")
    result = validate_reference_field(path)
    assert not result.ok
    assert "missing required column" in result.errors[0]


def test_nan_value_is_detected(tmp_path):
    text = _HEADER + "0,0,12.5,nan,0,0,0\n"
    path = _write(tmp_path, "ref.csv", text)
    result = validate_reference_field(path)
    assert not result.ok
    assert any("NaN/Inf" in e for e in result.errors)


def test_duplicate_sample_is_detected(tmp_path):
    text = _HEADER + "0,0,12.5,1,0,0,0\n" + "0,0,12.5,0.5,0,0,0\n"
    path = _write(tmp_path, "ref.csv", text)
    result = validate_reference_field(path)
    assert not result.ok
    assert any("duplicate sample" in e for e in result.errors)


def test_frequency_mismatch_is_detected(tmp_path):
    text = _HEADER + "0,0,10.0,1,0,0,0\n"
    path = _write(tmp_path, "ref.csv", text)
    result = validate_reference_field(path, expected_frequency_ghz=12.5)
    assert not result.ok
    assert any("no samples at the requested frequency" in e for e in result.errors)


def test_empty_file_is_detected(tmp_path):
    path = _write(tmp_path, "ref.csv", _HEADER)
    result = validate_reference_field(path)
    assert not result.ok


# ---------------------------------------------------------------------------- beam definition


def _beam(**overrides):
    beam = {"theta_deg": 20.0, "phi_deg": 30.0, "weight": 1.0, "initial_phase_deg": 0.0, "polarization": "x"}
    beam.update(overrides)
    return beam


def test_valid_beam_definition_passes():
    result = validate_beam_definition(_beam())
    assert result.ok, result.errors


def test_beam_missing_field_is_detected():
    beam = _beam()
    del beam["theta_deg"]
    result = validate_beam_definition(beam)
    assert not result.ok
    assert "missing beam field" in result.errors[0]


def test_beam_theta_out_of_range_is_detected():
    result = validate_beam_definition(_beam(theta_deg=120.0))
    assert not result.ok
    assert any("theta_deg" in e for e in result.errors)


def test_beam_negative_weight_is_detected():
    result = validate_beam_definition(_beam(weight=-1.0))
    assert not result.ok
    assert any("weight" in e for e in result.errors)


def test_beam_unsupported_polarization_is_detected():
    result = validate_beam_definition(_beam(polarization="circular"))
    assert not result.ok
    assert any("polarization" in e for e in result.errors)
