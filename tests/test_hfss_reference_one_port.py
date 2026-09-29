"""Phase 12C: one-port (grounded) HFSS CSV ingestion, additive to Phase 9's own two-port schema."""

import textwrap

import pytest

from src.phase9.hfss_reference import HfssReferenceError, load_hfss_reference


def _write(tmp_path, name, content):
    path = tmp_path / name
    path.write_text(textwrap.dedent(content))
    return path


def test_two_port_csv_is_unchanged(tmp_path):
    path = _write(tmp_path, "two_port.csv", """\
        cell_id,L_mm,W_mm,frequency_GHz,S11_real,S11_imag,S21_real,S21_imag
        CELL-01,1.0,1.0,10.0,0.1,0.2,0.3,0.4
        """)
    records = load_hfss_reference(path)
    assert len(records) == 1
    assert records[0].s21 == complex(0.3, 0.4)


def test_one_port_csv_gives_none_s21(tmp_path):
    path = _write(tmp_path, "one_port.csv", """\
        cell_id,L_mm,W_mm,frequency_GHz,S11_real,S11_imag
        G-01,6.0,6.0,12.5,-0.95,0.1
        """)
    records = load_hfss_reference(path)
    assert len(records) == 1
    assert records[0].s11 == complex(-0.95, 0.1)
    assert records[0].s21 is None


def test_partial_s21_columns_are_rejected(tmp_path):
    path = _write(tmp_path, "partial.csv", """\
        cell_id,L_mm,W_mm,frequency_GHz,S11_real,S11_imag,S21_real
        BAD-01,1.0,1.0,10.0,0.1,0.2,0.3
        """)
    with pytest.raises(HfssReferenceError, match="partial S21"):
        load_hfss_reference(path)


def test_missing_file_returns_none(tmp_path):
    assert load_hfss_reference(tmp_path / "does_not_exist.csv") is None


def test_missing_common_column_still_raises(tmp_path):
    path = _write(tmp_path, "bad.csv", """\
        cell_id,L_mm,frequency_GHz,S11_real,S11_imag
        BAD-01,1.0,10.0,0.1,0.2
        """)
    with pytest.raises(HfssReferenceError, match="missing required"):
        load_hfss_reference(path)
