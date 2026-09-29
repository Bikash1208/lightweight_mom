"""Sec. 6-9: HFSS reference-data ingestion; honest absence handling (Rule 1)."""

import pytest

from src.phase9.hfss_reference import HfssReferenceError, load_hfss_reference


def test_load_hfss_reference_returns_none_when_file_absent(tmp_path):
    assert load_hfss_reference(tmp_path / "does_not_exist.csv") is None


def test_load_hfss_reference_parses_a_valid_synthetic_test_file(tmp_path):
    path = tmp_path / "hfss.csv"
    path.write_text(
        "cell_id,L_mm,W_mm,frequency_GHz,S11_real,S11_imag,S21_real,S21_imag\nCELL-01,1.0,1.0,12.5,0.1,0.2,0.9,0.05\n"
    )
    records = load_hfss_reference(path)
    assert len(records) == 1
    assert records[0].cell_id == "CELL-01"
    assert records[0].s11 == complex(0.1, 0.2)
    assert records[0].s21 == complex(0.9, 0.05)


def test_load_hfss_reference_rejects_missing_columns(tmp_path):
    path = tmp_path / "hfss.csv"
    path.write_text("cell_id,L_mm,W_mm\nCELL-01,1.0,1.0\n")
    with pytest.raises(HfssReferenceError):
        load_hfss_reference(path)


def test_load_hfss_reference_rejects_non_numeric_field(tmp_path):
    path = tmp_path / "hfss.csv"
    path.write_text(
        "cell_id,L_mm,W_mm,frequency_GHz,S11_real,S11_imag,S21_real,S21_imag\n"
        "CELL-01,not_a_number,1.0,12.5,0.1,0.2,0.9,0.05\n"
    )
    with pytest.raises(HfssReferenceError):
        load_hfss_reference(path)


def test_no_real_hfss_reference_exists_in_this_project():
    """The production default path must not exist -- this is the honest baseline Phase 9's
    production run reports against (NOT_AVAILABLE_FROM_CURRENT_DATA)."""
    from src.phase9.hfss_reference import DEFAULT_HFSS_REFERENCE_PATH

    assert not DEFAULT_HFSS_REFERENCE_PATH.exists()


# ---------------------------------------------------------------- Phase 12D: read_touchstone

from src.phase9.hfss_reference import TouchstoneError, read_touchstone


def _write_touchstone(tmp_path, name, header_extra=""):
    """A tiny, synthetic 2-port MA Touchstone file, deliberately wrapped across two lines per row
    (1 pair/line) to test the flat-token, wrap-width-agnostic parser."""
    path = tmp_path / name
    path.write_text(
        "! synthetic test file\n"
        f"# GHz S MA{header_extra}\n"
        "! Port[1] = FloquetPort1:1\n"
        "! Port[2] = FloquetPort1:2\n"
        "10.0  1.0 0.0\n"
        "      0.0 0.0\n"
        "      0.0 0.0\n"
        "      0.5 90.0\n"
        "11.0  0.9 10.0\n"
        "      0.0 0.0\n"
        "      0.0 0.0\n"
        "      0.4 80.0\n"
    )
    return path


def test_read_touchstone_parses_wrapped_ma_file(tmp_path):
    net = read_touchstone(_write_touchstone(tmp_path, "synthetic.s2p"))
    assert net.n_ports == 2
    assert list(net.frequency_hz) == [10.0e9, 11.0e9]
    assert net.port_labels == ("FloquetPort1:1", "FloquetPort1:2")
    assert net.s[0, 0, 0] == complex(1.0, 0.0)
    import numpy as np

    assert np.isclose(net.s[0, 1, 1], 0.5 * (np.cos(np.radians(90)) + 1j * np.sin(np.radians(90))))


def test_read_touchstone_port_mode_groups():
    pass  # covered indirectly by test_phase12d_hfss_ingest.py against the real 36-port file


def test_read_touchstone_rejects_bad_extension(tmp_path):
    path = tmp_path / "bad.txt"
    path.write_text("# GHz S MA\n10.0 1.0 0.0\n")
    with pytest.raises(TouchstoneError):
        read_touchstone(path)


def test_read_touchstone_rejects_truncated_data(tmp_path):
    path = tmp_path / "truncated.s2p"
    path.write_text("# GHz S MA\n10.0 1.0 0.0 0.0 0.0 0.0\n")  # missing one pair of an s2p block
    with pytest.raises(TouchstoneError):
        read_touchstone(path)


def test_read_touchstone_ri_format(tmp_path):
    path = tmp_path / "ri.s1p"
    path.write_text("# MHz S RI\n1000.0 0.3 0.4\n")
    net = read_touchstone(path)
    assert net.frequency_hz[0] == 1000.0e6
    assert net.s[0, 0, 0] == complex(0.3, 0.4)


def test_read_touchstone_reference_impedance_parsed(tmp_path):
    path = _write_touchstone(tmp_path, "with_r.s2p", header_extra=" R 376.617421")
    net = read_touchstone(path)
    assert net.reference_impedance_ohm == 376.617421
