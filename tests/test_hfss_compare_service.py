"""Part 3: the HFSS-compare service (upload validation, comparison, ledger)."""

import pytest

from src.hfss_compare import service as svc
from src.phase9.hfss_reference import TouchstoneError
from src.phase12d.hfss_ingest import TOUCHSTONE_PATH

pytestmark = pytest.mark.skipif(not TOUCHSTONE_PATH.exists(), reason="real HFSS export not present in this checkout")


def test_validate_touchstone_upload_accepts_real_file():
    net = svc.validate_touchstone_upload(TOUCHSTONE_PATH)
    assert net.n_ports == 36


def test_validate_touchstone_upload_rejects_malformed_file(tmp_path):
    bad = tmp_path / "bad.s2p"
    bad.write_text("not a touchstone file\n")
    with pytest.raises(TouchstoneError):
        svc.validate_touchstone_upload(bad)


def test_compare_uploaded_file_reproduces_g1():
    summary, comparison = svc.compare_uploaded_file("G1", TOUCHSTONE_PATH)
    assert summary.meets_tolerance is True
    assert summary.mean_dphase_deg == pytest.approx(-0.2935, abs=1e-3)
    assert "HFSS_VALIDATED_OFF_RESONANCE" in summary.labels


def test_ledger_roundtrip(tmp_path):
    summary, _ = svc.compare_uploaded_file("G1", TOUCHSTONE_PATH)
    ledger_path = tmp_path / "ledger.json"
    svc.record_comparison(summary, ledger_path)
    reloaded = svc.load_ledger(ledger_path)
    assert len(reloaded) == 1
    assert reloaded[0]["case_id"] == "G1"

    svc.record_comparison(summary, ledger_path)
    assert len(svc.load_ledger(ledger_path)) == 2
