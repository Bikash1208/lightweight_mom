"""Phase 13 Part 3: the HFSS compare page's API -- uploads the real G1 Touchstone export and
compares it against the real MoM reference. Needs the optional 'app' extra AND the real HFSS
export (excluded from git, `hfss_data`-marked -- skipped entirely when absent, same convention
tests/test_phase12d_hfss_ingest.py uses)."""

import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from src.app.main import create_app  # noqa: E402
from src.phase12d import hfss_ingest as hi  # noqa: E402

pytestmark = [
    pytest.mark.hfss_data,
    pytest.mark.skipif(not hi.TOUCHSTONE_PATH.exists(), reason="real HFSS export not present in this checkout"),
]


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("LWMOM_APP_DATA_DIR", str(tmp_path / "lwmom_data"))
    monkeypatch.setattr("src.hfss_compare.service.DEFAULT_LEDGER_PATH", tmp_path / "ledger.json")
    return TestClient(create_app())


def test_cases_lists_g1_as_the_one_with_real_data(client):
    resp = client.get("/api/hfss_compare/cases")
    assert resp.status_code == 200
    body = resp.json()
    assert "G1" in body["case_ids"]
    assert body["case_with_real_hfss_data_on_disk"] == "G1"


def test_upload_g1_touchstone_reproduces_the_measured_agreement(client):
    with open(hi.TOUCHSTONE_PATH, "rb") as fh:
        resp = client.post(
            "/api/hfss_compare/upload",
            files={"touchstone_file": ("unit_HFSSDesign1.s36p", fh, "application/octet-stream")},
            data={"case_id": "G1"},
        )
    assert resp.status_code == 200, resp.text
    body = resp.json()
    summary = body["summary"]
    # The brief's own already-measured G1 result: mean phase diff ~ -0.29 deg, within tolerance.
    assert abs(summary["mean_dphase_deg"]) < 1.0
    assert summary["meets_tolerance"] is True
    assert "HFSS_VALIDATED_OFF_RESONANCE" in summary["labels"]
    assert len(body["overlay"]) > 0

    ledger_resp = client.get("/api/hfss_compare/ledger")
    assert len(ledger_resp.json()) == 1


def test_upload_rejects_a_malformed_touchstone(client, tmp_path):
    bad = tmp_path / "bad.s36p"
    bad.write_text("not a touchstone file\n")
    with open(bad, "rb") as fh:
        resp = client.post(
            "/api/hfss_compare/upload",
            files={"touchstone_file": ("bad.s36p", fh, "application/octet-stream")},
            data={"case_id": "G1"},
        )
    assert resp.status_code == 400
