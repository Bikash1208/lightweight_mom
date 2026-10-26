"""Phase 20 R5-5: the RIS designer through the app -- design job, named refusals, version/dither options, no heavy compute."""

import time

import pytest

pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from src.app.main import create_app  # noqa: E402
from src.ris import service  # noqa: E402


def _wait(client, job):
    for _ in range(600):
        st = client.get(f"/api/jobs/{job['id']}").json()
        if st["status"] in ("done", "failed", "cancelled"):
            return st
        time.sleep(0.2)
    raise AssertionError("job timed out")


def test_design_job_returns_metrics_labels_and_the_honest_coupling_status():
    client = TestClient(create_app())
    st = _wait(client, client.post("/api/ris/design", json={"theta_s_deg": 30.0, "sensitivity": True}).json())
    assert st["status"] == "done"
    res = st["result"]
    assert res["state"] == "OK" and res["family"] == service.FAMILY and res["validation_label"] == "MODEL_ONLY"
    assert res["metrics"]["gain_at_target_dbi"] == pytest.approx(27.63, abs=0.01)  # R5-3 table, plain 30 deg
    assert "LOCAL_RESPONSE_AT_INCIDENT_ANGLE_ASSUMED" in res["labels"] and "NOT_HFSS_VALIDATED" in res["labels"]
    assert res["coupling_for_this_coding"].startswith("COUPLING_NOT_VERIFIED_FOR_THIS_CODING")
    assert len(res["states"]) == 20 and len(res["states"][0]) == 20 and set(sum(res["states"], [])) <= {1, 2, 3, 4}
    assert set(res["tolerance_sensitivity"]["by_coff_pf"]) == {"0.0225", "0.025", "0.03"}
    assert res["tolerance_sensitivity"]["gain_spread_db"] < 0.3


@pytest.mark.parametrize(
    "body, fragment",
    [
        ({"frequency_ghz": 13.0}, "frequency"),
        ({"theta_s_deg": 55.0}, "0-50"),
        ({"phi_s_deg": 45.0}, "phi_s = 0 and 90"),
        ({"handedness": "lhcp"}, "LHCP"),
        ({"theta_i_deg": 40.0}, "incidence"),
        ({"theta_s_deg": 50.0, "frequency_ghz": 11.9}, "grating"),
    ],
)
def test_requests_outside_the_domain_end_in_the_named_state(body, fragment):
    client = TestClient(create_app())
    st = _wait(client, client.post("/api/ris/design", json=body).json())
    assert st["status"] == "done" and st["result"]["state"] == "NOT_AVAILABLE_FROM_CURRENT_DATA"
    assert any(fragment in why for why in st["result"]["reasons"])


def test_bad_options_are_rejected_and_auto_switches_to_dither_at_35_deg():
    client = TestClient(create_app())
    assert client.post("/api/ris/design", json={"rows": 1}).status_code == 400
    assert client.post("/api/ris/design", json={"variant": "bogus"}).status_code == 422
    res = service.run_design({"theta_s_deg": 40.0, "variant": "auto", "sensitivity": False})
    assert res["variant_used"] == "dither" and any("array-level" in n for n in res["notes"])
    assert service.run_design({"theta_s_deg": 20.0, "variant": "auto", "sensitivity": False})["variant_used"] == "plain"


def test_library_v1_is_selectable_and_agrees_inside_the_hfss_domain():
    a = service.run_design({"theta_s_deg": 30.0, "library": "v1", "sensitivity": False})["metrics"][
        "gain_at_target_dbi"
    ]
    b = service.run_design({"theta_s_deg": 30.0, "library": "v2", "sensitivity": False})["metrics"][
        "gain_at_target_dbi"
    ]
    assert a == pytest.approx(b, abs=1e-6)


def test_defaults_families_and_coupling_endpoints_and_page_served():
    client = TestClient(create_app())
    assert client.get("/api/ris/defaults").json()["frequency_ghz"] == 10.5
    assert service.FAMILY in client.get("/api/ris/families").json()
    assert client.get("/api/ris/coupling").json()["state"].startswith("COUPLING_")
    assert "RIS designer" in client.get("/ris.html").text


def test_library_labels_v2_validated_v1_fails_the_gate_near_phi_90():
    import json
    from pathlib import Path

    from src.hologram import cell_family_registry as registry
    from src.ris.calibrated_library import CalibratedLibrary

    r2 = service.run_design({"theta_s_deg": 30.0, "library": "v2", "sensitivity": False})
    r1 = service.run_design({"theta_s_deg": 30.0, "library": "v1", "sensitivity": False})
    assert r2["library_label"].startswith("R4_GATE_PASS") and r1["library_label"].startswith("FAILS_R4_GATE")
    assert (
        "LIBRARY_V1_FAILS_R4_GATE_NEAR_PHI_90" in r1["cell_response_labels"]
        and "LIBRARY_V1_FAILS_R4_GATE_NEAR_PHI_90" not in r2["cell_response_labels"]
    )
    assert CalibratedLibrary.load(version="v1").response(10.5, 0.0, 0.0, 1).notes[0].startswith("FAILS_R4_GATE")
    assert any("VALIDATED against the R4 gate" in b for b in registry.get(service.FAMILY).known_biases)
    # the recorded server gate results behind the labels
    d = Path(__file__).resolve().parents[1] / "results" / "phase20_r4"
    gate = {
        n: json.loads((d / f"validation_{n}.json").read_text())["gate_theta_le_30"]["result"]
        for n in ("v1_a", "v1_a_phi90", "v2_a", "v2_a_phi90")
    }
    assert gate == {"v1_a": "GATE_PASS", "v1_a_phi90": "GATE_FAIL", "v2_a": "GATE_PASS", "v2_a_phi90": "GATE_PASS"}


def test_page_offers_the_v1_label():
    client = TestClient(create_app())
    assert "fails the R4 gate" in client.get("/ris.html").text
