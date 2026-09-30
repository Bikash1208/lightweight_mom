"""Phase 15 Part D: the app's API surface for per-beam polarization/OAM, the vector pattern model,
multi-feed synthesis, and feed markers in Export -- end to end through the real HTTP app, matching
`test_app_hmsa.py`'s own fixtures/conventions (duplicated here per this project's own pattern of
per-file `client`/`grounded_library` fixtures, not a shared conftest)."""

import time
import uuid
from pathlib import Path

import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from src.app.main import create_app  # noqa: E402

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3}
_SMALL_MESH = {"cells_x": 8, "cells_y": 8, "floquet_order": 8}
_EXAMPLE_FLD = Path(__file__).resolve().parent.parent / "inputs" / "reference_fields" / "efl_phase_EzZ_example.fld"


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("LWMOM_APP_DATA_DIR", str(tmp_path / "lwmom_data"))
    return TestClient(create_app())


def _wait_for_job(client, job_id, timeout=60.0):
    deadline = time.time() + timeout
    while True:
        resp = client.get(f"/api/jobs/{job_id}")
        job = resp.json()
        if job["status"] in ("done", "failed", "cancelled"):
            return job
        if time.time() > deadline:
            raise TimeoutError(f"job {job_id} did not finish in {timeout}s: {job}")
        time.sleep(0.05)


@pytest.fixture
def grounded_library(client):
    library_id = f"test-phase15-{uuid.uuid4().hex[:8]}"
    body = {
        "library_id": library_id,
        "base_config": {
            "period_x": 8.0e-3,
            "period_y": 8.0e-3,
            "substrate": {**_SUBSTRATE, "backing": "pec"},
            "geometry": {"type": "rectangle", "length": 2.0e-3, "width": 2.0e-3},
            "excitation": {"frequency": 10.54e9},
            "mesh": _SMALL_MESH,
            "reference_planes": {"top_air_thickness": 5.0e-3},
        },
        "l_values_mm": [2.0, 3.0, 4.0],
        "w_values_mm": [2.0, 3.0, 4.0],
        "f_min_hz": 10.54e9,
        "f_max_hz": 10.54e9,
        "f_step_hz": 1.0e9,
        "f_design_hz": 10.54e9,
    }
    resp = client.post("/api/library/sweep", json=body)
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "done", job
    return library_id


def _upload_example_field(client):
    with open(_EXAMPLE_FLD, "rb") as fh:
        resp = client.post("/api/hmsa/upload_field", files={"phase_file": ("efl_phase_EzZ_example.fld", fh, "text/plain")})
    return resp.json()["upload_id"]


# ============================================================== Part A: polarization/OAM via API


def test_synthesize_with_desired_polarization_and_oam_charge(client, grounded_library):
    upload_id = _upload_example_field(client)
    body = {
        "upload_id": upload_id,
        "formulation": "B",
        "frequency_hz": 10.54e9,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0, "desired_polarization": "rhcp", "oam_charge": 0}],
        "period_mm": 8.0,
        "library_id": grounded_library,
    }
    resp = client.post("/api/hmsa/synthesize", json=body)
    job = _wait_for_job(client, resp.json()["id"], timeout=120.0)
    assert job["status"] == "done", job
    # feed position auto-detected from the example file's own NaN, labelled
    assert any(label.startswith("FEED_POSITION_AUTO_DETECTED(") for label in job["result"]["labels"])
    assert job["result"]["feeds_mm"] == [[0.0, 0.0]] or job["result"]["feeds_mm"] == [(0.0, 0.0)]


def test_natural_and_rhcp_beams_give_different_control_maps_via_api(client, grounded_library):
    upload_id = _upload_example_field(client)

    def _run(pol):
        body = {
            "upload_id": upload_id,
            "formulation": "B",
            "frequency_hz": 10.54e9,
            "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0, "desired_polarization": pol}],
            "period_mm": 8.0,
            "library_id": grounded_library,
        }
        resp = client.post("/api/hmsa/synthesize", json=body)
        return _wait_for_job(client, resp.json()["id"], timeout=120.0)

    natural = _run("natural")
    rhcp = _run("rhcp")
    assert natural["status"] == "done" and rhcp["status"] == "done"
    assert natural["result"]["control_map"] != rhcp["result"]["control_map"]


def test_explicit_feed_position_override(client, grounded_library):
    upload_id = _upload_example_field(client)
    body = {
        "upload_id": upload_id,
        "formulation": "B",
        "frequency_hz": 10.54e9,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        "period_mm": 8.0,
        "library_id": grounded_library,
        "feed_x_mm": 5.0,
        "feed_y_mm": -3.0,
    }
    resp = client.post("/api/hmsa/synthesize", json=body)
    job = _wait_for_job(client, resp.json()["id"], timeout=120.0)
    assert job["status"] == "done", job
    assert not any(label.startswith("FEED_POSITION_AUTO_DETECTED(") for label in job["result"]["labels"])
    assert job["result"]["feeds_mm"][0] == [5.0, -3.0] or job["result"]["feeds_mm"][0] == (5.0, -3.0)


def test_sector_multiplexing_via_api(client, grounded_library):
    upload_id = _upload_example_field(client)

    def _run(multiplexing):
        body = {
            "upload_id": upload_id,
            "formulation": "B",
            "frequency_hz": 10.54e9,
            "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}, {"theta_deg": 20.0, "phi_deg": 180.0}],
            "period_mm": 8.0,
            "library_id": grounded_library,
            "multiplexing": multiplexing,
        }
        resp = client.post("/api/hmsa/synthesize", json=body)
        return _wait_for_job(client, resp.json()["id"], timeout=120.0)

    full = _run("full_aperture")
    sectors = _run("sectors")
    assert full["status"] == "done" and sectors["status"] == "done"
    assert full["result"]["control_map"] != sectors["result"]["control_map"]


# =========================================================================== Part B: vector pattern


def test_pattern_predict_vector_model_via_api(client, grounded_library):
    upload_id = _upload_example_field(client)
    body = {
        "upload_id": upload_id,
        "formulation": "B",
        "frequency_hz": 10.54e9,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0, "desired_polarization": "rhcp"}],
        "period_mm": 8.0,
        "library_id": grounded_library,
    }
    synth_resp = client.post("/api/hmsa/synthesize", json=body)
    synth_job = _wait_for_job(client, synth_resp.json()["id"], timeout=120.0)
    assert synth_job["status"] == "done", synth_job

    pattern_resp = client.post(
        "/api/pattern/predict",
        json={
            "source_job_id": synth_resp.json()["id"],
            "frequency_hz": 10.54e9,
            "n_theta": 13,
            "n_phi": 9,
            "beams": [{"theta_deg": 20.0, "phi_deg": 0.0, "desired_polarization": "rhcp", "oam_charge": 0}],
            "pattern_model": "vector",
            "eps_r": 2.2,
            "tan_delta": 0.0009,
            "thickness_mm": 1.575,
        },
    )
    pattern_job = _wait_for_job(client, pattern_resp.json()["id"], timeout=60.0)
    assert pattern_job["status"] == "done", pattern_job
    result = pattern_job["result"]
    assert "VECTOR_PATTERN_MODEL" in result["labels"]
    assert len(result["per_beam"]) == 1
    entry = result["per_beam"][0]
    assert entry["desired_polarization"] == "rhcp"
    assert entry["beam_type_requested"] == "Pencil"
    assert entry["xpd_db"] is not None
    assert entry["phase_ring_png_url"]
    png_resp = client.get(entry["phase_ring_png_url"])
    assert png_resp.status_code == 200 and png_resp.content[:8] == b"\x89PNG\r\n\x1a\n"


def test_pattern_predict_defaults_to_scalar_and_is_unaffected(client, grounded_library):
    """Backward-compat: a request that omits pattern_model must behave exactly as before Phase 15."""
    upload_id = _upload_example_field(client)
    body = {
        "upload_id": upload_id,
        "formulation": "B",
        "frequency_hz": 10.54e9,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        "period_mm": 8.0,
        "library_id": grounded_library,
    }
    synth_resp = client.post("/api/hmsa/synthesize", json=body)
    synth_job = _wait_for_job(client, synth_resp.json()["id"], timeout=120.0)
    assert synth_job["status"] == "done", synth_job

    pattern_resp = client.post(
        "/api/pattern/predict",
        json={
            "source_job_id": synth_resp.json()["id"],
            "frequency_hz": 10.54e9,
            "n_theta": 13,
            "n_phi": 9,
            "beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        },
    )
    pattern_job = _wait_for_job(client, pattern_resp.json()["id"], timeout=60.0)
    assert pattern_job["status"] == "done", pattern_job
    result = pattern_job["result"]
    assert "SCALAR_PATTERN_MODEL_LEGACY" in result["labels"]
    assert result["per_beam"] == []
    assert result["result"]["main_beam_theta_deg"] is not None


def test_pattern_predict_vector_model_requires_substrate_params(client, grounded_library):
    upload_id = _upload_example_field(client)
    body = {
        "upload_id": upload_id,
        "formulation": "B",
        "frequency_hz": 10.54e9,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        "period_mm": 8.0,
        "library_id": grounded_library,
    }
    synth_resp = client.post("/api/hmsa/synthesize", json=body)
    synth_job = _wait_for_job(client, synth_resp.json()["id"], timeout=120.0)
    assert synth_job["status"] == "done", synth_job

    pattern_resp = client.post(
        "/api/pattern/predict",
        json={
            "source_job_id": synth_resp.json()["id"],
            "frequency_hz": 10.54e9,
            "n_theta": 13,
            "n_phi": 9,
            "beams": [],
            "pattern_model": "vector",
        },
    )
    pattern_job = _wait_for_job(client, pattern_resp.json()["id"], timeout=60.0)
    assert pattern_job["status"] == "failed", pattern_job


# =============================================================================== Part C: multi-feed


def test_multi_feed_synthesize_analytic_tm0_via_api(client, grounded_library):
    upload_id = _upload_example_field(client)
    body = {
        "formulation": "B",
        "frequency_hz": 10.54e9,
        "object_beams": [
            {"theta_deg": 20.0, "phi_deg": 0.0, "feed_index": 0},
            {"theta_deg": 15.0, "phi_deg": 180.0, "feed_index": 1},
        ],
        "period_mm": 8.0,
        "library_id": grounded_library,
        "feeds": [
            {"x_mm": -15.0, "y_mm": 0.0, "upload_id": upload_id},
            {"x_mm": 15.0, "y_mm": 0.0},
        ],
        "multi_feed_strategy": "superposition",
        "multi_feed_eps_r": 2.2,
        "multi_feed_thickness_mm": 1.575,
    }
    resp = client.post("/api/hmsa/synthesize", json=body)
    job = _wait_for_job(client, resp.json()["id"], timeout=120.0)
    assert job["status"] == "done", job
    assert any(label.startswith("MULTI_FEED_SUPERPOSITION(") for label in job["result"]["labels"])
    assert len(job["result"]["feeds_mm"]) == 2


def test_multi_feed_without_any_upload_needs_substrate_for_analytic_feeds(client, grounded_library):
    body = {
        "formulation": "B",
        "frequency_hz": 10.54e9,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        "period_mm": 8.0,
        "library_id": grounded_library,
        "feeds": [{"x_mm": -15.0, "y_mm": 0.0}, {"x_mm": 15.0, "y_mm": 0.0}],
    }
    resp = client.post("/api/hmsa/synthesize", json=body)
    job = _wait_for_job(client, resp.json()["id"], timeout=120.0)
    assert job["status"] == "failed", job


# =================================================================================== Export feeds


def test_export_marks_feed_locations(client, grounded_library):
    upload_id = _upload_example_field(client)
    body = {
        "upload_id": upload_id,
        "formulation": "B",
        "frequency_hz": 10.54e9,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        "period_mm": 8.0,
        "library_id": grounded_library,
        "feed_x_mm": 4.0,
        "feed_y_mm": -2.0,
    }
    synth_resp = client.post("/api/hmsa/synthesize", json=body)
    synth_job = _wait_for_job(client, synth_resp.json()["id"], timeout=120.0)
    assert synth_job["status"] == "done", synth_job

    export_resp = client.post(
        "/api/export/array",
        json={
            "source_job_id": synth_resp.json()["id"],
            "period_mm": 8.0,
            "eps_r": 2.2,
            "tan_delta": 0.0009,
            "thickness_mm": 1.575,
            "include_step": False,
        },
    )
    export_job = _wait_for_job(client, export_resp.json()["id"], timeout=60.0)
    assert export_job["status"] == "done", export_job
    readme_resp = client.get(export_job["result"]["readme_download_url"])
    assert "Feed 0: (4.000, -2.000) mm" in readme_resp.text
