"""Phase 13 Part 3: the HMSA designer page's API -- upload the real example .fld file, build a
small real grounded library, and synthesize Formulation B end to end. Formulation A is exercised
with a synthetic dispersion CSV (no real HFSS eigenmode export is bundled with the repo). Needs
the optional 'app' extra; skipped entirely otherwise."""

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
    library_id = f"test-hmsa-{uuid.uuid4().hex[:8]}"
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


def test_upload_field_parses_the_real_example_file(client):
    with open(_EXAMPLE_FLD, "rb") as fh:
        resp = client.post(
            "/api/hmsa/upload_field",
            files={"phase_file": ("efl_phase_EzZ_example.fld", fh, "text/plain")},
            data={"eps_r": "2.2", "thickness_mm": "1.575", "frequency_hz": "10.54e9"},
        )
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert body["nx"] == 151 and body["ny"] == 151
    assert body["nan_count"] == 1  # the single feed-point NaN
    assert body["regularity"]["regular"] is True
    assert body["tm0_comparison"] is not None
    assert body["upload_id"]


def test_upload_field_rejects_a_malformed_file(client, tmp_path):
    bad = tmp_path / "bad.fld"
    bad.write_text("not a valid fld file\nsecond line\n1 2 3\n")
    with open(bad, "rb") as fh:
        resp = client.post("/api/hmsa/upload_field", files={"phase_file": ("bad.fld", fh, "text/plain")})
    assert resp.status_code == 400


def test_synthesize_formulation_b_end_to_end(client, grounded_library):
    with open(_EXAMPLE_FLD, "rb") as fh:
        upload_resp = client.post(
            "/api/hmsa/upload_field", files={"phase_file": ("efl_phase_EzZ_example.fld", fh, "text/plain")}
        )
    upload_id = upload_resp.json()["upload_id"]

    body = {
        "upload_id": upload_id,
        "formulation": "B",
        "frequency_hz": 10.54e9,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        "period_mm": 8.0,
        "library_id": grounded_library,
    }
    resp = client.post("/api/hmsa/synthesize", json=body)
    job = _wait_for_job(client, resp.json()["id"], timeout=120.0)
    assert job["status"] == "done", job
    result = job["result"]
    assert result["formulation"] == "B_phase_difference"
    assert len(result["cells"]) == 151 * 151
    assert "LOCAL_PERIODICITY_APPROXIMATION" in result["labels"]
    assert "aperture_field_real" in result
    # The example file's single feed-point NaN must be filled (so synthesis can run at all) and
    # reported, never silently dropped.
    assert "FEED_REGION_MASKED" in result["labels"]
    assert len(result["masked_cells"]) == 1


def test_synthesize_unknown_library_is_404(client):
    with open(_EXAMPLE_FLD, "rb") as fh:
        upload_resp = client.post(
            "/api/hmsa/upload_field", files={"phase_file": ("efl_phase_EzZ_example.fld", fh, "text/plain")}
        )
    upload_id = upload_resp.json()["upload_id"]
    resp = client.post(
        "/api/hmsa/synthesize",
        json={
            "upload_id": upload_id,
            "formulation": "B",
            "frequency_hz": 10.54e9,
            "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
            "period_mm": 8.0,
            "library_id": "does-not-exist",
        },
    )
    assert resp.status_code == 404


def test_upload_dispersion_and_synthesize_formulation_a(client, tmp_path):
    csv_path = tmp_path / "dispersion.csv"
    csv_path.write_text("geometry_param,x_s_ohm\n2.0,50.0\n2.5,80.0\n3.0,110.0\n3.5,140.0\n4.0,170.0\n")
    with open(csv_path, "rb") as fh:
        disp_resp = client.post(
            "/api/hmsa/upload_dispersion",
            files={"csv_file": ("dispersion.csv", fh, "text/csv")},
            data={"frequency_hz": "10.54e9"},
        )
    assert disp_resp.status_code == 200, disp_resp.text
    dispersion_id = disp_resp.json()["dispersion_id"]
    assert disp_resp.json()["has_x_s"] is True

    with open(_EXAMPLE_FLD, "rb") as fh:
        upload_resp = client.post(
            "/api/hmsa/upload_field", files={"phase_file": ("efl_phase_EzZ_example.fld", fh, "text/plain")}
        )
    upload_id = upload_resp.json()["upload_id"]

    body = {
        "upload_id": upload_id,
        "formulation": "A",
        "frequency_hz": 10.54e9,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        "period_mm": 8.0,
        "dispersion_id": dispersion_id,
        "x_avg_ohm": 100.0,
        "modulation_index": 0.3,
    }
    resp = client.post("/api/hmsa/synthesize", json=body)
    job = _wait_for_job(client, resp.json()["id"], timeout=120.0)
    # Formulation A needs a COMPLEX reference field (amplitude available); the example file is
    # phase-only, so this is expected to fail with a clear, typed error -- not silently succeed
    # with a fabricated amplitude (hologram.hmsa_synthesis.formulation_a_control_map's own rule).
    assert job["status"] == "failed", job
