"""Phase 13 Part 3: reflectarray designer -> pattern prediction -> export, chained through the job
queue exactly as the app's own pages do (a small real grounded library, a small real array).
Needs the optional 'app'/'cad' extras; skipped entirely otherwise."""

import io
import json
import time
import uuid
import zipfile

import pytest

fastapi = pytest.importorskip("fastapi")
pytest.importorskip("cadquery")
from fastapi.testclient import TestClient  # noqa: E402

from src.app.main import create_app  # noqa: E402

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3}
_SMALL_MESH = {"cells_x": 8, "cells_y": 8, "floquet_order": 8}


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
    """Builds a small real grounded library (3x3 L/W grid, 1 frequency) via the actual sweep
    endpoint, exactly as the Library page would -- no shortcuts around the job queue."""
    library_id = f"test-refl-{uuid.uuid4().hex[:8]}"
    body = {
        "library_id": library_id,
        "base_config": {
            "period_x": 8.0e-3,
            "period_y": 8.0e-3,
            "substrate": {**_SUBSTRATE, "backing": "pec"},
            "geometry": {"type": "rectangle", "length": 2.0e-3, "width": 2.0e-3},
            "excitation": {"frequency": 12.5e9},
            "mesh": _SMALL_MESH,
            "reference_planes": {"top_air_thickness": 5.0e-3},
        },
        "l_values_mm": [2.0, 3.0, 4.0],
        "w_values_mm": [2.0, 3.0, 4.0],
        "f_min_hz": 12.5e9,
        "f_max_hz": 12.5e9,
        "f_step_hz": 1.0e9,
        "f_design_hz": 12.5e9,
    }
    resp = client.post("/api/library/sweep", json=body)
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "done", job
    return library_id


def test_reflectarray_design_selects_cells_and_reports_labels(client, grounded_library):
    resp = client.post(
        "/api/reflectarray/design",
        json={
            "library_id": grounded_library,
            "frequency_hz": 12.5e9,
            "nx": 3,
            "ny": 3,
            "period_mm": 8.0,
            "feed_xyz_mm": [0.0, 0.0, 50.0],
            "theta0_deg": 0.0,
            "phi0_deg": 0.0,
            "taper_q": 1.0,
        },
    )
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "done", job
    result = job["result"]
    assert len(result["assignments"]) == 9
    assert "LOCAL_PERIODICITY_APPROXIMATION" in result["labels"]
    for a in result["assignments"]:
        assert a["chosen_l_mm"] is not None
    assert result["phase_error_stats"]["n_cells"] == 9


def test_reflectarray_design_at_unavailable_frequency_is_a_named_error(client, grounded_library):
    """PHASE13B A3(a): the grounded_library fixture only has a record at 12.5 GHz."""
    resp = client.post(
        "/api/reflectarray/design",
        json={
            "library_id": grounded_library,
            "frequency_hz": 10.0e9,
            "nx": 2,
            "ny": 2,
            "period_mm": 8.0,
            "feed_xyz_mm": [0.0, 0.0, 50.0],
            "theta0_deg": 0.0,
            "phi0_deg": 0.0,
            "taper_q": 1.0,
        },
    )
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "failed", job
    assert job["error"]["code"] == "LIBRARY_FREQUENCY_UNAVAILABLE"
    assert "12.5" in job["error"]["message"]


def _design_job_id(client, grounded_library):
    resp = client.post(
        "/api/reflectarray/design",
        json={
            "library_id": grounded_library,
            "frequency_hz": 12.5e9,
            "nx": 3,
            "ny": 3,
            "period_mm": 8.0,
            "feed_xyz_mm": [0.0, 0.0, 50.0],
            "theta0_deg": 0.0,
            "phi0_deg": 0.0,
            "taper_q": 1.0,
        },
    )
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "done", job
    return resp.json()["id"]


def test_pattern_prediction_chains_off_a_design_job(client, grounded_library):
    design_id = _design_job_id(client, grounded_library)
    resp = client.post(
        "/api/pattern/predict",
        json={"source_job_id": design_id, "frequency_hz": 12.5e9, "n_theta": 19, "n_phi": 13},
    )
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "done", job
    result = job["result"]["result"]
    assert "main_beam_theta_deg" in result
    assert "directivity_estimate_dbi" in result


def test_pattern_prediction_404s_on_unknown_source(client):
    resp = client.post("/api/pattern/predict", json={"source_job_id": "nope", "frequency_hz": 12.5e9})
    assert resp.status_code == 404


def test_export_produces_a_real_step_file(client, grounded_library):
    design_id = _design_job_id(client, grounded_library)
    resp = client.post(
        "/api/export/array",
        json={
            "source_job_id": design_id,
            "period_mm": 8.0,
            "eps_r": 2.2,
            "tan_delta": 0.0009,
            "thickness_mm": 1.575,
            "include_step": True,
        },
    )
    job = _wait_for_job(client, resp.json()["id"], timeout=60.0)
    assert job["status"] == "done", job
    result = job["result"]
    assert result["n_cells_per_side"] == 3
    assert result["total_cells"] == 9
    assert result["step_download_url"] is not None

    step_resp = client.get(result["step_download_url"])
    assert step_resp.status_code == 200
    assert len(step_resp.content) > 0

    csv_resp = client.get(result["layout_csv_download_url"])
    assert csv_resp.status_code == 200
    assert b"i,j,length_m,width_m" in csv_resp.content


def test_export_package_bundles_a_real_downloadable_zip(client, grounded_library):
    """13C Step 3: one "Download package" (zip) with the STEP/layout/README and the project
    record, optionally including the pattern plots too."""
    design_id = _design_job_id(client, grounded_library)
    export_resp = client.post(
        "/api/export/array",
        json={
            "source_job_id": design_id,
            "period_mm": 8.0,
            "eps_r": 2.2,
            "tan_delta": 0.0009,
            "thickness_mm": 1.575,
            "include_step": False,
        },
    )
    export_job_id = export_resp.json()["id"]
    export_job = _wait_for_job(client, export_resp.json()["id"], timeout=60.0)
    assert export_job["status"] == "done", export_job

    pattern_resp = client.post(
        "/api/pattern/predict",
        json={"source_job_id": design_id, "frequency_hz": 12.5e9, "n_theta": 15, "n_phi": 9},
    )
    pattern_job_id = pattern_resp.json()["id"]
    pattern_job = _wait_for_job(client, pattern_resp.json()["id"], timeout=60.0)
    assert pattern_job["status"] == "done", pattern_job

    pkg_resp = client.post(
        "/api/export/package", json={"export_job_id": export_job_id, "pattern_job_id": pattern_job_id}
    )
    pkg_job = _wait_for_job(client, pkg_resp.json()["id"], timeout=30.0)
    assert pkg_job["status"] == "done", pkg_job
    download_url = pkg_job["result"]["download_url"]

    zip_resp = client.get(download_url)
    assert zip_resp.status_code == 200

    with zipfile.ZipFile(io.BytesIO(zip_resp.content)) as zf:
        names = set(zf.namelist())
        assert "layout.csv" in names
        assert "layout.json" in names
        assert "HFSS_README.md" in names
        assert "project_record.json" in names
        assert "pattern_uv.png" in names
        record = json.loads(zf.read("project_record.json"))
        assert record["export"]["total_cells"] == 9
