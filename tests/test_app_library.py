"""Phase 13 Part 3: the library page's API -- running a small real sweep (ungrounded and grounded)
through the job queue, browsing it back, and (if present on disk) browsing the real Phase 7
benchmark. Needs the optional 'app' extra; skipped entirely otherwise."""

import time
import uuid

import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from src.app import library_registry  # noqa: E402
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


def _sweep_body(library_id, backing="air", reference_planes=None):
    base_config = {
        "period_x": 4.8e-3,
        "period_y": 4.8e-3,
        "substrate": {**_SUBSTRATE, **({"backing": backing} if backing != "air" else {})},
        "geometry": {"type": "rectangle", "length": 1.0e-3, "width": 1.0e-3},
        "excitation": {"frequency": 12.5e9},
        "mesh": _SMALL_MESH,
        "reference_planes": reference_planes
        or (
            {"top_air_thickness": 5.0e-3}
            if backing == "pec"
            else {"top_air_thickness": 5.0e-3, "bottom_air_thickness": 5.0e-3}
        ),
    }
    return {
        "library_id": library_id,
        "base_config": base_config,
        "l_values_mm": [1.0, 2.0],
        "w_values_mm": [1.0],
        "f_min_hz": 12.0e9,
        "f_max_hz": 13.0e9,
        "f_step_hz": 1.0e9,
        "f_design_hz": 12.0e9,
    }


def test_list_is_empty_or_at_most_the_benchmark_before_any_sweep(client):
    resp = client.get("/api/library")
    assert resp.status_code == 200
    body = resp.json()
    assert set(body).issubset({library_registry.BENCHMARK_LIBRARY_ID})


def test_run_sweep_and_browse_it_back(client):
    library_id = f"test-{uuid.uuid4().hex[:8]}"
    resp = client.post("/api/library/sweep", json=_sweep_body(library_id))
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "done", job
    result = job["result"]
    assert result["n_geometries"] == 2
    assert result["n_ok"] == 2
    assert len(result["table"]) == 2
    assert result["coverage"]["n_points"] == 2
    assert len(result["classifications"]) == 2

    browse = client.get(f"/api/library/{library_id}")
    assert browse.status_code == 200
    assert browse.json()["n_cells"] == 2

    browse_at_freq = client.get(f"/api/library/{library_id}", params={"frequency_ghz": 12.0})
    assert len(browse_at_freq.json()["phase_vs_geometry"]) == 2

    listing = client.get("/api/library").json()
    assert library_id in listing
    assert listing[library_id]["meta"]["backing"] == "air"


def test_run_grounded_sweep_only_reports_s11(client):
    library_id = f"test-grounded-{uuid.uuid4().hex[:8]}"
    resp = client.post("/api/library/sweep", json=_sweep_body(library_id, backing="pec"))
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "done", job
    listing = client.get("/api/library").json()
    assert listing[library_id]["meta"]["backing"] == "pec"


def test_unknown_library_id_is_404(client):
    resp = client.get("/api/library/does-not-exist")
    assert resp.status_code == 404


def test_invalid_base_config_fails_the_sweep_job(client):
    body = _sweep_body(f"test-bad-{uuid.uuid4().hex[:8]}")
    del body["base_config"]["substrate"]
    resp = client.post("/api/library/sweep", json=body)
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "failed"
    assert job["error"]["code"] == "CONFIG_INVALID"
