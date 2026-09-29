"""Phase 13 Part 3: the unit-cell page's API -- real solves through the job queue (fast mesh and
Richardson-extrapolated), STEP export, and config-error propagation. Needs the optional 'app'
extra; skipped entirely otherwise."""

import time

import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from src.app.main import create_app  # noqa: E402

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3}
_SMALL_MESH = {"cells_x": 8, "cells_y": 8, "floquet_order": 8}


def _config(**overrides):
    cfg = {
        "period_x": 4.8e-3,
        "period_y": 4.8e-3,
        "substrate": _SUBSTRATE,
        "geometry": {"type": "rectangle", "length": 2.0e-3, "width": 2.0e-3},
        "excitation": {"frequency": 12.5e9},
        "mesh": _SMALL_MESH,
        "reference_planes": {"top_air_thickness": 5.0e-3, "bottom_air_thickness": 5.0e-3},
    }
    cfg.update(overrides)
    return cfg


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("LWMOM_APP_DATA_DIR", str(tmp_path / "lwmom_data"))
    return TestClient(create_app())


def _wait_for_job(client, job_id, timeout=30.0):
    deadline = time.time() + timeout
    while True:
        resp = client.get(f"/api/jobs/{job_id}")
        job = resp.json()
        if job["status"] in ("done", "failed", "cancelled"):
            return job
        if time.time() > deadline:
            raise TimeoutError(f"job {job_id} did not finish in {timeout}s: {job}")
        time.sleep(0.05)


def test_solve_fast_mode_returns_a_record_with_labels(client):
    resp = client.post("/api/unit_cell/solve", json={"config": _config(), "mesh_mode": "fast"})
    assert resp.status_code == 200
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "done", job
    point = job["result"]["points"][0]
    assert point["record"]["frequency_hz"] == 12.5e9
    assert "MESH_NOT_CONVERGED" in point["labels"]
    assert "ZERO_THICKNESS_METAL_MODEL" in point["labels"]
    assert point["diagnostics"]["n_unknowns"] > 0


def test_solve_richardson_mode_reports_convergence_fit(client):
    resp = client.post("/api/unit_cell/solve", json={"config": _config(), "mesh_mode": "richardson"})
    job = _wait_for_job(client, resp.json()["id"], timeout=60.0)
    assert job["status"] == "done", job
    point = job["result"]["points"][0]
    assert point["convergence"]["phase"]["status"] in (
        "CONVERGED",
        "NOT_MONOTONIC_NO_FIT",
        "FIT_FAILED",
        "TOO_FEW_POINTS",
    )
    assert any(label in ("MESH_CONVERGED", "MESH_NOT_CONVERGED") for label in point["labels"])


def test_solve_sweeps_multiple_frequencies(client):
    resp = client.post(
        "/api/unit_cell/solve",
        json={"config": _config(), "frequencies_hz": [11.0e9, 12.0e9, 13.0e9], "mesh_mode": "fast"},
    )
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "done", job
    freqs = [p["record"]["frequency_hz"] for p in job["result"]["points"]]
    assert freqs == [11.0e9, 12.0e9, 13.0e9]


def test_solve_grounded_config_omits_two_port_fields(client):
    grounded = _config(
        substrate={**_SUBSTRATE, "backing": "pec"},
        reference_planes={"top_air_thickness": 5.0e-3},
    )
    resp = client.post("/api/unit_cell/solve", json={"config": grounded})
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "done", job
    record = job["result"]["points"][0]["record"]
    assert any(k.startswith("S11_") for k in record)
    assert not any(k.startswith("S21_") for k in record)


def test_invalid_config_fails_with_config_invalid_error_code(client):
    bad = _config()
    del bad["substrate"]
    resp = client.post("/api/unit_cell/solve", json={"config": bad})
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "failed"
    assert job["error"]["code"] == "CONFIG_INVALID"


def test_same_request_hits_the_job_cache(client):
    body = {"config": _config(), "mesh_mode": "fast"}
    resp1 = client.post("/api/unit_cell/solve", json=body)
    id1 = resp1.json()["id"]
    _wait_for_job(client, id1)
    resp2 = client.post("/api/unit_cell/solve", json=body)
    assert resp2.json()["id"] == id1


def test_export_step_produces_a_downloadable_file(client):
    pytest.importorskip("cadquery")  # STEP export needs the optional 'cad' extra
    resp = client.post("/api/unit_cell/export_step", json={"config": _config()})
    job = _wait_for_job(client, resp.json()["id"], timeout=30.0)
    assert job["status"] == "done", job
    url = job["result"]["download_url"]
    file_resp = client.get(url)
    assert file_resp.status_code == 200
    assert file_resp.content[:1] != b""  # a real STEP file was written and served
