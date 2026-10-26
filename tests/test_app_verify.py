"""Phase 13 Part 3: the Verify page's API -- cost estimate (including the refused/infinite-cost
case, which must come back as valid JSON, not a bare `Infinity` token) and a real, small
full-supercell MoM run. Needs the optional 'app' extra; skipped entirely otherwise."""

import time

import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from src.app.main import create_app  # noqa: E402


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("LWMOM_APP_DATA_DIR", str(tmp_path / "lwmom_data"))
    return TestClient(create_app())


def _wait_for_job(client, job_id, timeout=60.0):
    deadline = time.monotonic() + timeout
    while True:
        resp = client.get(f"/api/jobs/{job_id}")
        job = resp.json()
        if job["status"] in ("done", "failed", "cancelled"):
            return job
        if time.monotonic() > deadline:
            raise TimeoutError(f"job {job_id} did not finish in {timeout}s: {job}")
        time.sleep(0.05)


def test_cost_estimate_at_a_measured_point(client):
    resp = client.get("/api/verify/cost_estimate", params={"n_cells": 6})
    assert resp.status_code == 200
    body = resp.json()
    assert body["is_measured"] is True
    assert body["refused"] is False
    assert body["estimated_runtime_s"] == pytest.approx(6.997816)


def test_cost_estimate_refused_is_valid_json_not_bare_infinity(client):
    resp = client.get("/api/verify/cost_estimate", params={"n_cells": 25})
    assert resp.status_code == 200
    body = resp.json()  # would raise if the body were literal `Infinity` (invalid JSON)
    assert body["refused"] is True
    assert body["estimated_runtime_s"] is None
    assert body["refusal_reason"] is not None


def test_run_refuses_a_known_infeasible_size_without_force(client):
    resp = client.post("/api/verify/run", json={"n_cells": 25, "force": False})
    assert resp.status_code == 409


@pytest.mark.slow
def test_run_verification_at_the_smallest_measured_size(client):
    resp = client.post("/api/verify/run", json={"n_cells": 6})
    job = _wait_for_job(client, resp.json()["id"], timeout=30.0)
    assert job["status"] == "done", job
    result = job["result"]
    assert result["n_cells"] == 6
    assert result["cost"]["is_measured"] is True
    assert "predicted_vs_measured_intended_fraction_diff" in result
    # PHASE13B A10: the app's summary line reads this field directly -- must be a real number,
    # not missing/None (the blank "(predicted , |...|" bug).
    assert isinstance(result["predicted_intended_fraction"], float)
    assert result["aperture_result"]["intended_order"] is not None
