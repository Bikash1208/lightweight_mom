"""Phase 13 Part 3: the FastAPI app factory -- health endpoint, static frontend, job queue wiring.
Needs the optional 'app' extra (fastapi/uvicorn); skipped entirely otherwise."""

import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from src.app.main import create_app  # noqa: E402


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("LWMOM_APP_DATA_DIR", str(tmp_path / "lwmom_data"))
    return TestClient(create_app())


def test_health_endpoint_reports_ok(client):
    resp = client.get("/api/health")
    assert resp.status_code == 200
    assert resp.json() == {"status": "ok"}


def test_static_index_is_served_at_root(client):
    resp = client.get("/")
    assert resp.status_code == 200
    assert "lwmom" in resp.text.lower()


def test_app_factory_attaches_a_job_queue():
    app = create_app()
    from src.app.jobs import JobQueue

    assert isinstance(app.state.jobs, JobQueue)


def test_unknown_job_id_returns_404(client):
    resp = client.get("/api/jobs/does-not-exist")
    assert resp.status_code == 404
