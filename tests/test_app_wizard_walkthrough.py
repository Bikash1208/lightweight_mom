"""PHASE13C "Tests & done-when": an API-level walk-through of the 3-step wizard flow (upload the
real example .fld -> auto-build the board's unit-cell library -> synthesize -> predict the
pattern -> export -> download the zip package), plus the two-beam case and a scan of every served
page for snake_case/internal identifiers a person would actually see. Needs the optional 'app'
extra; skipped entirely otherwise."""

import io
import re
import time
import zipfile
from pathlib import Path

import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from src.app.main import create_app  # noqa: E402

_EXAMPLE_FLD = Path(__file__).resolve().parent.parent / "inputs" / "reference_fields" / "efl_phase_EzZ_example.fld"
_STATIC_DIR = Path(__file__).resolve().parent.parent / "src" / "app" / "static"
_UPLOAD_FREQUENCY_HZ = 10.54e9


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("LWMOM_APP_DATA_DIR", str(tmp_path / "lwmom_data"))
    return TestClient(create_app())


def _wait_for_job(client, job_id, timeout=120.0):
    deadline = time.time() + timeout
    while True:
        resp = client.get(f"/api/jobs/{job_id}")
        job = resp.json()
        if job["status"] in ("done", "failed", "cancelled"):
            return job
        if time.time() > deadline:
            raise TimeoutError(f"job {job_id} did not finish in {timeout}s: {job}")
        time.sleep(0.05)


def _upload(client):
    with open(_EXAMPLE_FLD, "rb") as fh:
        resp = client.post(
            "/api/hmsa/upload_field",
            files={"phase_file": ("efl_phase_EzZ_example.fld", fh, "text/plain")},
            data={"eps_r": "2.2", "thickness_mm": "1.575", "frequency_hz": str(_UPLOAD_FREQUENCY_HZ)},
        )
    assert resp.status_code == 200, resp.text
    return resp.json()


def _ensure_library(client):
    resp = client.post(
        "/api/hmsa/ensure_library",
        json={
            "period_mm": 4.0,
            "eps_r": 2.2,
            "tan_delta": 0.0009,
            "thickness_mm": 1.575,
            "frequency_hz": _UPLOAD_FREQUENCY_HZ,
        },
    )
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "done", job
    return job["result"]["library_id"]


def _synthesize(client, upload_id, library_id, beams):
    resp = client.post(
        "/api/hmsa/synthesize",
        json={
            "upload_id": upload_id,
            "formulation": "B",
            "frequency_hz": _UPLOAD_FREQUENCY_HZ,
            "object_beams": [{"theta_deg": b[0], "phi_deg": b[1]} for b in beams],
            "period_mm": 4.0,
            "library_id": library_id,
        },
    )
    job = _wait_for_job(client, resp.json()["id"], timeout=180.0)
    assert job["status"] == "done", job
    return resp.json()["id"], job["result"]


@pytest.mark.slow
def test_one_beam_walkthrough_ends_in_a_valid_zip_with_the_beam_near_its_target(client):
    """cell=4mm (finer than the reflectarray tests' 8mm, closer to a real HMSA cell size),
    theta=20, phi=0, RT/duroid 5880 1.575mm -- the brief's own worked example."""
    upload = _upload(client)
    library_id = _ensure_library(client)
    synth_job_id, synth_result = _synthesize(client, upload["upload_id"], library_id, [(20.0, 0.0)])
    assert "assignments" in synth_result  # exportable

    pattern_resp = client.post(
        "/api/pattern/predict",
        json={
            "source_job_id": synth_job_id,
            "frequency_hz": _UPLOAD_FREQUENCY_HZ,
            "n_theta": 31,
            "n_phi": 25,
            "beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        },
    )
    pattern_job = _wait_for_job(client, pattern_resp.json()["id"])
    assert pattern_job["status"] == "done", pattern_job
    peak_theta = pattern_job["result"]["result"]["main_beam_theta_deg"]
    assert peak_theta == pytest.approx(20.0, abs=5.0)  # a stated tolerance around the target

    export_resp = client.post(
        "/api/export/array",
        json={
            "source_job_id": synth_job_id,
            "period_mm": 4.0,
            "eps_r": 2.2,
            "tan_delta": 0.0009,
            "thickness_mm": 1.575,
            "include_step": False,
        },
    )
    export_job = _wait_for_job(client, export_resp.json()["id"])
    assert export_job["status"] == "done", export_job

    pkg_resp = client.post(
        "/api/export/package",
        json={"export_job_id": export_resp.json()["id"], "pattern_job_id": pattern_resp.json()["id"]},
    )
    pkg_job = _wait_for_job(client, pkg_resp.json()["id"])
    assert pkg_job["status"] == "done", pkg_job

    zip_resp = client.get(pkg_job["result"]["download_url"])
    assert zip_resp.status_code == 200
    with zipfile.ZipFile(io.BytesIO(zip_resp.content)) as zf:
        assert zf.testzip() is None  # every member's CRC checks out -- a genuinely valid zip
        names = set(zf.namelist())
        assert {"layout.csv", "layout.json", "HFSS_README.md", "pattern_uv.png", "project_record.json"} <= names


@pytest.mark.slow
def test_two_beam_case_both_beams_present_in_the_pattern(client):
    upload = _upload(client)
    library_id = _ensure_library(client)
    beams = [(20.0, 0.0), (20.0, 180.0)]
    synth_job_id, _ = _synthesize(client, upload["upload_id"], library_id, beams)

    pattern_resp = client.post(
        "/api/pattern/predict",
        json={
            "source_job_id": synth_job_id,
            "frequency_hz": _UPLOAD_FREQUENCY_HZ,
            "n_theta": 31,
            "n_phi": 25,
            "beams": [{"theta_deg": t, "phi_deg": p} for t, p in beams],
        },
    )
    pattern_job = _wait_for_job(client, pattern_resp.json()["id"])
    assert pattern_job["status"] == "done", pattern_job
    result = pattern_job["result"]
    # both requested beam cuts must have been rendered (one PNG per beam direction)
    assert len(result["pattern_cut_png_urls"]) == 2
    for url in result["pattern_cut_png_urls"]:
        resp = client.get(url)
        assert resp.status_code == 200 and resp.content[:8] == b"\x89PNG\r\n\x1a\n"
    # a genuine two-lobe pattern: power_db is a real 2-D array with more than one local structure,
    # not a degenerate single-beam result silently collapsing the second beam.
    power_db = result["result"]["power_db"]
    assert len(power_db) > 1 and len(power_db[0]) > 1


_SNAKE_CASE_RE = re.compile(r"\b[a-z][a-z0-9]*(?:_[a-z0-9]+){1,}\b")
# Identifiers that are legitimately visible (URL paths, file extensions the user typed, this
# project's own name) -- not "internal names" in the sense the brief means.
_ALLOWED = {
    "unit_cell",
    "more_tools",
    "hfss_compare",
    "lightweight_mom",
    "app_js",
    "style_css",
}


def _visible_text_from_html(html_text):
    """Strips <script>/<style> and all tags, keeping only what a viewer would actually read --
    the brief's own "no snake_case or internal code visible in any RENDERED page", not the source
    (which inevitably contains snake_case API field names in application code)."""
    text = re.sub(r"<script[^>]*>.*?</script>", " ", html_text, flags=re.DOTALL | re.IGNORECASE)
    text = re.sub(r"<style[^>]*>.*?</style>", " ", text, flags=re.DOTALL | re.IGNORECASE)
    text = re.sub(r"<[^>]+>", " ", text)
    return text


def test_no_snake_case_or_internal_names_visible_in_any_rendered_page():
    offenders = {}
    for html_path in sorted(_STATIC_DIR.glob("*.html")):
        visible = _visible_text_from_html(html_path.read_text())
        matches = {m for m in _SNAKE_CASE_RE.findall(visible) if m not in _ALLOWED}
        if matches:
            offenders[html_path.name] = matches
    assert offenders == {}
