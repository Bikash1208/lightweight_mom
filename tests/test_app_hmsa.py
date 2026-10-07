"""Phase 13 Part 3: the HMSA designer page's API -- upload the real example .fld file, build a
small real grounded library, and synthesize Formulation B end to end. Formulation A is exercised
with a synthetic dispersion CSV (no real HFSS eigenmode export is bundled with the repo). Needs
the optional 'app' extra; skipped entirely otherwise."""

import base64
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
    # 13C Step 1: "show immediately: a PNG of the reference phase map".
    assert body["phase_map_png_base64"]
    assert base64.b64decode(body["phase_map_png_base64"])[:8] == b"\x89PNG\r\n\x1a\n"


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
    # PHASE13B A2: resampled onto the 8mm-period cell grid (300mm aperture / 8mm -> 37x37), not
    # one "cell" per raw 2mm field sample (151x151).
    assert len(result["cells"]) == 37 * 37
    assert result["cells"][0]["x_m"] != pytest.approx(-0.15)  # not the raw field's own first sample
    assert "LOCAL_PERIODICITY_APPROXIMATION" in result["labels"]
    assert "aperture_field_real" in result
    # The example file's single feed-point NaN must be filled (so synthesis can run at all) and
    # reported, never silently dropped.
    assert "FEED_REGION_MASKED" in result["labels"]
    assert len(result["masked_cells"]) == 1

    # 13C Step 3: "the hologram/layout preview" PNG.
    layout_resp = client.get(result["layout_png_url"])
    assert layout_resp.status_code == 200
    assert layout_resp.content[:8] == b"\x89PNG\r\n\x1a\n"

    # PHASE13B B2: Export was not wired up from an HMSA design -- it must chain exactly like a
    # reflectarray design job (real L/W per cell, not a placeholder).
    export_resp = client.post(
        "/api/export/array",
        json={
            "source_job_id": resp.json()["id"],
            "period_mm": 8.0,
            "eps_r": 2.2,
            "tan_delta": 0.0009,
            "thickness_mm": 1.575,
            "include_step": False,  # avoid the optional 'cad' extra just to prove the wiring
        },
    )
    export_job = _wait_for_job(client, export_resp.json()["id"], timeout=60.0)
    assert export_job["status"] == "done", export_job
    assert export_job["result"]["n_cells_per_side"] == 37
    assert export_job["result"]["total_cells"] == 37 * 37
    assert export_job["result"]["layout_csv_download_url"]

    # Chain into Pattern prediction and confirm its own plots are produced too (13C Step 3).
    pattern_resp = client.post(
        "/api/pattern/predict",
        json={
            "source_job_id": resp.json()["id"],
            "frequency_hz": 10.54e9,
            "n_theta": 15,
            "n_phi": 9,
            "beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        },
    )
    pattern_job = _wait_for_job(client, pattern_resp.json()["id"], timeout=60.0)
    assert pattern_job["status"] == "done", pattern_job
    pattern_result = pattern_job["result"]
    assert pattern_result["pattern_uv_png_url"]
    assert len(pattern_result["pattern_cut_png_urls"]) == 1
    uv_resp = client.get(pattern_result["pattern_uv_png_url"])
    assert uv_resp.status_code == 200 and uv_resp.content[:8] == b"\x89PNG\r\n\x1a\n"


def test_synthesize_at_a_frequency_the_library_never_solved_is_a_named_error(client, grounded_library):
    """PHASE13B A3(a): the grounded_library fixture only has a record at 10.54 GHz; requesting
    12.5 GHz must fail with LIBRARY_FREQUENCY_UNAVAILABLE naming what IS available, not silently
    synthesize against zero candidates."""
    with open(_EXAMPLE_FLD, "rb") as fh:
        upload_resp = client.post(
            "/api/hmsa/upload_field", files={"phase_file": ("efl_phase_EzZ_example.fld", fh, "text/plain")}
        )
    upload_id = upload_resp.json()["upload_id"]
    body = {
        "upload_id": upload_id,
        "formulation": "B",
        "frequency_hz": 12.5e9,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        "period_mm": 8.0,
        "library_id": grounded_library,
    }
    resp = client.post("/api/hmsa/synthesize", json=body)
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "failed", job
    assert job["error"]["code"] == "LIBRARY_FREQUENCY_UNAVAILABLE"
    assert "10.54" in job["error"]["message"]


def test_synthesize_frequency_mismatch_with_the_upload_is_labeled_not_silent(client, grounded_library):
    """PHASE13B A3(b): the reviewer uploaded a 10.54 GHz field and synthesized at a different
    frequency with no explanation of the mismatch. Upload at 10.0 GHz, synthesize at 10.54 GHz
    (a frequency the library DOES have) so the run succeeds but must carry a loud label."""
    with open(_EXAMPLE_FLD, "rb") as fh:
        upload_resp = client.post(
            "/api/hmsa/upload_field",
            files={"phase_file": ("efl_phase_EzZ_example.fld", fh, "text/plain")},
            data={"eps_r": "2.2", "thickness_mm": "1.575", "frequency_hz": str(10.0e9)},
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
    assert any(label.startswith("REFERENCE_FIELD_FREQUENCY_MISMATCH(") for label in job["result"]["labels"])


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
    # PHASE13B B4: Formulation A's kernel needs a complex field, but the example file (like every
    # real reference-field export this project has) is phase-only -- it must run anyway with a
    # loudly-labeled uniform-amplitude assumption, not refuse outright.
    assert job["status"] == "done", job
    result = job["result"]
    assert "REFERENCE_AMPLITUDE_ASSUMED_UNIFORM" in result["labels"]
    # This dispersion table has no beta_rad_per_m column, so the B3 beta_bar rescale is silently
    # unavailable (no beta to rescale to) -- neither rescale label should appear.
    assert not any(label.startswith("REFERENCE_WAVE_") for label in result["labels"])


def test_formulation_a_rescales_to_beta_bar_by_default_when_the_table_has_it(client, tmp_path):
    """PHASE13B B3: with a dispersion table that DOES carry beta_rad_per_m, the default
    (rescale_to_beta_bar=True) must apply the correction and label it; opting out must instead
    report the accumulated-error estimate without changing the result."""
    csv_path = tmp_path / "dispersion.csv"
    csv_path.write_text(
        "geometry_param,beta_rad_per_m,x_s_ohm\n"
        "2.0,230.0,50.0\n2.5,228.0,80.0\n3.0,225.0,110.0\n3.5,222.0,140.0\n4.0,219.0,170.0\n"
    )
    with open(csv_path, "rb") as fh:
        disp_resp = client.post(
            "/api/hmsa/upload_dispersion",
            files={"csv_file": ("dispersion.csv", fh, "text/csv")},
            data={"frequency_hz": "10.54e9"},
        )
    dispersion_id = disp_resp.json()["dispersion_id"]

    with open(_EXAMPLE_FLD, "rb") as fh:
        upload_resp = client.post(
            "/api/hmsa/upload_field", files={"phase_file": ("efl_phase_EzZ_example.fld", fh, "text/plain")}
        )
    upload_id = upload_resp.json()["upload_id"]

    base_body = {
        "upload_id": upload_id,
        "formulation": "A",
        "frequency_hz": 10.54e9,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        "period_mm": 8.0,
        "dispersion_id": dispersion_id,
        "x_avg_ohm": 100.0,
        "modulation_index": 0.3,
    }

    resp_default = client.post("/api/hmsa/synthesize", json=base_body)
    job_default = _wait_for_job(client, resp_default.json()["id"], timeout=120.0)
    assert job_default["status"] == "done", job_default
    labels_default = job_default["result"]["labels"]
    assert any(label.startswith("REFERENCE_WAVE_RESCALED_TO_BETA_BAR(") for label in labels_default)

    resp_opt_out = client.post("/api/hmsa/synthesize", json={**base_body, "rescale_to_beta_bar": False})
    job_opt_out = _wait_for_job(client, resp_opt_out.json()["id"], timeout=120.0)
    assert job_opt_out["status"] == "done", job_opt_out
    labels_opt_out = job_opt_out["result"]["labels"]
    assert any(label.startswith("REFERENCE_WAVE_NOT_RESCALED(") for label in labels_opt_out)


@pytest.mark.slow
def test_ensure_dispersion_builds_reuses_and_feeds_formulation_a(client):
    """PHASE14 Part B: the "Surface-wave (eigenmode)" method's own auto-built dispersion table
    (dispersion.sheet_model) -- built once, reused on repeat, and directly usable by Formulation A
    with no separate upload step, using the X_avg/modulation depth it derives automatically."""
    body = {"period_mm": 4.0, "eps_r": 2.2, "tan_delta": 0.0009, "thickness_mm": 1.575, "frequency_hz": 10.54e9}
    resp1 = client.post("/api/hmsa/ensure_dispersion", json=body)
    job1 = _wait_for_job(client, resp1.json()["id"], timeout=120.0)
    assert job1["status"] == "done", job1
    result1 = job1["result"]
    assert result1["already_built"] is False
    assert result1["n_geometries"] > 0
    assert result1["x_avg_ohm"] is not None and result1["modulation_index"] is not None
    assert 0.0 < result1["modulation_index"] < 1.0
    dispersion_id = result1["dispersion_id"]

    resp2 = client.post("/api/hmsa/ensure_dispersion", json=body)
    job2 = _wait_for_job(client, resp2.json()["id"], timeout=30.0)
    assert job2["status"] == "done", job2
    assert job2["result"]["already_built"] is True
    assert job2["result"]["dispersion_id"] == dispersion_id

    with open(_EXAMPLE_FLD, "rb") as fh:
        upload_resp = client.post(
            "/api/hmsa/upload_field", files={"phase_file": ("efl_phase_EzZ_example.fld", fh, "text/plain")}
        )
    upload_id = upload_resp.json()["upload_id"]

    synth_body = {
        "upload_id": upload_id,
        "formulation": "A",
        "frequency_hz": 10.54e9,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        "period_mm": 4.0,
        "dispersion_id": dispersion_id,
        "x_avg_ohm": result1["x_avg_ohm"],
        "modulation_index": result1["modulation_index"],
    }
    synth_resp = client.post("/api/hmsa/synthesize", json=synth_body)
    synth_job = _wait_for_job(client, synth_resp.json()["id"], timeout=120.0)
    assert synth_job["status"] == "done", synth_job
    assert (
        "DISPERSION_FROM_SHEET_MODEL" not in synth_job["result"]["labels"]
    )  # not a design label; it's the table's own provenance, not carried onto the synth result -- confirms no accidental mislabeling
    assert any(label.startswith("REFERENCE_WAVE_RESCALED_TO_BETA_BAR(") for label in synth_job["result"]["labels"])


@pytest.mark.slow
def test_ensure_library_builds_once_and_reuses_on_repeat(client):
    """PHASE13C: the wizard never sends the user to a separate Library page -- /ensure_library
    must build a real square-patch library for a chosen board/frequency the first time, then
    reuse it (not re-solve) on a repeat request for the identical board."""
    body = {"period_mm": 8.0, "eps_r": 2.2, "tan_delta": 0.0009, "thickness_mm": 1.575, "frequency_hz": 12.5e9}
    resp1 = client.post("/api/hmsa/ensure_library", json=body)
    job1 = _wait_for_job(client, resp1.json()["id"], timeout=120.0)
    assert job1["status"] == "done", job1
    assert job1["result"]["already_built"] is False
    assert job1["result"]["n_geometries"] == 13  # PHASE13D: denser, mesh-aligned auto-library
    library_id = job1["result"]["library_id"]

    lib_resp = client.get(f"/api/library/{library_id}")
    assert lib_resp.status_code == 200
    assert lib_resp.json()["n_cells"] == 13

    resp2 = client.post("/api/hmsa/ensure_library", json=body)
    job2 = _wait_for_job(client, resp2.json()["id"], timeout=30.0)
    assert job2["status"] == "done", job2
    assert job2["result"]["already_built"] is True
    assert job2["result"]["library_id"] == library_id
    # PHASE13D: the achievable phase-lookup range is reported (both the fresh-build and the
    # already-built fast paths), so Step 2 can show it before synthesis runs.
    assert job1["result"]["phase_span_deg"] > 0
    assert job2["result"]["phase_span_deg"] == job1["result"]["phase_span_deg"]


@pytest.mark.slow
def test_ensure_library_accurate_tier_builds_a_separate_finer_library(client):
    """PHASE13D item 2 ("Words"): the Advanced "Accuracy" control -- "accurate" must be a
    genuinely different (finer-mesh) library, not just a relabeled fast one."""
    body_fast = {"period_mm": 8.0, "eps_r": 2.2, "tan_delta": 0.0009, "thickness_mm": 1.575, "frequency_hz": 12.5e9}
    body_accurate = {**body_fast, "accuracy": "accurate"}
    resp_fast = client.post("/api/hmsa/ensure_library", json=body_fast)
    job_fast = _wait_for_job(client, resp_fast.json()["id"], timeout=120.0)
    resp_accurate = client.post("/api/hmsa/ensure_library", json=body_accurate)
    job_accurate = _wait_for_job(client, resp_accurate.json()["id"], timeout=120.0)
    assert job_fast["status"] == "done" and job_accurate["status"] == "done"
    assert job_fast["result"]["library_id"] != job_accurate["result"]["library_id"]
