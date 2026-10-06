"""Phase 15 Part C gap-fill / Part E: the app's API surface for multi-port sub-apertures and the
Coded (1-bit/N-bit) design method -- end to end through the real HTTP app, matching
`test_app_phase15_polarization_vortex.py`'s own fixtures/conventions."""

import time
import uuid

import numpy as np
import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from src.app.main import create_app  # noqa: E402

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3}
_SMALL_MESH = {"cells_x": 8, "cells_y": 8, "floquet_order": 8}
_FREQ_HZ = 10.54e9


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
    library_id = f"test-phase15-ports-{uuid.uuid4().hex[:8]}"
    body = {
        "library_id": library_id,
        "base_config": {
            "period_x": 8.0e-3, "period_y": 8.0e-3,
            "substrate": {**_SUBSTRATE, "backing": "pec"},
            "geometry": {"type": "rectangle", "length": 2.0e-3, "width": 2.0e-3},
            "excitation": {"frequency": _FREQ_HZ},
            "mesh": _SMALL_MESH,
            "reference_planes": {"top_air_thickness": 5.0e-3},
        },
        "l_values_mm": [2.0, 3.0, 4.0, 5.0, 6.0],
        "w_values_mm": [2.0, 3.0, 4.0, 5.0, 6.0],
        "f_min_hz": _FREQ_HZ, "f_max_hz": _FREQ_HZ, "f_step_hz": 1.0e9, "f_design_hz": _FREQ_HZ,
    }
    resp = client.post("/api/library/sweep", json=body)
    job = _wait_for_job(client, resp.json()["id"])
    assert job["status"] == "done", job
    return library_id


# =================================================================================== Ports (Part C)


def _synthesize_ports_1xn(client, grounded_library):
    body = {
        "formulation": "B",
        "frequency_hz": _FREQ_HZ,
        "period_mm": 8.0,
        "layout": {"preset": "1xn", "n": 2, "cell_side_mm": 40.0, "gap_mm": 0.0},
        "ports": [
            {"beams": [{"theta_deg": 20.0, "phi_deg": 0.0}]},
            {"beams": [{"theta_deg": 20.0, "phi_deg": 180.0}]},
        ],
        "library_id": grounded_library,
        "port_eps_r": _SUBSTRATE["eps_r"],
        "port_thickness_mm": _SUBSTRATE["thickness"] * 1e3,
    }
    resp = client.post("/api/hmsa/synthesize_ports", json=body)
    return _wait_for_job(client, resp.json()["id"], timeout=60.0)


def test_synthesize_ports_1xn_layout_end_to_end(client, grounded_library):
    job = _synthesize_ports_1xn(client, grounded_library)
    assert job["status"] == "done", job
    result = job["result"]
    assert len(result["ports"]) == 2
    assert "MULTI_PORT_SHARED_METHOD" in result["labels"]
    assert "PORT_ISOLATION_NOT_FULL_WAVE" in result["labels"]
    assert len(result["isolation"]) == 1  # 2 ports -> 1 pair
    assert len(result["feeds_mm"]) == 2
    assert result["lw_grid_m"] is not None


def test_synthesize_ports_rejects_mismatched_beam_count(client, grounded_library):
    body = {
        "formulation": "B", "frequency_hz": _FREQ_HZ, "period_mm": 8.0,
        "layout": {"preset": "1xn", "n": 2, "cell_side_mm": 40.0},
        "ports": [{"beams": [{"theta_deg": 20.0, "phi_deg": 0.0}]}],  # only 1, layout has 2
        "library_id": grounded_library,
        "port_eps_r": _SUBSTRATE["eps_r"], "port_thickness_mm": _SUBSTRATE["thickness"] * 1e3,
    }
    resp = client.post("/api/hmsa/synthesize_ports", json=body)
    job = _wait_for_job(client, resp.json()["id"], timeout=60.0)
    assert job["status"] == "failed", job


def test_export_chains_off_a_ports_synthesis_job(client, grounded_library):
    synth_job = _synthesize_ports_1xn(client, grounded_library)
    assert synth_job["status"] == "done", synth_job
    export_body = {
        "source_job_id": synth_job["id"], "period_mm": 8.0,
        "eps_r": _SUBSTRATE["eps_r"], "tan_delta": _SUBSTRATE["tan_delta"],
        "thickness_mm": _SUBSTRATE["thickness"] * 1e3, "include_step": False,
    }
    resp = client.post("/api/export/array", json=export_body)
    job = _wait_for_job(client, resp.json()["id"], timeout=60.0)
    assert job["status"] == "done", job
    assert job["result"]["layout_csv_download_url"]
    readme = client.get(job["result"]["readme_download_url"])
    assert "Port 0: monopole at" in readme.text
    assert "S-parameter matrix" in readme.text


# ================================================================================== Coded (Part E)


def _upload_example_field(client):
    from pathlib import Path

    example = Path(__file__).resolve().parent.parent / "inputs" / "reference_fields" / "efl_phase_EzZ_example.fld"
    with open(example, "rb") as fh:
        resp = client.post("/api/hmsa/upload_field", files={"phase_file": ("efl_phase_EzZ_example.fld", fh, "text/plain")})
    return resp.json()["upload_id"]


def _synthesize_coded_one_bit_square_only(client, grounded_library):
    upload_id = _upload_example_field(client)
    body = {
        "upload_id": upload_id,
        "frequency_hz": _FREQ_HZ,
        "period_mm": 8.0,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        "n_bits": 1,
        "eps_r": _SUBSTRATE["eps_r"], "tan_delta": _SUBSTRATE["tan_delta"],
        "thickness_mm": _SUBSTRATE["thickness"] * 1e3,
        "library_id": grounded_library,
        "two_state": {"source": "square_only"},
    }
    resp = client.post("/api/hmsa/synthesize_coded", json=body)
    return _wait_for_job(client, resp.json()["id"], timeout=60.0)


def test_synthesize_coded_one_bit_square_only(client, grounded_library):
    job = _synthesize_coded_one_bit_square_only(client, grounded_library)
    assert job["status"] == "done", job
    result = job["result"]
    assert result["n_bits"] == 1
    states = result["two_state"]["states"]
    assert [st["geometry"]["kind"] for st in states] == ["rectangle", "rectangle"]
    assert all(st["beta_rad_per_m"] > 0.0 for st in states)
    assert result["two_state"]["beta_bar_rad_per_m"] > 0.0
    assert "CODED_DESIGN_METHOD" in result["labels"]
    assert any(label.startswith("TWO_STATE_PHASE_SEPARATION(") for label in result["labels"])
    assert any(label.startswith("TWO_STATE_SURFACE_WAVE(") for label in result["labels"])
    assert result["quantization"]["textbook_db"] == pytest.approx(3.9224, abs=1e-3)
    assert result["pattern_checks"]["main_beam_theta_deg"] is not None
    assert result["lw_grid_m"] is not None


def test_export_chains_off_a_one_bit_coded_job(client, grounded_library):
    synth_job = _synthesize_coded_one_bit_square_only(client, grounded_library)
    assert synth_job["status"] == "done", synth_job
    export_body = {
        "source_job_id": synth_job["id"], "period_mm": 8.0,
        "eps_r": _SUBSTRATE["eps_r"], "tan_delta": _SUBSTRATE["tan_delta"],
        "thickness_mm": _SUBSTRATE["thickness"] * 1e3, "include_step": False,
    }
    resp = client.post("/api/export/array", json=export_body)
    job = _wait_for_job(client, resp.json()["id"], timeout=60.0)
    assert job["status"] == "done", job
    assert job["result"]["layout_csv_download_url"]


def test_synthesize_coded_uploaded_two_state(client):
    upload_id = _upload_example_field(client)
    body = {
        "upload_id": upload_id,
        "frequency_hz": _FREQ_HZ,
        "period_mm": 8.0,
        "object_beams": [{"theta_deg": 0.0, "phi_deg": 0.0}],
        "n_bits": 1,
        "eps_r": _SUBSTRATE["eps_r"], "tan_delta": _SUBSTRATE["tan_delta"],
        "thickness_mm": _SUBSTRATE["thickness"] * 1e3,
        "two_state": {
            "source": "uploaded", "uploaded_label_0": "full patch", "uploaded_value_0_real": 1.0,
            "uploaded_value_0_imag": 0.0, "uploaded_label_1": "slotted patch",
            "uploaded_value_1_real": -0.98, "uploaded_value_1_imag": 0.0,
        },
    }
    resp = client.post("/api/hmsa/synthesize_coded", json=body)
    job = _wait_for_job(client, resp.json()["id"], timeout=60.0)
    # a state at the patch-plane resonance (Gamma = +1, 0 deg) is refused, so no valid pair remains
    assert job["status"] == "failed", job
    assert "NEAR_PATCH_PLANE_RESONANCE" in str(job)


def test_synthesize_coded_n_bit_matches_against_continuous_library(client, grounded_library):
    upload_id = _upload_example_field(client)
    body = {
        "upload_id": upload_id,
        "frequency_hz": _FREQ_HZ,
        "period_mm": 8.0,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        "n_bits": 2,
        "library_id": grounded_library,
    }
    resp = client.post("/api/hmsa/synthesize_coded", json=body)
    job = _wait_for_job(client, resp.json()["id"], timeout=60.0)
    assert job["status"] == "done", job
    result = job["result"]
    assert result["n_bits"] == 2
    assert result["assignments"] is not None
    assert result["lw_grid_m"] is None


@pytest.mark.slow
def test_ensure_annular_slot_library_builds_a_real_library(client):
    body = {
        "period_mm": 6.0, "eps_r": 2.2, "tan_delta": 0.0009, "thickness_mm": 1.575,
        "frequency_hz": 12.5e9, "a_mm": 4.0, "b_mm": 4.0, "r_mm": 1.0, "accuracy": "fast",
    }
    resp = client.post("/api/hmsa/ensure_annular_slot_library", json=body)
    job = _wait_for_job(client, resp.json()["id"], timeout=300.0)
    assert job["status"] == "done", job
    assert job["result"]["n_geometries"] == 8


@pytest.mark.slow
def test_synthesize_coded_square_and_annular_two_state(client, grounded_library):
    annular_body = {
        "period_mm": 8.0, "eps_r": _SUBSTRATE["eps_r"], "tan_delta": _SUBSTRATE["tan_delta"],
        "thickness_mm": _SUBSTRATE["thickness"] * 1e3, "frequency_hz": _FREQ_HZ,
        "a_mm": 4.0, "b_mm": 4.0, "r_mm": 1.0, "accuracy": "fast",
    }
    resp = client.post("/api/hmsa/ensure_annular_slot_library", json=annular_body)
    annular_job = _wait_for_job(client, resp.json()["id"], timeout=300.0)
    assert annular_job["status"] == "done", annular_job
    annular_library_id = annular_job["result"]["annular_slot_library_id"]

    upload_id = _upload_example_field(client)
    body = {
        "upload_id": upload_id,
        "frequency_hz": _FREQ_HZ,
        "period_mm": 8.0,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}],
        "n_bits": 1,
        "eps_r": _SUBSTRATE["eps_r"], "tan_delta": _SUBSTRATE["tan_delta"],
        "thickness_mm": _SUBSTRATE["thickness"] * 1e3,
        "library_id": grounded_library,
        "two_state": {"source": "square_and_annular", "annular_slot_library_id": annular_library_id},
    }
    resp = client.post("/api/hmsa/synthesize_coded", json=body)
    job = _wait_for_job(client, resp.json()["id"], timeout=60.0)
    assert job["status"] == "done", job
    result = job["result"]

    # also search the square library ALONE, as a baseline
    square_only_body = {**body, "two_state": {"source": "square_only"}}
    resp2 = client.post("/api/hmsa/synthesize_coded", json=square_only_body)
    job2 = _wait_for_job(client, resp2.json()["id"], timeout=60.0)
    assert job2["status"] == "done", job2

    # adding the annular-slot candidates to the search pool can never make the best achievable
    # separation WORSE (a real, measured property of a brute-force best-of-N search over a
    # superset) -- whether the winning pair actually includes an annular cell this time is a real
    # library-dependent outcome, not asserted either way.
    combined_sep = result["two_state"]["achieved_separation_deg"]
    square_only_sep = job2["result"]["two_state"]["achieved_separation_deg"]
    assert combined_sep >= square_only_sep - 1e-9


# ============================================================= per-port pattern prediction (Part C)


def test_pattern_predict_with_port_index(client, grounded_library):
    synth_job = _synthesize_ports_1xn(client, grounded_library)
    assert synth_job["status"] == "done", synth_job
    pattern_body = {
        "source_job_id": synth_job["id"], "frequency_hz": _FREQ_HZ, "pattern_model": "vector",
        "eps_r": _SUBSTRATE["eps_r"], "tan_delta": _SUBSTRATE["tan_delta"],
        "thickness_mm": _SUBSTRATE["thickness"] * 1e3,
        "beams": [{"theta_deg": 20.0, "phi_deg": 0.0}], "port_index": 0,
        "n_theta": 31, "n_phi": 37,
    }
    resp = client.post("/api/pattern/predict", json=pattern_body)
    job = _wait_for_job(client, resp.json()["id"], timeout=60.0)
    assert job["status"] == "done", job
    assert job["result"]["result"]["directivity_estimate_dbi"] is not None
    assert len(job["result"]["per_beam"]) == 1


def test_pattern_predict_port_beams_report_achieved_direction_and_xpd(client, grounded_library):
    synth_job = _synthesize_ports_1xn(client, grounded_library)
    assert synth_job["status"] == "done", synth_job
    pattern_body = {
        "source_job_id": synth_job["id"], "frequency_hz": _FREQ_HZ, "pattern_model": "vector",
        "eps_r": _SUBSTRATE["eps_r"], "tan_delta": _SUBSTRATE["tan_delta"],
        "thickness_mm": _SUBSTRATE["thickness"] * 1e3,
        "beams": [{"theta_deg": 20.0, "phi_deg": 0.0, "desired_polarization": "rhcp"}], "port_index": 0,
        "n_theta": 31, "n_phi": 37,
    }
    resp = client.post("/api/pattern/predict", json=pattern_body)
    job = _wait_for_job(client, resp.json()["id"], timeout=60.0)
    assert job["status"] == "done", job
    entry = job["result"]["per_beam"][0]
    assert entry["achieved_theta_deg"] == pytest.approx(20.0, abs=10.0)
    assert entry["level_db"] <= 0.0
    assert entry["integrated_xpd_db"] is not None


def test_pattern_predict_rejects_out_of_range_port_index(client, grounded_library):
    synth_job = _synthesize_ports_1xn(client, grounded_library)
    assert synth_job["status"] == "done", synth_job
    pattern_body = {
        "source_job_id": synth_job["id"], "frequency_hz": _FREQ_HZ, "pattern_model": "vector",
        "eps_r": _SUBSTRATE["eps_r"], "tan_delta": _SUBSTRATE["tan_delta"],
        "thickness_mm": _SUBSTRATE["thickness"] * 1e3, "port_index": 99,
    }
    resp = client.post("/api/pattern/predict", json=pattern_body)
    assert resp.status_code == 422


_PAPER_PAIR_26GHZ_PATCH_PLANE = {  # (phase deg, magnitude dB): empty cell, slotted cell (examples/make_paper_board_26ghz.py)
    "empty": (48.955768, -0.010430),
    "slotted": (-135.156476, -0.008195),
}


def _paper_pair_uploaded():
    def gamma(phase_deg, mag_db):
        return 10.0 ** (mag_db / 20.0) * np.exp(1j * np.radians(phase_deg))

    empty = gamma(*_PAPER_PAIR_26GHZ_PATCH_PLANE["empty"])
    slotted = gamma(*_PAPER_PAIR_26GHZ_PATCH_PLANE["slotted"])
    return {
        "source": "uploaded", "uploaded_label_0": "empty",
        "uploaded_value_0_real": empty.real, "uploaded_value_0_imag": empty.imag,
        "uploaded_label_1": "2.17 mm square, annular slot",
        "uploaded_value_1_real": slotted.real, "uploaded_value_1_imag": slotted.imag,
        "uploaded_geometry_0": {"kind": "empty"},
        "uploaded_geometry_1": {"kind": "annular_slot", "a_mm": 2.17, "b_mm": 2.17, "R_mm": 0.75, "r_mm": 0.46},
    }


def _synthesize_ports_coded_paper_pair(client, beta_mode):
    body = {
        "formulation": "coded", "frequency_hz": 26e9, "period_mm": 2.57,
        "layout": {"preset": "1xn", "n": 1, "cell_side_mm": 10.0, "gap_mm": 0.0},
        "ports": [{"beams": [{"theta_deg": 20.0, "phi_deg": 0.0, "desired_polarization": "lhcp"}]}],
        "port_eps_r": 2.2, "port_thickness_mm": 1.575, "port_tan_delta": 0.0009,
        "two_state": _paper_pair_uploaded(), "beta_mode": beta_mode,
    }
    resp = client.post("/api/hmsa/synthesize_ports", json=body)
    return _wait_for_job(client, resp.json()["id"], timeout=60.0)


def test_synthesize_ports_coded_paper_beta_mode_accepts_the_paper_pair(client):
    job = _synthesize_ports_coded_paper_pair(client, "paper")
    assert job["status"] == "done", job
    two_state = job["result"]["two_state"]
    assert two_state["beta_mode"] == "paper"
    assert two_state["beta_used_rad_per_m"] == pytest.approx(689.27, abs=0.5)
    assert 640.0 < two_state["beta_model_rad_per_m"] < 700.0
    assert two_state["state_betas_rad_per_m"][0] < two_state["state_betas_rad_per_m"][1]


def test_synthesize_ports_coded_state_beta_mode_refuses_the_paper_pair_on_beta_spread(client):
    job = _synthesize_ports_coded_paper_pair(client, "state")
    assert job["status"] == "failed"
    assert "NO_VALID_PAIR" in str(job)


def test_synthesize_ports_pointing_correction_reports_before_and_after(client):
    body = {
        "formulation": "coded", "frequency_hz": 26e9, "period_mm": 2.57,
        "layout": {"preset": "1xn", "n": 1, "cell_side_mm": 10.0, "gap_mm": 0.0},
        "ports": [{"beams": [{"theta_deg": 20.0, "phi_deg": 0.0, "desired_polarization": "lhcp"}]}],
        "port_eps_r": 2.2, "port_thickness_mm": 1.575, "port_tan_delta": 0.0009,
        "two_state": _paper_pair_uploaded(), "beta_mode": "paper",
        "pointing_correction": True, "pointing_iterations": 2,
    }
    resp = client.post("/api/hmsa/synthesize_ports", json=body)
    job = _wait_for_job(client, resp.json()["id"], timeout=60.0)
    assert job["status"] == "done", job
    (entry,) = job["result"]["refinement"]
    report = entry["pointing"]
    assert 0 <= report["passes"] <= 2
    assert report["after"][0]["miss_deg"] <= report["before"][0]["miss_deg"]


def test_synthesize_ports_matches_requested_beam_levels_by_default(client):
    body = {
        "formulation": "coded", "frequency_hz": 26e9, "period_mm": 2.57,
        "layout": {"preset": "1xn", "n": 1, "cell_side_mm": 10.0, "gap_mm": 0.0},
        "ports": [{"beams": [
            {"theta_deg": 20.0, "phi_deg": 0.0, "desired_polarization": "lhcp", "weight_real": 1.0},
            {"theta_deg": 20.0, "phi_deg": 180.0, "desired_polarization": "lhcp", "weight_real": 1.0},
        ]}],
        "port_eps_r": 2.2, "port_thickness_mm": 1.575, "port_tan_delta": 0.0009,
        "two_state": _paper_pair_uploaded(), "beta_mode": "paper",
    }
    resp = client.post("/api/hmsa/synthesize_ports", json=body)
    job = _wait_for_job(client, resp.json()["id"], timeout=60.0)
    assert job["status"] == "done", job
    balance = job["result"]["refinement"][0]["balance"]
    assert len(balance["before_db"]) == 2
    assert balance["max_error_after_db"] <= balance["max_error_before_db"]
