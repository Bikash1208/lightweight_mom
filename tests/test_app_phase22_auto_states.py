"""Phase 22: the automatic coded-cell states through the API (single aperture and multi-port), with an
analytic stand-in for the MoM square solve so the tests stay fast."""

import numpy as np
import pytest

fastapi = pytest.importorskip("fastapi")

from src.hologram import auto_states as auto_states_module  # noqa: E402
from tests.test_app_phase15_ports_and_coded import _upload_example_field, _wait_for_job, client  # noqa: E402,F401
from tests.test_app_phase21_multiport import F, _beam  # noqa: E402

_SUB = {"eps_r": 2.2, "tan_delta": 0.0009, "thickness_mm": 1.575}


@pytest.fixture(autouse=True)
def fake_mom(monkeypatch):
    """phase(L) = -42 - 260 (L - 5.85), limited to [-110, 150]; first-order mesh error 50 deg/mm * dx."""
    calls = []

    def solver(side_cells, n, substrate, frequency_hz, period_m):
        calls.append((side_cells, n))
        dx_mm = period_m / n * 1e3
        phase = float(np.clip(-42.0 - 260.0 * (side_cells * dx_mm - 5.85), -110.0, 150.0)) + 50.0 * dx_mm
        return 0.996 * np.exp(1j * np.radians(phase))

    monkeypatch.setattr(auto_states_module, "mom_square_gamma", solver)
    return calls


def _auto_states(client, **overrides):
    body = {**_SUB, "period_mm": 6.4, "frequency_hz": F, "n_bits": 1}
    body.update(overrides)
    resp = client.post("/api/hmsa/auto_states", json=body)
    assert resp.status_code == 200, resp.text
    return _wait_for_job(client, resp.json()["id"], timeout=60.0)


def test_auto_states_endpoint_returns_the_table_and_caches(client, fake_mom):
    job = _auto_states(client)
    assert job["status"] == "done", job
    result = job["result"]
    s0, s1 = result["states"]
    assert s0["geometry"] == {"kind": "empty"} and s0["phase_deg_patch_plane"] == pytest.approx(138.14, abs=0.05)
    assert s1["geometry"]["L_mm"] == pytest.approx(5.85, abs=0.05)
    assert s1["separation_deg"] == pytest.approx(180.0, abs=5.0)
    assert "NEAR_RESONANCE" in s1["labels"] and s1["sensitivity_deg_per_0p01mm"] > 1.0
    assert result["n_solves"] == 14 and result["n_cached"] == 0
    n_calls = len(fake_mom)
    again = _auto_states(client)["result"]
    assert again["n_solves"] == 0 and again["n_cached"] == 14 and len(fake_mom) == n_calls


def test_auto_states_endpoint_validates_bits_and_fixed_side(client):
    body = {**_SUB, "period_mm": 6.4, "frequency_hz": F}
    assert client.post("/api/hmsa/auto_states", json={**body, "n_bits": 4}).status_code == 422
    assert client.post("/api/hmsa/auto_states", json={**body, "n_bits": 2, "fixed_L_mm": 5.8}).status_code == 422
    fixed = _auto_states(client, fixed_L_mm=5.8)["result"]
    assert fixed["states"][1]["geometry"]["L_mm"] == 5.8


def _ports_body(**overrides):
    body = {
        "formulation": "coded",
        "frequency_hz": F,
        "period_mm": 6.4,
        "layout": {"preset": "grid", "rows": 2, "cols": 2, "cell_side_mm": 32.0, "gap_mm": 0.0,
                   "snap_side_to_cells": True},
        "ports": [{"beams": [_beam(30.0, phi)]} for phi in (45.0, 135.0, 225.0, 315.0)],
        "port_eps_r": 2.2,
        "port_thickness_mm": 1.575,
        "port_tan_delta": 0.0009,
        "two_state": {"source": "auto"},
        "beta_mode": "model",
    }
    body.update(overrides)
    return body


def test_multi_port_one_bit_with_the_automatic_pair(client):
    job = _wait_for_job(client, client.post("/api/hmsa/synthesize_ports", json=_ports_body()).json()["id"], timeout=120)
    assert job["status"] == "done", job
    result = job["result"]
    auto = result["two_state"]["auto"]
    assert auto["states"][1]["geometry"]["L_mm"] == pytest.approx(5.85, abs=0.05)
    assert result["two_state"]["achieved_separation_deg"] == pytest.approx(180.0, abs=5.0)
    assert all(p["patch_count"] > 0 and p["patches_mm"]["cells"][0]["geometry"][0] == pytest.approx(5.85, abs=0.05)
               for p in result["ports"])


def test_multi_port_two_bit_with_automatic_states_needs_no_library(client):
    body = _ports_body(n_bits=2, layout={"preset": "grid", "rows": 1, "cols": 2, "cell_side_mm": 32.0, "gap_mm": 0.0},
                       ports=[{"beams": [_beam(25.0, 0.0)]}, {"beams": [_beam(25.0, 90.0)]}])
    job = _wait_for_job(client, client.post("/api/hmsa/synthesize_ports", json=body).json()["id"], timeout=120)
    assert job["status"] == "done", job
    result = job["result"]
    assert [s["level"] for s in result["auto_states"]["states"]] == [0, 1, 2, 3]
    assert len(result["state_table"]) == 4 and all(p["patch_count"] > 0 for p in result["ports"])


def test_single_aperture_coded_with_the_automatic_pair(client):
    upload_id = _upload_example_field(client)
    body = {
        "upload_id": upload_id, "frequency_hz": F, "period_mm": 6.4,
        "object_beams": [{"theta_deg": 20.0, "phi_deg": 0.0}], "n_bits": 1,
        "eps_r": 2.2, "tan_delta": 0.0009, "thickness_mm": 1.575, "two_state": {"source": "auto"},
    }
    job = _wait_for_job(client, client.post("/api/hmsa/synthesize_coded", json=body).json()["id"], timeout=120)
    assert job["status"] == "done", job
    assert job["result"]["two_state"]["auto"]["states"][1]["geometry"]["L_mm"] == pytest.approx(5.85, abs=0.05)


def test_auto_states_endpoint_min_gap_setting(client):
    narrow = _auto_states(client, n_bits=2, min_gap_mm=0.4)["result"]
    assert narrow["min_gap_mm"] == 0.4
    assert max(s["geometry"].get("L_mm", 0) for s in narrow["states"]) <= 6.0 + 1e-9
    assert any("minimum gap 0.4 mm" in s["message"] for s in narrow["states"])
    body = {**_SUB, "period_mm": 6.4, "frequency_hz": F, "min_gap_mm": 3.5}
    assert client.post("/api/hmsa/auto_states", json=body).status_code == 422


def test_fixed_side_reports_its_separation_and_warns_when_the_pair_is_weak(client):
    weak = _auto_states(client, fixed_L_mm=4.0)["result"]
    assert weak["fixed_L_mm"] == 4.0
    assert weak["states"][1]["separation_deg"] < 150.0
    assert len(weak["warnings"]) == 1
    assert "Weak pair: 4.00 mm gives" in weak["warnings"][0] and "leave the fixed side blank to search" in weak["warnings"][0]
    good = _auto_states(client, fixed_L_mm=5.85)["result"]
    assert good["states"][1]["separation_deg"] == pytest.approx(180.0, abs=5.0) and good["warnings"] == []
    searched = _auto_states(client)["result"]
    assert searched["fixed_L_mm"] is None and searched["warnings"] == []


def test_state_labels_format_the_side_with_two_decimals(client):
    result = _auto_states(client, n_bits=2)["result"]
    import re

    for st in result["states"][1:]:
        assert re.fullmatch(r"square patch \d+\.\d{2} mm", st["label"]), st["label"]
