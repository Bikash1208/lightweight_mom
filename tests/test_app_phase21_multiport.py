"""Phase 21: flexible multi-port layouts and the coded (1-bit / N-bit) method for multi-port --
the API surface (grid layout, keep-out, per-port report, N-bit) and the acceptance reproduction of
the hand-built four-port 300.8 mm square-cell board (`hfss_exports/four_port_1bit_square_300mm_gap{0,28}`)."""

import numpy as np
import pytest

fastapi = pytest.importorskip("fastapi")

from src.app.routers.hmsa import HmsaSynthesizePortsRequest, TwoStateSourceRequest, _run_synthesize_ports_job  # noqa: E402
from tests.test_app_phase15_ports_and_coded import (  # noqa: E402,F401
    _FREQ_HZ,
    _wait_for_job,
    client,
    grounded_library,
)

F = 10.54e9
# patch-plane reflection of the reference board's pair: the empty grounded slab and a 5.85 mm square
_G_EMPTY = np.exp(1j * np.radians(138.14)) * 10 ** (-0.0005 / 20)
_G_SQUARE = np.exp(1j * np.radians(-42.0)) * 10 ** (-0.033 / 20)


def _square_pair():
    return TwoStateSourceRequest(
        source="uploaded",
        uploaded_label_0="empty (bare substrate)",
        uploaded_value_0_real=_G_EMPTY.real,
        uploaded_value_0_imag=_G_EMPTY.imag,
        uploaded_label_1="square patch 5.85 mm",
        uploaded_value_1_real=_G_SQUARE.real,
        uploaded_value_1_imag=_G_SQUARE.imag,
        uploaded_geometry_0={"kind": "empty"},
        uploaded_geometry_1={"kind": "rectangle", "L_mm": 5.85, "W_mm": 5.85},
    ).model_dump()


def _beam(theta, phi):
    return {
        "theta_deg": theta,
        "phi_deg": phi,
        "weight_real": 1.0,
        "weight_imag": 0.0,
        "initial_phase_deg": 0.0,
        "desired_polarization": "rhcp",
        "oam_charge": 0,
    }


class _Job:
    def __init__(self, job_id):
        self.id = job_id
        self.progress = 0.0
        self.progress_message = ""


def _grid_request(side_mm, gap_mm, **overrides):
    body = {
        "formulation": "coded",
        "frequency_hz": F,
        "period_mm": 6.4,
        "layout": {"preset": "grid", "rows": 2, "cols": 2, "cell_side_mm": side_mm, "gap_mm": gap_mm,
                   "snap_side_to_cells": True},
        "ports": [{"beams": [_beam(30.0, phi)]} for phi in (45.0, 135.0, 225.0, 315.0)],
        "port_eps_r": 2.2,
        "port_thickness_mm": 1.575,
        "port_tan_delta": 0.0009,
        "two_state": _square_pair(),
        "beta_mode": "model",
    }
    body.update(overrides)
    return HmsaSynthesizePortsRequest(**body).model_dump()


def test_coded_one_bit_grid_small_board(client, tmp_path):
    """Fast: 2x2 grid of 5x5-cell ports, 1-bit coded, uploaded square pair."""
    result = _run_synthesize_ports_job(_Job("p21small"), _grid_request(32.0, 0.0))
    assert result["n_bits"] == 1 and len(result["ports"]) == 4
    # port order matches the reference board: P1 (+,+), P2 (-,+), P3 (-,-), P4 (+,-)
    signs = [(p["feed_x_mm"] > 0, p["feed_y_mm"] > 0) for p in result["ports"]]
    assert signs == [(True, True), (False, True), (False, False), (True, False)]
    for p in result["ports"]:
        assert p["n_cells"] == 25  # 32 mm / 6.4 mm, all five rows and columns kept
        assert p["patch_count"] == len(p["patches_mm"]["cells"])
        assert 0 <= p["keepout_cells_changed"] <= 1  # feed at a cell centre: at most the 1 cell under it
        rep = p["report"]
        assert rep["beams"][0]["requested_theta_deg"] == 30.0
        assert rep["beams"][0]["pointing_error_deg"] < 10.0
        assert rep["beams"][0]["cross_pol_level_db"] is not None
        assert p["patches_mm"]["cells"][0]["geometry"] == [5.85, 5.85]
    assert result["two_state"]["beta_mode"] == "model"


def test_feed_keepout_empties_the_cell_under_the_feed_and_reports_the_count(client):
    on = _run_synthesize_ports_job(_Job("p21ko_on"), _grid_request(32.0, 0.0))
    off = _run_synthesize_ports_job(_Job("p21ko_off"), _grid_request(32.0, 0.0, feed_keepout=False))
    for p_on, p_off in zip(on["ports"], off["ports"], strict=True):
        assert p_off["keepout_cells_changed"] is None
        assert p_on["patch_count"] == p_off["patch_count"] - p_on["keepout_cells_changed"]
        centre = (p_on["patches_mm"]["nx"] // 2, p_on["patches_mm"]["ny"] // 2)
        assert not any((c["j"], c["i"]) == centre for c in p_on["patches_mm"]["cells"])


def test_feed_keepout_on_a_cell_corner_clears_up_to_four_cells(client):
    even = _run_synthesize_ports_job(_Job("p21ko_even"), _grid_request(38.4, 0.0, feed_keepout=True))  # 6x6 cells
    assert all(0 <= p["keepout_cells_changed"] <= 4 for p in even["ports"])
    assert any(p["keepout_cells_changed"] > 0 for p in even["ports"])


def test_coded_two_bit_multi_port_uses_the_library_and_reports_the_state_table(client, grounded_library):
    body = {
        "formulation": "coded",
        "n_bits": 2,
        "frequency_hz": _FREQ_HZ,
        "period_mm": 8.0,
        "library_id": grounded_library,
        "layout": {"preset": "grid", "rows": 1, "cols": 2, "cell_side_mm": 40.0, "gap_mm": 0.0},
        "ports": [{"beams": [_beam(20.0, 0.0)]}, {"beams": [_beam(25.0, 90.0)]}],
        "port_eps_r": 2.2,
        "port_thickness_mm": 1.575,
        "port_tan_delta": 0.0009,
    }
    resp = client.post("/api/hmsa/synthesize_ports", json=body)
    job = _wait_for_job(client, resp.json()["id"], timeout=120.0)
    assert job["status"] == "done", job
    result = job["result"]
    assert result["n_bits"] == 2 and result["two_state"] is None
    table = result["state_table"]
    assert [r["target_phase_deg"] for r in table] == [0.0, 90.0, 180.0, 270.0]
    assert sum(r["n_cells"] for r in table) == sum(p["n_cells"] for p in result["ports"])
    assert all(p["report"] is not None and p["patch_count"] > 0 for p in result["ports"])


def test_coded_n_bit_multi_port_needs_a_library(client):
    body = {
        "formulation": "coded",
        "n_bits": 2,
        "frequency_hz": _FREQ_HZ,
        "period_mm": 8.0,
        "layout": {"preset": "grid", "rows": 1, "cols": 1, "cell_side_mm": 40.0},
        "ports": [{"beams": [_beam(20.0, 0.0)]}],
        "port_eps_r": 2.2,
        "port_thickness_mm": 1.575,
        "port_tan_delta": 0.0009,
    }
    job = _wait_for_job(client, client.post("/api/hmsa/synthesize_ports", json=body).json()["id"])
    assert job["status"] == "failed" and "n_bits > 1 needs library_id" in str(job)


@pytest.mark.parametrize("gap_mm", [0.0, 28.0])
def test_acceptance_four_port_square_300mm_board(client, gap_mm):
    """The reference board of `hfss_exports/four_port_1bit_square_300mm_gap{0,28}`, built from the
    app's own grid layout: 47 x 47 cells per port, beams at theta = 30 deg within 1 deg, 27.1-27.8 dBi,
    cross-pol -16 to -20 dB, 4416 patches."""
    result = _run_synthesize_ports_job(_Job(f"p21acc{int(gap_mm)}"), _grid_request(300.0, gap_mm))
    assert all(p["n_cells"] == 47 * 47 for p in result["ports"])
    assert sum(p["patch_count"] for p in result["ports"]) == 4416
    for p, phi in zip(result["ports"], (45.0, 135.0, 225.0, 315.0), strict=True):
        beam = p["report"]["beams"][0]
        assert beam["theta_deg"] == pytest.approx(30.0, abs=1.0)
        assert beam["pointing_error_deg"] <= 1.0
        assert 27.0 <= beam["directivity_dbi"] <= 27.9
        assert -20.5 <= beam["cross_pol_level_db"] <= -16.0
    assert result["ports"][0]["feed_x_mm"] == pytest.approx(150.4 + gap_mm / 2.0)
