"""Phase 22 acceptance: the 4-port 300 mm grid run with the AUTOMATIC cell pair (real MoM) reproduces the
Phase 21 numbers. Slow: a fresh pair costs about 80 s of MoM."""

import pytest

fastapi = pytest.importorskip("fastapi")

from src.app.routers.hmsa import HmsaSynthesizePortsRequest, _run_synthesize_ports_job  # noqa: E402
from tests.test_app_phase21_multiport import F, _beam, _Job  # noqa: E402


@pytest.mark.slow
def test_four_port_300mm_grid_with_the_automatic_pair(tmp_path, monkeypatch):
    monkeypatch.setenv("LWMOM_APP_DATA_DIR", str(tmp_path / "appdata"))
    body = HmsaSynthesizePortsRequest(
        formulation="coded",
        frequency_hz=F,
        period_mm=6.4,
        layout={"preset": "grid", "rows": 2, "cols": 2, "cell_side_mm": 300.0, "gap_mm": 0.0, "snap_side_to_cells": True},
        ports=[{"beams": [_beam(30.0, phi)]} for phi in (45.0, 135.0, 225.0, 315.0)],
        port_eps_r=2.2,
        port_thickness_mm=1.575,
        port_tan_delta=0.0009,
        two_state={"source": "auto"},
        beta_mode="model",
    ).model_dump()
    result = _run_synthesize_ports_job(_Job("p22acc"), body)
    auto = result["two_state"]["auto"]
    assert auto["states"][1]["geometry"]["L_mm"] == pytest.approx(5.85, abs=0.1)
    assert auto["states"][1]["phase_deg_patch_plane"] == pytest.approx(-42.0, abs=5.0)
    assert result["two_state"]["achieved_separation_deg"] == pytest.approx(180.0, abs=5.0)
    for p in result["ports"]:
        assert p["n_cells"] == 47 * 47
        beam = p["report"]["beams"][0]
        assert beam["theta_deg"] == pytest.approx(30.0, abs=1.0)
        assert beam["pointing_error_deg"] <= 1.0
        assert 27.0 <= beam["directivity_dbi"] <= 28.0
