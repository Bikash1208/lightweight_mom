"""Phase 27 item 3: the tensor-CP design service and its app endpoint."""

import time

import pytest

from src.hologram import cell_family_registry as registry
from src.hologram import tensor_cp_service as svc


def test_single_default_reproduces_the_hfss_anchored_numbers():
    r = svc.run_design({"board": "single"})
    assert r["state"] == "OK" and r["n_cells"] == 1430
    p = r["ports"][0]
    assert p["gain_rhcp_or_lhcp_dbi"] == pytest.approx(23.03, abs=0.02)  # the one HFSS run: 23.03 dBi
    assert p["xpd_db"] == pytest.approx(16.9, abs=0.1) and p["sll_db"] == pytest.approx(-22.8, abs=0.1)
    assert p["sources"]["gain"] == "HFSS_ANCHORED" and r["validation_label"] == "SOLVER_SUPERCELL"
    assert any("SLL optimistic" in b for b in r["known_biases"])
    assert r["power_budget"]["aperture_radiated"] == pytest.approx(0.62 * r["eta_conv"])


def test_lhcp_is_the_mirror_of_rhcp_for_the_single_hexagon():
    rh = svc.run_design({"board": "single", "beams_deg": [[30.0, 225.0]]})["ports"][0]
    lh = svc.run_design({"board": "single", "handedness": "lhcp", "beams_deg": [[30.0, 135.0]]})["ports"][
        0
    ]  # y -> -y mirror of phi 225
    assert lh["gain_rhcp_or_lhcp_dbi"] == pytest.approx(rh["gain_rhcp_or_lhcp_dbi"], abs=0.05)
    assert lh["xpd_db"] == pytest.approx(rh["xpd_db"], abs=0.2) and lh["ar_peak_db"] == pytest.approx(
        rh["ar_peak_db"], abs=0.1
    )


def test_outside_domain_returns_named_state_and_never_runs():
    r = svc.run_design({"board": "single", "substrate": {"h_mm": 3.175}})
    assert r["state"] == registry.NOT_AVAILABLE and "substrate" in r["reasons"][0] and "ports" not in r
    r = svc.run_design({"board": "single", "frequency_ghz": 10.2})
    assert r["state"] == registry.NOT_AVAILABLE and "frequency" in r["reasons"][0]


def test_warnings_for_beams_pin_height_and_hand():
    r = svc.run_design(
        {"board": "single", "beams_deg": [[40.0, 225.0]], "pin_height_mm": 5.0, "handedness": "lhcp", "quick": True}
    )
    w = " ".join(r["warnings"])
    assert "BEAM_OUTSIDE_HFSS_ENVELOPE" in w and "PIN_HEIGHT_5_MM" in w and "LHCP_NOT_HFSS_CHECKED" in w
    assert r["launch_source"] == "MODEL_ONLY" and r["launch_fraction"] != 0.62


def test_request_validation():
    with pytest.raises(ValueError):
        svc.normalise_request({"board": "hex4_static", "beams_deg": [[30, 45]]})
    with pytest.raises(ValueError):
        svc.normalise_request({"board": "single", "handedness": "lp"})


def test_hex4_combined_quick_has_one_beam_row_and_four_ports_alone(tmp_path):
    r = svc.run_design({"board": "hex4_combined", "quick": True}, out_dir=tmp_path)
    assert r["state"] == "OK" and r["n_cells"] > 5600 and len(r["ports"]) == 1
    assert len(r["ports"][0]["each_port_alone_dbi"]) == 4
    assert any("QUICK_GRID" in w for w in r["warnings"])
    assert (tmp_path / "cell_map.png").exists()


def test_app_endpoint_runs_a_job_and_serves_the_registry():
    pytest.importorskip("fastapi")
    from fastapi.testclient import TestClient

    from src.app.main import create_app

    client = TestClient(create_app())
    fams = client.get("/api/tensor_cp/families").json()
    assert "tensor_slotted_circle_rt5880_6p4" in fams
    assert client.get("/api/tensor_cp/defaults?board=hex4_static").json()["beams_deg"][3] == [30.0, 315.0]
    job = client.post("/api/tensor_cp/design", json={"board": "single", "quick": True}).json()
    for _ in range(300):
        st = client.get(f"/api/jobs/{job['id']}").json()
        if st["status"] in ("done", "failed"):
            break
        time.sleep(0.2)
    assert st["status"] == "done", st
    res = st["result"]
    assert res["state"] == "OK" and res["files"]["cell_map_png"] == "cell_map.png"
    assert client.get(f"/api/files/{job['id']}/cell_map.png").status_code == 200
    assert (
        client.post("/api/tensor_cp/design", json={"board": "hex4_static", "beams_deg": [[30, 45]]}).status_code == 400
    )


@pytest.mark.slow
def test_hex4_combined_full_grid_matches_the_phase26_board_b_numbers():
    r = svc.run_design({"board": "hex4_combined"})
    p = r["ports"][0]
    assert (
        p["gain_rhcp_or_lhcp_dbi"] == pytest.approx(28.6, abs=0.1)
        and p["gain_range_over_edge_power_scale_dbi"][0] < p["gain_range_over_edge_power_scale_dbi"][1] + 1e-9
    )
    assert p["ar_peak_db"] == pytest.approx(2.5, abs=0.1) and p["xpd_db"] == pytest.approx(16.9, abs=0.1)


@pytest.mark.slow
def test_hex4_static_full_grid_matches_the_phase26_board_a_numbers():
    r = svc.run_design({"board": "hex4_static"})
    g = [p["gain_rhcp_or_lhcp_dbi"] for p in r["ports"]]
    assert g == pytest.approx([23.11, 22.97, 23.11, 22.97], abs=0.05)
