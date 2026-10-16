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


# ----------------------------------------------------------------------------------- item 4: exports
from pathlib import Path  # noqa: E402

from src.export.tensor_cp_export import export_app_result  # noqa: E402

_ROOT = Path(__file__).resolve().parents[1]


def _dxf_entities(path):
    """Geometry of a DXF as comparable tuples. ezdxf's text output is not byte-stable (timestamps, GUIDs, and the
    order of two internal objects change between runs), so entities are compared, not bytes."""
    import ezdxf

    out = []
    for e in ezdxf.readfile(path).modelspace():
        if e.dxftype() == "CIRCLE":
            out.append(
                ("CIRCLE", e.dxf.layer, round(e.dxf.center.x, 6), round(e.dxf.center.y, 6), round(e.dxf.radius, 6))
            )
        else:
            out.append((e.dxftype(), e.dxf.layer, tuple((round(x, 6), round(y, 6)) for x, y, *_ in e.get_points())))
    return out


def test_app_export_of_the_single_hexagon_matches_the_committed_phase26_folder(tmp_path):
    r = svc.run_design({"board": "single"})
    names = export_app_result(r, tmp_path)
    assert {
        "hfss_build_board.py",
        "layout.csv",
        "layout.json",
        "layout.dxf",
        "layout.png",
        "predictions.json",
        "README.md",
    } <= set(names)
    ref = _ROOT / "hfss_exports/hex150_tensor_cp_single"
    for f in ("hfss_build_board.py", "layout.csv", "layout.json", "predictions.json", "layout.png"):
        assert (tmp_path / f).read_bytes() == (ref / f).read_bytes(), f
    assert _dxf_entities(tmp_path / "layout.dxf") == _dxf_entities(ref / "layout.dxf")


def test_app_export_hex4_has_drive_port_option_and_four_ports(tmp_path):
    r = svc.run_design({"board": "hex4_combined", "quick": True})
    export_app_result(r, tmp_path, drive_port=3)
    src = (tmp_path / "hfss_build_board.py").read_text()
    assert "DRIVE_PORT = 3" in src and src.count("AssignLumpedPort") == 1 and "FEEDS = [(64.952, 112.500)" in src
    compile(src, "hfss_build_board.py", "exec")
    import json

    pred = json.loads((tmp_path / "predictions.json").read_text())
    assert pred["board"] == "hex4_combined" and pred["state"] == "PREDICTION_NOT_HFSS_VALIDATED"


def test_export_refuses_an_unavailable_design(tmp_path):
    with pytest.raises(ValueError):
        export_app_result(svc.run_design({"board": "single", "frequency_ghz": 9.0}), tmp_path)


def test_app_export_endpoint_produces_downloadable_files():
    pytest.importorskip("fastapi")
    from fastapi.testclient import TestClient

    from src.app.main import create_app

    client = TestClient(create_app())

    def wait(job):
        for _ in range(600):
            st = client.get(f"/api/jobs/{job['id']}").json()
            if st["status"] in ("done", "failed"):
                return st
            time.sleep(0.2)
        raise AssertionError("job timed out")

    d = wait(client.post("/api/tensor_cp/design", json={"board": "single", "quick": True}).json())
    assert d["status"] == "done"
    e = wait(client.post("/api/tensor_cp/export", json={"job_id": d["id"]}).json())
    assert e["status"] == "done" and "hfss_build_board.py" in e["result"]["files"]
    r = client.get(e["result"]["download"]["hfss_build_board.py"])
    assert r.status_code == 200 and b"InsertDesign" in r.content
    assert client.post("/api/tensor_cp/export", json={"job_id": "nope"}).status_code == 404
    assert client.post("/api/tensor_cp/export", json={"job_id": d["id"], "drive_port": 7}).status_code == 400
