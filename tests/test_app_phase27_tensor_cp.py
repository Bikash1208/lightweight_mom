"""Phase 27 close-out: the tensor-CP flow end to end through the app on the single hexagon -- design job, export job,
HFSS comparison with the committed Phase 26 HFSS files (23.03 dBi, XPD 23.0 dB, S11 -15.0 dB at 10.54 GHz)."""

import json
import time
from pathlib import Path

import pytest

pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from src.app.main import create_app  # noqa: E402
from tests._platform_tol import assert_matches_reference, assert_text_matches_reference  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
HFSS = ROOT / "hfss_results/hex150_tensor_cp_single"


def _wait(client, job):
    for _ in range(900):
        st = client.get(f"/api/jobs/{job['id']}").json()
        if st["status"] in ("done", "failed", "cancelled"):
            return st
        time.sleep(0.2)
    raise AssertionError("job timed out")


def test_single_hexagon_design_export_compare_end_to_end():
    client = TestClient(create_app())
    design = _wait(client, client.post("/api/tensor_cp/design", json={"board": "single"}).json())
    assert design["status"] == "done" and design["result"]["state"] == "OK"
    p = design["result"]["ports"][0]
    assert p["gain_rhcp_or_lhcp_dbi"] == pytest.approx(23.03, abs=0.02) and p["sources"]["gain"] == "HFSS_ANCHORED"

    export = _wait(client, client.post("/api/tensor_cp/export", json={"job_id": design["id"]}).json())
    assert export["status"] == "done"
    files = export["result"]["files"]
    assert {"hfss_build_board.py", "layout.csv", "layout.dxf", "predictions.json", "README.md"} <= set(files)
    built = client.get(export["result"]["download"]["hfss_build_board.py"]).content
    assert_text_matches_reference(
        built, (ROOT / "hfss_exports/hex150_tensor_cp_single/hfss_build_board.py").read_bytes(), "hfss_build_board.py"
    )

    uploads = {
        "gain_csv": ("g.csv", (HFSS / "realized_gain_RHCP_LHCP_10p5GHz_1deg.csv").read_bytes()),
        "params_txt": ("p.txt", (HFSS / "antenna_params_10p5GHz.txt").read_bytes()),
        "s_param": ("S11.csv", (HFSS / "S11_8_12GHz.csv").read_bytes()),
    }
    r = client.post("/api/tensor_cp/compare", data={"board": "single", "design_job_id": design["id"]}, files=uploads)
    assert r.status_code == 200, r.text
    rep = r.json()["report"]
    assert rep["measured"]["peak_rhcp_gain_dbi"] == pytest.approx(23.03, abs=0.01)
    assert rep["measured"]["xpd_at_peak_db"] == pytest.approx(23.0, abs=0.1)
    assert rep["s_parameter"]["s_db_at_design"] == pytest.approx(-15.0, abs=0.05)
    assert rep["checks"]["peak_gain"]["state"] == "GAIN_WITHIN_PREDICTED_RANGE"
    used = json.loads(Path(r.json()["predictions_used"]).read_text())
    assert_matches_reference(
        used,
        json.loads((ROOT / "hfss_exports/hex150_tensor_cp_single/predictions.json").read_text()),
        "predictions.json",
    )


def test_unavailable_request_ends_in_the_named_state_and_cannot_be_exported():
    client = TestClient(create_app())
    st = _wait(
        client, client.post("/api/tensor_cp/design", json={"board": "single", "substrate": {"h_mm": 3.175}}).json()
    )
    assert st["status"] == "done" and st["result"]["state"] == "NOT_AVAILABLE_FROM_CURRENT_DATA"
    assert client.post("/api/tensor_cp/export", json={"job_id": st["id"]}).status_code == 409
