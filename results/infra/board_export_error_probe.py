"""Scratch (outside the repo): replay the board_export calls of the failing test and print the job's own error."""
import json

from tests.test_app_phase21_multiport import F, _grid_request  # noqa: F401
from tests.test_app_phase15_ports_and_coded import _wait_for_job, client, grounded_library  # noqa: F401


def test_show_job_error(client):
    body = _grid_request(32.0, 0.0)
    submitted = client.post("/api/hmsa/synthesize_ports", json=body).json()
    synth = _wait_for_job(client, submitted["id"], timeout=120.0)
    print("SYNTH STATUS", synth["status"])
    payload = {"source_job_id": submitted["id"], "period_mm": 6.4, "eps_r": 2.2, "tan_delta": 0.0009,
               "thickness_mm": 1.575, "frequency_hz": F, "name": "pins_test"}
    job = _wait_for_job(client, client.post("/api/export/board", json=payload).json()["id"], timeout=120.0)
    print("JOB_JSON_BEGIN")
    print(json.dumps(job, indent=1, default=str)[:6000])
    print("JOB_JSON_END")
