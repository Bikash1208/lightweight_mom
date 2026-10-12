"""Phase 23 bit-identical gate: hashes the complete `/synthesize_ports` job result (every number, every
array, byte for byte) and the exact-board export of several RECTANGLE layouts. Run on the code before and
after the polygon-port change; the two JSON files must be identical.

    python scripts/research/phase23_rect_snapshot.py <out.json>
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.app.routers.hmsa import (  # noqa: E402
    HmsaSynthesizePortsRequest,
    TwoStateSourceRequest,
    _run_synthesize_ports_job,
)
from src.export.exact_board import exact_board_from_ports_result  # noqa: E402

F = 10.54e9
_G0 = np.exp(1j * np.radians(138.14)) * 10 ** (-0.0005 / 20)
_G1 = np.exp(1j * np.radians(-42.0)) * 10 ** (-0.033 / 20)


class _Job:
    def __init__(self, job_id):
        self.id = job_id
        self.progress = 0.0
        self.progress_message = ""


def _pair():
    return TwoStateSourceRequest(
        source="uploaded",
        uploaded_label_0="empty",
        uploaded_value_0_real=_G0.real,
        uploaded_value_0_imag=_G0.imag,
        uploaded_label_1="square 5.85",
        uploaded_value_1_real=_G1.real,
        uploaded_value_1_imag=_G1.imag,
        uploaded_geometry_0={"kind": "empty"},
        uploaded_geometry_1={"kind": "rectangle", "L_mm": 5.85, "W_mm": 5.85},
    ).model_dump()


def _beam(theta, phi, pol="rhcp"):
    return {
        "theta_deg": theta,
        "phi_deg": phi,
        "weight_real": 1.0,
        "weight_imag": 0.0,
        "initial_phase_deg": 0.0,
        "desired_polarization": pol,
        "oam_charge": 0,
    }


def _body(layout, ports, **kw):
    body = {
        "formulation": "coded",
        "frequency_hz": F,
        "period_mm": 6.4,
        "layout": layout,
        "ports": ports,
        "port_eps_r": 2.2,
        "port_thickness_mm": 1.575,
        "port_tan_delta": 0.0009,
        "two_state": _pair(),
        "beta_mode": "model",
    }
    body.update(kw)
    return HmsaSynthesizePortsRequest(**body).model_dump()


FOUR = [{"beams": [_beam(30.0, p)]} for p in (45.0, 135.0, 225.0, 315.0)]
CASES = {
    "grid2x2_gap0": _body(
        {"preset": "grid", "rows": 2, "cols": 2, "cell_side_mm": 64.0, "gap_mm": 0.0, "snap_side_to_cells": True}, FOUR
    ),
    "grid2x2_gap28": _body(
        {"preset": "grid", "rows": 2, "cols": 2, "cell_side_mm": 64.0, "gap_mm": 28.0, "snap_side_to_cells": True},
        FOUR,
        pointing_correction=True,
    ),
    "grid1x3_two_beams": _body(
        {"preset": "grid", "rows": 1, "cols": 3, "cell_side_mm": 51.2, "gap_mm": 6.4, "snap_side_to_cells": True},
        [{"beams": [_beam(20.0, 0.0), _beam(25.0, 180.0)]} for _ in range(3)],
    ),
    "onexn": _body({"preset": "1xn", "n": 3, "cell_side_mm": 44.8, "gap_mm": 0.0}, FOUR[:3], feed_keepout=False),
    "mishra": _body(
        {"preset": "2x2_plus_1_large", "small_side_mm": 38.4, "gap_mm": 3.2},
        [{"beams": [_beam(25.0, 72.0 * k)]} for k in range(5)],
    ),
    "custom": _body(
        {"preset": "custom", "regions_mm": [[0, 51.2, 0, 44.8], [51.2, 96.0, 0, 44.8, 70.0, 20.0]]}, FOUR[:2]
    ),
    "nbit2": _body(
        {"preset": "grid", "rows": 1, "cols": 2, "cell_side_mm": 51.2, "gap_mm": 0.0, "snap_side_to_cells": True},
        FOUR[:2],
        n_bits=2,
        two_state={"source": "auto"},
    )
    if False
    else None,
}


def _digest(obj, h):
    if isinstance(obj, np.ndarray):
        h.update(b"nd" + str(obj.dtype).encode() + str(obj.shape).encode() + np.ascontiguousarray(obj).tobytes())
    elif isinstance(obj, dict):
        for k in sorted(obj, key=str):
            if k in ("patch_layout_png_url",):
                continue
            h.update(b"k" + str(k).encode())
            _digest(obj[k], h)
    elif isinstance(obj, (list, tuple)):
        h.update(b"l" + str(len(obj)).encode())
        for v in obj:
            _digest(v, h)
    elif isinstance(obj, (float, np.floating)):
        h.update(b"f" + float(obj).hex().encode())
    elif isinstance(obj, complex):
        h.update(b"c" + obj.real.hex().encode() + obj.imag.hex().encode())
    else:
        h.update(b"o" + repr(obj).encode())


def main(out):
    out_d = {}
    for name, body in CASES.items():
        if body is None:
            continue
        res = _run_synthesize_ports_job(_Job(f"p23snap_{name}"), body)
        h = hashlib.sha256()
        _digest(res, h)
        eb = exact_board_from_ports_result(res, body["period_mm"], body["port_thickness_mm"])
        h2 = hashlib.sha256()
        for p in eb.patches:
            h2.update(f"{p.port},{p.i},{p.j},{p.x0},{p.y0},{p.x1},{p.y1};".encode())
        h2.update(repr((eb.feeds_mm, eb.port_rects_mm, eb.x_min, eb.x_max, eb.y_min, eb.y_max)).encode())
        out_d[name] = {
            "result_sha256": h.hexdigest(),
            "exact_board_sha256": h2.hexdigest(),
            "n_patches": len(eb.patches),
            "n_ports": len(res["ports"]),
            "n_cells": [p["n_cells"] for p in res["ports"]],
        }
        print(
            name, out_d[name]["result_sha256"][:16], out_d[name]["exact_board_sha256"][:16], len(eb.patches), flush=True
        )
    Path(out).write_text(json.dumps(out_d, indent=1, sort_keys=True))


if __name__ == "__main__":
    main(sys.argv[1])
