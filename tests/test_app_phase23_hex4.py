"""Phase 23 A: a hex4 board end to end through `/synthesize_ports` (1-bit coded, uploaded square pair)
and the exact export -- cells only inside the hexagons, the polygon path reproduces the rectangle path
when the polygon IS the rectangle, trimmed STEP, DXF PORTS layer, README outlines."""

import numpy as np
import pytest

pytest.importorskip("fastapi")

from src.app.routers.hmsa import (  # noqa: E402
    HmsaSynthesizePortsRequest,
    TwoStateSourceRequest,
    _run_synthesize_ports_job,
)
from src.export.aedt_script import generate_aedt_script  # noqa: E402
from src.export.board_readme import build_board_readme  # noqa: E402
from src.export.exact_board import (  # noqa: E402
    exact_board_from_ports_result,
    verify_step_exact,
    write_exact_step,
    write_layout_dxf,
    write_layout_json,
)
from src.hologram.ports import PortRegion  # noqa: E402

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


def _beam(phi):
    return {
        "theta_deg": 30.0,
        "phi_deg": phi,
        "weight_real": 1.0,
        "weight_imag": 0.0,
        "initial_phase_deg": 0.0,
        "desired_polarization": "rhcp",
        "oam_charge": 0,
    }


def _request(layout, n_ports=4, **kw):
    body = {
        "formulation": "coded",
        "frequency_hz": F,
        "period_mm": 6.4,
        "layout": layout,
        "ports": [{"beams": [_beam(45.0 + 90.0 * k)]} for k in range(n_ports)],
        "port_eps_r": 2.2,
        "port_thickness_mm": 1.575,
        "port_tan_delta": 0.0009,
        "two_state": _pair(),
        "beta_mode": "model",
    }
    body.update(kw)
    return HmsaSynthesizePortsRequest(**body).model_dump()


@pytest.fixture(scope="module")
def hex_result():
    return _run_synthesize_ports_job(_Job("p23hexsmall"), _request({"preset": "hex4", "side_mm": 40.0, "gap_mm": 0.0}))


def test_hex4_job_keeps_only_cells_inside_each_hexagon(hex_result):
    res = hex_result
    assert len(res["ports"]) == 4 and res["layout_preset"] == "hex4" and res["trim_board"] is False
    from src.hologram.ports import regular_hexagon_vertices

    for p in res["ports"]:
        assert len(p["outline_mm"]) == 6
        verts = [(x * 1e-3, y * 1e-3) for x, y in p["outline_mm"]]
        region = PortRegion(
            0,
            min(v[0] for v in verts),
            max(v[0] for v in verts),
            min(v[1] for v in verts),
            max(v[1] for v in verts),
            p["feed_x_mm"] * 1e-3,
            p["feed_y_mm"] * 1e-3,
            vertices_m=tuple(verts),
        )
        x_m, y_m = p["x_m"], p["y_m"]
        inside = region.contains_array(x_m, y_m)
        field = np.asarray(p["aperture_field_real"]) + 1j * np.asarray(p["aperture_field_imag"])
        assert np.all(field[~inside] == 0.0) and np.all(np.abs(field[inside]) > 0.0)
        assert p["n_cells"] == int(inside.sum())
        for c in p["patches_mm"]["cells"]:  # every exported patch sits on a cell inside the hexagon
            assert inside[c["i"], c["j"]]
    assert regular_hexagon_vertices(0, 0, 1.0)  # (import used above)


def test_hex4_cells_share_one_lattice_and_never_overlap(hex_result):
    board = exact_board_from_ports_result(hex_result, 6.4, 1.575)
    # one shared lattice: every patch centre is at an odd multiple of half a period (origin on a cell corner)
    for p in board.patches[:400]:
        cx = float((p.x0 + p.x1) / 2)
        cy = float((p.y0 + p.y1) / 2)
        assert abs((cx / 6.4) - round(cx / 6.4 - 0.5) - 0.5) < 1e-6
        assert abs((cy / 6.4) - round(cy / 6.4 - 0.5) - 0.5) < 1e-6
    centres = {(float(p.x0 + p.x1), float(p.y0 + p.y1)) for p in board.patches}
    assert len(centres) == len(board.patches)  # no cell built twice (square patches on one lattice cannot overlap)


def test_polygon_equal_to_a_rectangle_reproduces_the_rectangle_design():
    """A polygon region that IS the rectangle takes the polygon code path (mask, `region=`) and must give
    exactly the rectangle's design."""
    from src.hologram import ports as ports_mod

    layout = {"preset": "grid", "rows": 1, "cols": 2, "cell_side_mm": 51.2, "gap_mm": 0.0, "snap_side_to_cells": True}
    body = _request(layout, n_ports=2)
    plain = _run_synthesize_ports_job(_Job("p23rect_plain"), body)

    original = ports_mod.layout_preset_grid

    def as_polygons(*a, **k):
        out = []
        for r in original(*a, **k):
            out.append(
                PortRegion(
                    r.port_id,
                    r.x_min_m,
                    r.x_max_m,
                    r.y_min_m,
                    r.y_max_m,
                    r.feed_x_m,
                    r.feed_y_m,
                    vertices_m=r.outline_m,
                )
            )
        return tuple(out)

    import src.app.routers.hmsa as hmsa

    hmsa.layout_preset_grid = as_polygons
    try:
        poly = _run_synthesize_ports_job(_Job("p23rect_poly"), body)
    finally:
        hmsa.layout_preset_grid = original
    for a, b in zip(plain["ports"], poly["ports"], strict=True):
        assert a["n_cells"] == b["n_cells"] and a["patch_count"] == b["patch_count"]
        assert np.array_equal(a["aperture_field_real"], b["aperture_field_real"])
        assert np.array_equal(a["aperture_field_imag"], b["aperture_field_imag"])
        assert a["patches_mm"] == b["patches_mm"]
        assert a["report"]["beams"] == b["report"]["beams"]
    assert plain["lw_grid_m"] == poly["lw_grid_m"]


def test_exact_export_of_a_hex_board(hex_result, tmp_path):
    board = exact_board_from_ports_result(hex_result, 6.4, 1.575)
    assert len(board.port_outlines_mm) == 4 and not board.trim_board
    step = write_exact_step(board, tmp_path / "board.step")
    assert verify_step_exact(step) > 0
    dxf = write_layout_dxf(board, tmp_path / "layout.dxf").read_text()
    assert "PORTS" in dxf and dxf.count("\nPORTS\n") == 24  # four hexagons, six edges each
    js = write_layout_json(board, tmp_path / "layout.json").read_text()
    assert '"port_outlines"' in js
    readme = build_board_readme(board, hex_result, "HEX", 2.2, 0.0009, F, ["board.step"])
    assert "Sub-aperture outlines" in readme and "bounding rectangle" in readme
    script = generate_aedt_script(board, "HEX", 2.2, 0.0009, F)
    assert "OUTLINES" not in script  # rectangular substrate when not trimmed


def test_trimmed_board_follows_the_outline(hex_result, tmp_path):
    trimmed = dict(hex_result, trim_board=True)
    board = exact_board_from_ports_result(trimmed, 6.4, 1.575)
    plain = exact_board_from_ports_result(hex_result, 6.4, 1.575)
    assert board.trim_board
    assert board.width_mm <= plain.width_mm and board.height_mm <= plain.height_mm
    step = write_exact_step(board, tmp_path / "trim.step")
    assert verify_step_exact(step) > 0
    script = generate_aedt_script(board, "HEX", 2.2, 0.0009, F)
    assert "OUTLINES = [[(" in script and script.count("poly_sheet(") >= 2
    readme = build_board_readme(board, trimmed, "HEX", 2.2, 0.0009, F, ["board.step"])
    assert "trimmed to these outlines" in readme
    dxf = write_layout_dxf(board, tmp_path / "trim.dxf").read_text()
    assert dxf.count("\nBOARD\n") == 24


def test_hex4_needs_exactly_one_size():
    from src.errors import LwmomError as ProblemError

    with pytest.raises(ProblemError):
        _run_synthesize_ports_job(_Job("p23bad"), _request({"preset": "hex4", "gap_mm": 0.0}))
    with pytest.raises(ProblemError):
        _run_synthesize_ports_job(
            _Job("p23bad2"), _request({"preset": "hex4", "side_mm": 40.0, "flat_to_flat_mm": 69.0, "gap_mm": 0.0})
        )
