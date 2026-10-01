"""Part 3: full-array export (layout files always; STEP export gated behind the optional `cad`
extra, and tested at a realistic-but-small array size, reporting file size/time per the brief)."""

import json
import pathlib

import pytest

from src.export import array_export as ae
from src.materials import Metal, Substrate

pytest.importorskip("cadquery", reason="STEP export needs the optional cad extra")


def _grid(n=4, length_m=2.0e-3, width_m=2.0e-3):
    return [[(length_m, width_m) for _ in range(n)] for _ in range(n)]


@pytest.mark.cad
def test_export_array_writes_layout_readme_and_step(tmp_path):
    substrate = Substrate("RT5880", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    metal = Metal("pec")
    result = ae.export_array(_grid(), 8.0e-3, 8.0e-3, substrate, metal, tmp_path)
    assert result.n_cells_per_side == 4
    assert result.total_cells == 16
    assert result.step_path is not None
    assert result.file_size_bytes > 0
    assert result.export_time_s >= 0.0

    layout = json.loads((tmp_path / "layout.json").read_text())
    assert layout["n_cells_per_side"] == 4
    assert layout["total_cells"] == 16
    assert len(layout["cells"]) == 16
    assert (tmp_path / "HFSS_README.md").exists()


@pytest.mark.cad
def test_non_square_grid_raises_clear_error(tmp_path):
    substrate = Substrate("RT5880", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    metal = Metal("pec")
    bad_grid = [[(2e-3, 2e-3)] * 4, [(2e-3, 2e-3)] * 3]  # ragged, not square
    with pytest.raises(ValueError, match="square"):
        ae.export_array(bad_grid, 8.0e-3, 8.0e-3, substrate, metal, tmp_path)


def test_layout_files_written_without_cad(tmp_path):
    """Layout/README must not require cadquery at all."""
    csv_path, json_path = ae.write_layout_files(_grid(n=2), tmp_path)
    assert csv_path.exists() and json_path.exists()
    readme_path = ae.write_hfss_readme(tmp_path)
    assert readme_path.exists()
    assert "Floquet port" in readme_path.read_text()


# ================================================================== PHASE15 Part C: feed markers


def test_readme_without_feeds_mm_is_byte_identical_to_before(tmp_path):
    with_none = ae.write_hfss_readme(tmp_path / "a").read_text()
    default = ae.write_hfss_readme(tmp_path / "b", feeds_mm=None).read_text()
    assert with_none == default == ae._HFSS_README_TEMPLATE


def test_readme_with_feeds_mm_lists_each_feed_location(tmp_path):
    text = ae.write_hfss_readme(tmp_path, feeds_mm=[(-2.0, 0.0), (2.0, 0.5)]).read_text()
    assert text.startswith(ae._HFSS_README_TEMPLATE)
    assert "Feed 0: (-2.000, 0.000) mm" in text
    assert "Feed 1: (2.000, 0.500) mm" in text


@pytest.mark.cad
def test_step_export_without_feeds_mm_is_unaffected(tmp_path):
    # STEP files carry a cadquery-generated timestamp/GUID, so two exports of the SAME geometry
    # are not byte-identical -- size (no extra "Feeds" body/geometry) is the property that matters.
    substrate = Substrate("RT5880", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    path_a = tmp_path / "a.step"
    path_b = tmp_path / "b.step"
    ae.write_hfss_step(_grid(n=2), 4e-3, 4e-3, substrate.thickness, path_a)
    ae.write_hfss_step(_grid(n=2), 4e-3, 4e-3, substrate.thickness, path_b, feeds_mm=None)
    assert abs(path_a.stat().st_size - path_b.stat().st_size) < 100


@pytest.mark.cad
def test_step_export_with_feeds_mm_adds_a_feeds_body(tmp_path):
    substrate = Substrate("RT5880", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    no_feeds = tmp_path / "no_feeds.step"
    with_feeds = tmp_path / "with_feeds.step"
    ae.write_hfss_step(_grid(n=2), 4e-3, 4e-3, substrate.thickness, no_feeds)
    ae.write_hfss_step(_grid(n=2), 4e-3, 4e-3, substrate.thickness, with_feeds, feeds_mm=[(-1.0, 0.0), (1.0, 0.0)])
    assert with_feeds.stat().st_size > no_feeds.stat().st_size
    assert "Feeds" in with_feeds.read_text(errors="ignore")


@pytest.mark.cad
def test_export_array_threads_feeds_mm_through(tmp_path):
    substrate = Substrate("RT5880", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    metal = Metal("pec")
    result = ae.export_array(_grid(), 8.0e-3, 8.0e-3, substrate, metal, tmp_path, feeds_mm=[(-4.0, 0.0), (4.0, 0.0)])
    assert result.step_path is not None
    readme_text = pathlib.Path(result.readme_path).read_text()
    assert "Feed 0: (-4.000, 0.000) mm" in readme_text


# ============================================================== PHASE15 Part C: multi-port export


def _rect_grid(n_y, n_x, length_m=2.0e-3, width_m=2.0e-3):
    return [[(length_m, width_m) for _ in range(n_x)] for _ in range(n_y)]


def test_write_layout_dxf_rect_accepts_a_non_square_grid(tmp_path):
    path = ae.write_layout_dxf_rect(_rect_grid(2, 4), 8.0e-3, 8.0e-3, tmp_path)
    assert path.exists()
    assert "PATCHES" in path.read_text()


def test_write_layout_dxf_still_rejects_non_square():
    with pytest.raises(ValueError, match="square"):
        ae.write_layout_dxf([[(2e-3, 2e-3)] * 4, [(2e-3, 2e-3)] * 3], 8e-3, 8e-3, pathlib.Path("/tmp"))


def test_readme_with_ports_info_names_one_port_per_monopole_and_the_s_matrix(tmp_path):
    text = ae.write_hfss_readme(
        tmp_path, feeds_mm=[(-2.0, 0.0), (2.0, 0.0)], ports_info=[(0, -2.0, 0.0), (1, 2.0, 0.0)]
    ).read_text()
    assert "Port 0: monopole at (-2.000, 0.000) mm" in text
    assert "S-parameter matrix" in text


# ===================================================== PHASE15 Part E: annular-slot coded cells


def _mixed_grid():
    annular = {"kind": "annular_slot", "a_m": 2e-3, "b_m": 2e-3, "R_m": 0.7e-3, "r_m": 0.4e-3}
    plain = (2e-3, 2e-3)
    return [[plain, annular], [annular, plain]]


def test_write_layout_files_reports_kind_and_slot_radii(tmp_path):
    csv_path, json_path = ae.write_layout_files(_mixed_grid(), tmp_path)
    text = csv_path.read_text()
    assert "kind" in text.splitlines()[0]
    assert "annular_slot" in text and "rectangle" in text
    layout = json.loads(json_path.read_text())
    kinds = {c["kind"] for c in layout["cells"]}
    assert kinds == {"rectangle", "annular_slot"}
    annular_cell = next(c for c in layout["cells"] if c["kind"] == "annular_slot")
    assert annular_cell["R_m"] == pytest.approx(0.7e-3)
    assert annular_cell["r_m"] == pytest.approx(0.4e-3)


def test_write_layout_dxf_rect_draws_circles_for_annular_slot_cells(tmp_path):
    path = ae.write_layout_dxf_rect(_mixed_grid(), 4e-3, 4e-3, tmp_path)
    text = path.read_text()
    assert text.count("CIRCLE") == 4  # 2 annular cells x (outer + inner) circle each


@pytest.mark.cad
def test_write_hfss_step_builds_extra_faces_for_annular_slot_cells(tmp_path):
    substrate_thickness_m = 1.575e-3
    rect_only = tmp_path / "rect_only.step"
    mixed = tmp_path / "mixed.step"
    rect_grid = [[(2e-3, 2e-3), (2e-3, 2e-3)], [(2e-3, 2e-3), (2e-3, 2e-3)]]
    info_rect = ae.write_hfss_step(rect_grid, 4e-3, 4e-3, substrate_thickness_m, rect_only)
    info_mixed = ae.write_hfss_step(_mixed_grid(), 4e-3, 4e-3, substrate_thickness_m, mixed)
    assert info_rect["n_patches"] == 4
    # 2 plain patches (1 face each) + 2 annular-slot cells (2 faces each: outer ring + inner disk)
    assert info_mixed["n_patches"] == 6
    assert mixed.stat().st_size > rect_only.stat().st_size


@pytest.mark.cad
def test_export_multi_port_array_handles_a_rectangular_board(tmp_path):
    substrate = Substrate("RT5880", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    metal = Metal("pec")
    grid = _rect_grid(2, 4)  # 2 rows x 4 columns: a rectangular, non-square board
    result = ae.export_multi_port_array(
        grid, 8.0e-3, 8.0e-3, substrate, metal, tmp_path,
        feeds_mm=[(-12.0, 0.0), (12.0, 0.0)], ports_info=[(0, -12.0, 0.0), (1, 12.0, 0.0)],
    )
    assert result.n_cells_x == 4
    assert result.n_cells_y == 2
    assert result.total_cells == 8
    assert result.step_path is not None
    assert result.file_size_bytes > 0
    readme_text = pathlib.Path(result.readme_path).read_text()
    assert "Port 0: monopole at (-12.000, 0.000) mm" in readme_text
    layout = json.loads(pathlib.Path(result.layout_json_path).read_text())
    assert len(layout["cells"]) == 8
