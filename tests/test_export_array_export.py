"""Part 3: full-array export (layout files always; STEP export gated behind the optional `cad`
extra, and tested at a realistic-but-small array size, reporting file size/time per the brief)."""

import json

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
