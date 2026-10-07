"""Phase 12C Task E: HFSS-ready exports. Skips without the optional cadquery dependency."""

import pytest

cadquery = pytest.importorskip("cadquery")

from src.phase12c import exports  # noqa: E402


def test_export_single_cells_writes_a_step_file_per_l(tmp_path):
    # mesh_cells_per_side=4 is too coarse to resolve L=1.6 mm at P=8 mm (no mesh cell centre falls
    # inside the patch); use values that resolve at this deliberately small test mesh.
    results = exports.export_single_cells(tmp_path, l_values_mm=(3.2, 4.0), mesh_cells_per_side=4)
    assert len(results) == 2
    for r in results:
        assert (tmp_path / f"single_cell_L{r.l_mm:.1f}mm.step").exists()
        assert 0.0 < r.s11_mag <= 1.0 + 1e-9


def test_export_supercells_writes_uniform_gradient_and_strong_contrast(tmp_path):
    results = exports.export_supercells(tmp_path, mesh_cells_per_subcell=4)
    names = {r.name for r in results}
    assert names == {"uniform_3x3", "smooth_gradient_3x3", "strong_contrast_3x3"}
    for r in results:
        assert (tmp_path / f"{r.name}.step").exists()
        assert r.n_metal_solids == 9  # 3x3
        assert r.has_ground is True


def test_synthetic_hfss_reference_is_clearly_labelled_and_one_port(tmp_path):
    path, records = exports.write_synthetic_hfss_reference_test(tmp_path / "synthetic.csv")
    assert "SYNTHETIC" in records[0].cell_id
    assert all(r.s21 is None for r in records)


def test_write_export_manifest_records_the_synthetic_disclaimer(tmp_path):
    singles = exports.export_single_cells(tmp_path, l_values_mm=(3.2,), mesh_cells_per_side=4)
    supers = exports.export_supercells(tmp_path, mesh_cells_per_subcell=4)
    path, records = exports.write_synthetic_hfss_reference_test(tmp_path / "synthetic.csv")
    manifest_path = exports.write_export_manifest(singles, supers, path, tmp_path / "manifest.json")
    import json

    manifest = json.loads(manifest_path.read_text())
    assert "No real HFSS results" in manifest["note"]
