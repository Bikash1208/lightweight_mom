"""Phase 12C: the working design point, loaded from config_phase12c_grounded_8mm.yaml."""

from src.phase12c import design_point as dp


def test_design_point_is_grounded_8mm():
    assert dp.BACKING == "pec"
    assert dp.PERIOD_M == 8.0e-3
    assert dp.SUBSTRATE.eps_r == 2.2
    assert dp.SUBSTRATE.thickness == 1.575e-3
    assert dp.DESIGN_FREQUENCY_HZ == 12.5e9


def test_reference_planes_have_no_bottom_port():
    assert dp.REFERENCE_PLANES.bottom_air_thickness is None
    assert dp.REFERENCE_PLANES.grounded is True


def test_unit_cell_at_varies_only_geometry():
    a = dp.unit_cell_at(2.0e-3)
    b = dp.unit_cell_at(3.0e-3, 4.0e-3)
    assert a.geometry.length == a.geometry.width == 2.0e-3
    assert (b.geometry.length, b.geometry.width) == (3.0e-3, 4.0e-3)
    assert a.lattice is dp.LATTICE and b.lattice is dp.LATTICE
    assert a.substrate is dp.SUBSTRATE and a.metal is dp.METAL
    assert a.backing == "pec" == b.backing


def test_config_file_round_trips_through_the_real_loader():
    from src.config import load_config

    cfg = load_config(dp.CONFIG_PATH)
    assert cfg.unit_cell.backing == "pec"
    assert cfg.unit_cell.lattice.period_x == dp.PERIOD_M
