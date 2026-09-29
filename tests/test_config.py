import copy
import math
from pathlib import Path

import pytest
import yaml

from src.config import ConfigError, config_from_dict, load_config

CONFIG_YAML = Path(__file__).resolve().parent.parent / "configs" / "config.yaml"


def base():
    return {
        "period_x": 8.0e-3,
        "period_y": 8.0e-3,
        "substrate": {"name": "sub", "eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3},
        "metal": {"model": "pec"},
        "geometry": {"type": "rectangle", "length": 5.0e-3, "width": 2.0e-3},
        "excitation": {"frequency": 10.0e9},
        "mesh": {"cells_x": 16, "cells_y": 16, "floquet_order": 32},
        "reference_planes": {"top_air_thickness": 10.0e-3, "bottom_air_thickness": 10.0e-3},
    }


def modified(**changes):
    d = base()
    for path, value in changes.items():
        keys = path.split("__")
        node = d
        for k in keys[:-1]:
            node = node[k]
        if value is KeyError:
            del node[keys[-1]]
        else:
            node[keys[-1]] = value
    return d


def test_shipped_config_yaml_loads_and_converts_to_si():
    cfg = load_config(CONFIG_YAML)
    cell, exc = cfg.unit_cell, cfg.excitation
    assert cell.lattice.period_x == 8e-3 and cell.lattice.period_y == 8e-3
    assert cell.substrate.eps_r == 2.2
    assert cell.substrate.tan_delta == 0.0009
    assert cell.substrate.thickness == 1.575e-3
    assert cell.geometry.kind == "rectangle"
    assert cell.geometry.length == 5e-3 and cell.geometry.width == 2e-3
    assert exc.frequency == 10e9
    assert exc.theta == 0.0
    assert (cfg.mesh.cells_x, cfg.mesh.cells_y, cfg.mesh.floquet_order) == (16, 16, 32)
    rp = cfg.reference_planes
    assert rp.top_air_thickness == 10e-3 and rp.bottom_air_thickness == 10e-3
    assert rp.z_top == 10e-3
    assert rp.z_bottom == pytest.approx(-(1.575e-3 + 10e-3), rel=1e-14)


def test_pyyaml_reads_8e_minus_3_as_string_and_loader_copes():
    assert isinstance(yaml.safe_load("x: 8e-3")["x"], str)  # documents the PyYAML gotcha
    d = yaml.safe_load(
        "period_x: 8e-3\nperiod_y: 8e-3\n"
        "substrate: {eps_r: 2.2, tan_delta: 0.0009, thickness: 1.575e-3}\n"
        "geometry: {type: rectangle, length: 5e-3, width: 2e-3}\n"
        "excitation: {frequency: 10e9}\n"
        "mesh: {cells_x: 16, cells_y: 16, floquet_order: 32}\n"
        "reference_planes: {top_air_thickness: 10e-3, bottom_air_thickness: 10e-3}\n"
    )
    cfg = config_from_dict(d)
    assert cfg.unit_cell.lattice.period_x == 8e-3
    assert cfg.unit_cell.geometry.length == 5e-3
    assert cfg.excitation.frequency == 10e9


def test_unit_suffixes_convert_at_input_boundary():
    d = modified(
        period_x="8 mm",
        period_y="0.8 cm",
        substrate__thickness="1.575 mm",
        geometry__length="5000 um",
        excitation__frequency="10 GHz",
    )
    cfg = config_from_dict(d)
    assert cfg.unit_cell.lattice.period_x == pytest.approx(8e-3, rel=1e-14)
    assert cfg.unit_cell.lattice.period_y == pytest.approx(8e-3, rel=1e-14)
    assert cfg.unit_cell.substrate.thickness == pytest.approx(1.575e-3, rel=1e-14)
    assert cfg.unit_cell.geometry.length == pytest.approx(5e-3, rel=1e-14)
    assert cfg.excitation.frequency == pytest.approx(10e9, rel=1e-14)


def test_angles_in_degrees_become_radians():
    d = modified(geometry__rotation_deg=90.0, excitation__phi_deg=90.0, excitation__polarization="x")
    cfg = config_from_dict(d)
    assert cfg.unit_cell.geometry.rotation == pytest.approx(math.pi / 2, rel=1e-15)
    assert cfg.excitation.phi == pytest.approx(math.pi / 2, rel=1e-15)


def test_reference_planes_accept_units_and_zero_thickness():
    d = modified(reference_planes__top_air_thickness="12 mm", reference_planes__bottom_air_thickness=0.0)
    rp = config_from_dict(d).reference_planes
    assert rp.z_top == pytest.approx(12e-3, rel=1e-14)
    assert rp.z_bottom == pytest.approx(-1.575e-3, rel=1e-14)  # plane on the substrate bottom face


def test_optional_fields_default_sensibly():
    cfg = config_from_dict(base())
    g = cfg.unit_cell.geometry
    assert (g.center_x, g.center_y, g.rotation) == (0.0, 0.0, 0.0)
    assert cfg.excitation.polarization == "x"
    assert cfg.unit_cell.metal.model == "pec"
    assert config_from_dict({k: v for k, v in base().items() if k != "metal"}).unit_cell.metal.model == "pec"


def test_loader_does_not_mutate_input():
    d = base()
    snapshot = copy.deepcopy(d)
    config_from_dict(d)
    assert d == snapshot


@pytest.mark.parametrize(
    "changes, match",
    [
        # no electrical defaults: every substrate value must be given
        ({"substrate__eps_r": KeyError}, "eps_r"),
        ({"substrate__tan_delta": KeyError}, "tan_delta"),
        ({"substrate__thickness": KeyError}, "thickness"),
        ({"excitation__frequency": KeyError}, "frequency"),
        ({"period_x": KeyError}, "period_x"),
        # typos must not be silently ignored
        ({"periodx": 8e-3}, "unknown"),
        ({"substrate__epsr": 2.2}, "unknown"),
        ({"geometry__lenght": 5e-3}, "unknown"),
        # unit mistakes
        ({"period_x": 8}, "metres"),
        ({"period_x": "8 Mm"}, "unit"),
        ({"excitation__frequency": "10 ghz"}, "unit"),
        ({"excitation__frequency": -1e9}, "frequency"),
        ({"excitation__frequency": 0}, "frequency"),
        ({"period_x": float("nan")}, "finite"),
        # physical validity
        ({"substrate__eps_r": -2.2}, "eps_r"),
        ({"substrate__tan_delta": -0.1}, "tan_delta"),
        ({"substrate__thickness": 0.0}, "thickness"),
        ({"geometry__type": "circle"}, "geometry type"),
        ({"geometry__length": 9e-3}, "outside the unit cell"),
        ({"geometry__length": -5e-3}, "geometry"),
        # Phase-3 sections: mesh and explicit reference planes are mandatory
        ({"mesh": KeyError}, "mesh"),
        ({"reference_planes": KeyError}, "reference_planes"),
        ({"mesh__cells_x": KeyError}, "cells_x"),
        ({"mesh__floquet_order": KeyError}, "floquet_order"),
        ({"mesh__cells_x": 1}, "cells_x"),
        ({"mesh__cells_y": 16.0}, "integer"),
        ({"mesh__cells_y": True}, "integer"),
        ({"mesh__floquet_order": -1}, "floquet_order"),
        ({"mesh__nx": 8}, "unknown"),
        ({"reference_planes__top_air_thickness": KeyError}, "top_air_thickness"),
        ({"reference_planes__bottom_air_thickness": -1e-3}, ">= 0"),
        ({"reference_planes__top_air_thickness": 10}, "metres"),
        ({"reference_planes__height": 1e-3}, "unknown"),
        # Phase-1 scope gates
        ({"metal__model": "copper"}, "pec"),
        ({"excitation__theta_deg": 30.0}, "normal incidence"),
        ({"excitation__polarization": "z"}, "polarization"),
    ],
)
def test_invalid_configs_are_rejected(changes, match):
    with pytest.raises(ConfigError, match=match):
        config_from_dict(modified(**changes))


def test_root_must_be_mapping():
    for bad in (None, [], "text", 5):
        with pytest.raises(ConfigError):
            config_from_dict(bad)


def test_load_config_error_paths(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_config(tmp_path / "missing.yaml")
    bad = tmp_path / "bad.yaml"
    bad.write_text("period_x: [unclosed\n")
    with pytest.raises(ConfigError, match="invalid YAML"):
        load_config(bad)
    empty = tmp_path / "empty.yaml"
    empty.write_text("")
    with pytest.raises(ConfigError):
        load_config(empty)


def test_config_error_is_value_error():
    assert issubclass(ConfigError, ValueError)
