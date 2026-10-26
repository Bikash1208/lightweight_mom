"""Phase 20 R1 item 1: oblique incidence in the TE/TM basis -- bare-slab Fresnel/Airy, energy conservation, reciprocity,
symmetry, the theta = 0 limit, and the bit-identical normal-incidence gate."""

import json
import math
import warnings
from pathlib import Path

import numpy as np
import pytest

from src.basis_functions import MeshSettings
from src.constants import MU0
from src.excitation import PlaneWaveExcitation, ReferencePlanes
from src.geometry import Geometry, UnitCell
from src.green_function import LayeredGreenFunction
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.sparameters import (
    TE,
    TM,
    ObliqueSParameters,
    compute_s_parameters,
    compute_s_parameters_oblique,
    oblique_stack_s_parameters,
    s_parameters_from_config,
    sweep_s_parameters_oblique,
)
from tests._platform_tol import ON_REFERENCE_PLATFORM, assert_matches_reference

F = 10e9
SUB = Substrate("rt5880", 2.2, 0.0, 1.575e-3)
LOSSY = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
ROOT = Path(__file__).resolve().parents[1]


def _planes(backing, top=2e-3):
    return ReferencePlanes(top, None if backing == "pec" else 2e-3, SUB.thickness)


def _cell(backing, size=5e-3, period=8e-3, sub=SUB, rotation=0.0):
    return UnitCell(
        RectangularLattice(period, period),
        Geometry(kind="rectangle", length=size, width=size, rotation=rotation),
        sub,
        Metal(),
        backing=backing,
    )


def _airy(theta, backing, top=2e-3):
    """Analytic reflection (and transmission) of the bare slab, power-normalised, TE/TM, referenced to z_top."""
    g = LayeredGreenFunction(SUB, F, backing=backing)
    kt = g.k0 * math.sin(theta)
    a, b = g.normal_wavenumbers(kt, 0.0)
    w, eps, h = g.omega, SUB.eps_r, SUB.thickness
    out = {}
    for mode, z_air, z_sub in (
        ("TE", w * MU0 / a, w * MU0 / b),
        ("TM", a / (w * 8.8541878128e-12), b / (w * 8.8541878128e-12 * eps)),
    ):
        if backing == "pec":
            z_in = 1j * z_sub * np.tan(b * h)
            r = (z_in - z_air) / (z_in + z_air)
            out[mode] = (r * np.exp(-2j * a * top), None)
        else:
            z_l = z_air  # air below
            z_in = z_sub * (z_l + 1j * z_sub * np.tan(b * h)) / (z_sub + 1j * z_l * np.tan(b * h))
            r = (z_in - z_air) / (z_in + z_air)
            # transmission: energy-based magnitude, phase from the TL cascade
            out[mode] = (r * np.exp(-2j * a * top), None)
    return out


@pytest.mark.parametrize("backing", ["pec", "air"])
@pytest.mark.parametrize("theta_deg", [0.0, 30.0, 60.0])
@pytest.mark.parametrize("phi_deg", [0.0, 37.0])
def test_bare_slab_reflection_equals_fresnel_airy_te_and_tm(backing, theta_deg, phi_deg):
    theta, phi = math.radians(theta_deg), math.radians(phi_deg)
    g = LayeredGreenFunction(SUB, F, backing=backing)
    s = oblique_stack_s_parameters(g, _planes(backing), theta, phi)
    want = _airy(theta, backing)
    assert abs(s.s11[TE, TE] - want["TE"][0]) < 1e-10 and abs(s.s11[TM, TM] - want["TM"][0]) < 1e-10
    assert abs(s.s11[TM, TE]) < 1e-12 and abs(s.s11[TE, TM]) < 1e-12  # no polarisation conversion without metal
    if backing == "pec":
        assert abs(abs(s.s11[TE, TE]) - 1.0) < 1e-12  # lossless grounded: total reflection
    else:
        p = s.power_out_per_incident()
        assert np.allclose(p, 1.0, atol=1e-12)  # |S11|^2 + |S21|^2 = 1 per mode
    assert s.incoming_residual < 1e-10


@pytest.mark.parametrize("theta_deg", [0.0, 30.0, 60.0])
@pytest.mark.parametrize("backing", ["pec", "air"])
def test_energy_conservation_with_a_patch_at_oblique_incidence(theta_deg, backing):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = compute_s_parameters_oblique(
            _cell(backing), F, MeshSettings(16, 16, 10), _planes(backing), math.radians(theta_deg), math.radians(25.0)
        )
    assert not s.higher_orders_propagate
    assert np.allclose(s.power_out_per_incident(), 1.0, atol=2e-9)
    assert s.incoming_residual < 1e-8


@pytest.mark.parametrize("backing", ["pec", "air"])
def test_phi_plus_180_for_a_centrosymmetric_cell_and_square_patch_90_degree_rotation(backing):
    theta = math.radians(40.0)
    mesh = MeshSettings(16, 16, 10)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = compute_s_parameters_oblique(_cell(backing), F, mesh, _planes(backing), theta, 0.0)
        b = compute_s_parameters_oblique(_cell(backing), F, mesh, _planes(backing), theta, math.pi)
        c = compute_s_parameters_oblique(_cell(backing), F, mesh, _planes(backing), theta, math.pi / 2)
    assert np.allclose(a.s11, b.s11, atol=1e-9)  # (theta, phi) <-> (theta, phi + 180): centrosymmetric cell
    assert np.allclose(
        a.s11, c.s11, atol=1e-9
    )  # square patch on a square lattice: TE/TM at phi = 0 equal those at phi = 90
    if backing == "air":
        assert np.allclose(a.s21, c.s21, atol=1e-9)


def test_theta_zero_te_tm_are_y_and_x_of_the_normal_incidence_solver():
    mesh = MeshSettings(16, 16, 10)
    uc = UnitCell(
        RectangularLattice(10e-3, 8e-3),
        Geometry(kind="rectangle", length=6e-3, width=5e-3),
        LOSSY,
        Metal(),
        backing="pec",
    )
    pl = ReferencePlanes(2e-3, None, LOSSY.thickness)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        old = compute_s_parameters(uc, F, mesh, pl)
        new = compute_s_parameters_oblique(uc, F, mesh, pl, 0.0, 0.0)
    # [received, incident]: TE = y (index 1), TM = x (index 0)
    assert abs(new.s11[TE, TE] - old.s11[1, 1]) < 1e-9 and abs(new.s11[TM, TM] - old.s11[0, 0]) < 1e-9
    assert abs(new.s11[TM, TE] - old.s11[0, 1]) < 1e-9 and abs(new.s11[TE, TM] - old.s11[1, 0]) < 1e-9


def test_ungrounded_theta_zero_matches_the_normal_incidence_solver_for_all_four_blocks():
    mesh = MeshSettings(16, 16, 10)
    uc = _cell("air", sub=LOSSY)
    pl = ReferencePlanes(2e-3, 2e-3, LOSSY.thickness)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        old = compute_s_parameters(uc, F, mesh, pl)
        new = compute_s_parameters_oblique(uc, F, mesh, pl, 0.0, 0.0)
    idx = {TE: 1, TM: 0}
    for fam in ("s11", "s12", "s21", "s22"):
        for r in (TE, TM):
            for i in (TE, TM):
                assert abs(getattr(new, fam)[r, i] - getattr(old, fam)[idx[r], idx[i]]) < 1e-9


def test_reference_plane_shift_uses_kz0_and_is_exact():
    theta = math.radians(50.0)
    mesh = MeshSettings(16, 16, 10)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = compute_s_parameters_oblique(_cell("pec"), F, mesh, _planes("pec", 2e-3), theta, 0.2)
        b = compute_s_parameters_oblique(_cell("pec"), F, mesh, _planes("pec", 5e-3), theta, 0.2)
    assert np.allclose(a.at_planes(_planes("pec", 5e-3)).s11, b.s11, atol=1e-9)


def test_grating_order_warning_and_record_at_an_angle_that_makes_a_higher_order_propagate():
    uc = _cell("pec", size=8e-3, period=25e-3)  # 25 mm period at 10 GHz: orders propagate even at normal incidence
    with pytest.warns(UserWarning, match="higher Floquet order"):
        s = compute_s_parameters_oblique(uc, F, MeshSettings(16, 16, 10), _planes("pec"), math.radians(20.0), 0.0)
    assert s.higher_orders_propagate and s.to_record()["higher_orders_propagate"] is True
    res, recs = sweep_s_parameters_oblique(
        _cell("pec"), [9e9, 10e9], MeshSettings(16, 16, 10), _planes("pec"), math.radians(30.0), 0.0
    )
    assert (
        len(res) == 2
        and all(isinstance(r, ObliqueSParameters) for r in res)
        and recs[0]["higher_orders_propagate"] is False
    )


def test_config_allows_oblique_incidence_and_dispatches():
    from src.config import config_from_dict

    cfg = config_from_dict(
        {
            "period_x": 8.0e-3,
            "period_y": 8.0e-3,
            "substrate": {"name": "sub", "eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3, "backing": "pec"},
            "metal": {"model": "pec"},
            "geometry": {"type": "rectangle", "length": 5.0e-3, "width": 5.0e-3},
            "excitation": {"frequency": 10.0e9, "theta_deg": 30.0, "phi_deg": 10.0, "polarization": "te"},
            "mesh": {"cells_x": 16, "cells_y": 16, "floquet_order": 10},
            "reference_planes": {"top_air_thickness": 10.0e-3},
        }
    )
    assert math.isclose(cfg.excitation.theta, math.radians(30.0)) and cfg.excitation.polarization == "te"
    assert np.isclose(np.dot(cfg.excitation.e_field_direction, np.array(cfg.excitation.wavevector)), 0.0, atol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = s_parameters_from_config(cfg)
    assert isinstance(s, ObliqueSParameters) and math.isclose(s.theta, cfg.excitation.theta)


def test_te_tm_excitation_vectors_are_perpendicular_to_k_and_match_the_x_y_labels_at_normal_incidence():
    for pol in ("te", "tm"):
        e = PlaneWaveExcitation(F, theta=math.radians(35.0), phi=math.radians(70.0), polarization=pol)
        assert (
            abs(np.dot(e.e_field_direction, e.wavevector)) < 1e-9
            and abs(np.linalg.norm(e.e_field_direction) - 1.0) < 1e-12
        )
    assert np.allclose(PlaneWaveExcitation(F, polarization="te").e_field_direction, [0, 1, 0])
    assert np.allclose(PlaneWaveExcitation(F, polarization="tm").e_field_direction, [1, 0, 0])


# ------------------------------------------------------------------------------ the hard gate
def test_normal_incidence_results_are_bit_identical_to_the_phase20_base_commit():
    """`tests/data/r1_theta0_reference.json` was written at the base commit of branch phase20-ris (hex floats)."""
    from tests.data._gen_r1_theta0_reference import compute

    ref = json.loads((ROOT / "tests/data/r1_theta0_reference.json").read_text())
    got = compute()
    assert got.keys() == ref.keys()

    def unhex(o):
        if isinstance(o, str):
            return float.fromhex(o)
        return [unhex(x) for x in o]

    for key in ref:
        for blk in ref[key]:
            g, r = json.loads(json.dumps(got[key][blk])), ref[key][blk]
            if ON_REFERENCE_PLATFORM:
                assert g == r, (key, blk)  # hex strings: bit for bit
            else:
                assert_matches_reference(unhex(g), unhex(r), f"{key}/{blk}")
