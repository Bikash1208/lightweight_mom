"""Phase 27 item 2: the cell-family registry -- domains of validity, named unavailable states, N-state RIS schema."""

import json

import pytest

from src.hologram import cell_family_registry as reg


def test_tensor_family_domain_and_named_unavailable_state():
    ok = reg.check_request(
        "tensor_slotted_circle_rt5880_6p4", substrate=(2.2, 0.0009, 1.575), f_ghz=10.54, period_mm=6.4
    )
    assert (
        ok.state == reg.AVAILABLE
        and ok.validation_label == "SOLVER_SUPERCELL"
        and any("SLL optimistic" in b for b in ok.known_biases)
    )
    bad = reg.check_request("tensor_slotted_circle_rt5880_6p4", substrate=(2.2, 0.0009, 3.175), f_ghz=10.2)
    assert bad.state == reg.NOT_AVAILABLE == "NOT_AVAILABLE_FROM_CURRENT_DATA"
    assert len(bad.reasons) == 2  # substrate and frequency, both reported
    assert (
        reg.check_request("tensor_slotted_circle_rt5880_6p4", theta_deg=30.0).state == reg.NOT_AVAILABLE
    )  # normal-incidence table only


def test_band_family_accepts_inside_and_refuses_outside_band():
    assert reg.check_request("square_patch_grounded_rt5880_8mm", f_ghz=12.5).state == reg.AVAILABLE
    assert reg.check_request("square_patch_grounded_rt5880_8mm", f_ghz=9.0).state == reg.NOT_AVAILABLE


def test_builtin_families_are_schema_valid_and_json_serialisable():
    fams = reg.list_families()
    assert {
        "tensor_slotted_circle_rt5880_6p4",
        "onebit_empty_square585_rt5880_6p4",
        "square_patch_grounded_rt5880_8mm",
    } <= set(fams)
    for f in fams.values():
        reg.validate_family(f)
        json.dumps(f.to_dict())
    one_bit = fams["onebit_empty_square585_rt5880_6p4"]
    assert (
        one_bit.data_kind == "discrete_states"
        and len(one_bit.states) == 2
        and one_bit.validation_label == "NOT_HFSS_VALIDATED"
    )


def _ris_family(n_states=4, f=(9.0, 10.0, 11.0)):
    grid = reg.ResponseGrid(
        "reflection_2x2",
        f,
        theta_deg=(0.0, 30.0, 60.0),
        phi_deg=(0.0, 90.0),
        basis="linear_xy",
        data_path="ris/state.npz",
    )
    states = tuple(
        reg.StateSpec(f"s{k}", f"{90 * k} deg", grid, bias={"d1": k & 1, "d2": (k >> 1) & 1}) for k in range(n_states)
    )
    return reg.CellFamily(
        name="synthetic_2bit_ris",
        geometry_type="synthetic",
        substrate=reg.SubstrateSpec("x", 3.0, 0.002, 1.5),
        backing="pec",
        period_mm=10.0,
        data_kind="discrete_states",
        validation_label="MODEL_ONLY",
        frequencies_ghz=f,
        incidence_max_theta_deg=60.0,
        states=states,
    )


def test_synthetic_four_state_ris_family_registers_without_schema_change():
    fam = reg.register(_ris_family())
    try:
        assert len(reg.get("synthetic_2bit_ris").states) == 4
        assert (
            reg.check_request("synthetic_2bit_ris", f_ghz=10.0, theta_deg=45.0, period_mm=10.0).state == reg.AVAILABLE
        )
        out = reg.check_request("synthetic_2bit_ris", f_ghz=10.0, theta_deg=75.0)
        assert out.state == reg.NOT_AVAILABLE and "incidence" in out.reasons[0]
        assert json.loads(json.dumps(fam.to_dict()))["states"][3]["bias"] == {"d1": 1, "d2": 1}
    finally:
        reg.unregister("synthetic_2bit_ris")


def test_malformed_entries_are_rejected():
    fam = _ris_family()
    with pytest.raises(ValueError):
        reg.validate_family(_ris_family(n_states=1))
    with pytest.raises(ValueError):  # response grid that misses a family frequency
        reg.validate_family(
            _ris_family(f=(9.0, 10.0, 11.0)).__class__(**{**fam.__dict__, "frequencies_ghz": (9.0, 10.0, 11.0, 12.0)})
        )
    with pytest.raises(ValueError):
        reg.validate_family(reg.CellFamily(**{**fam.__dict__, "validation_label": "TRUST_ME"}))
    with pytest.raises(ValueError):  # states on a continuous family
        reg.validate_family(reg.CellFamily(**{**fam.__dict__, "data_kind": "scalar_gamma"}))
    with pytest.raises(KeyError):
        reg.check_request("no_such_family")


def test_ris_families_are_registered_model_only_discrete_states_with_their_domain_and_biases():
    cp = reg.get("ris_cp2_octagon_madp_rt5880_15")
    lp = reg.get("ris_lp_1bit_from_cp2_madp_rt5880_15")
    for fam in (cp, lp):
        reg.validate_family(fam)
        assert fam.data_kind == "discrete_states" and fam.validation_label == "MODEL_ONLY" and len(fam.states) == 4
        assert fam.period_mm == 15.0 and fam.incidence_max_theta_deg == 45.0 and fam.frequency_band_ghz == (8.0, 12.0)
        assert all(s.response.kind == "reflection_2x2" for s in fam.states)
        assert any("flip-chip" in b for b in fam.known_biases) and any("random codings" in b for b in fam.known_biases)
        assert any("EXTRAPOLATED_BEYOND_R3_TARGETS" in b for b in fam.known_biases) and fam.diffraction_check
        json.dumps(fam.to_dict())
    assert [s.bias["diode_1"] for s in cp.states] == ["on", "off", "off", "off"]  # one-hot
    assert reg.check_request(cp.name, f_ghz=10.1, theta_deg=30.0, period_mm=15.0).state == reg.AVAILABLE
    out = reg.check_request(cp.name, f_ghz=10.0, theta_deg=60.0, period_mm=6.4)
    assert out.state == reg.NOT_AVAILABLE and len(out.reasons) == 2  # 60 deg and the wrong period, both reported
    assert reg.check_request(cp.name, f_ghz=12.5).state == reg.NOT_AVAILABLE


def test_ris_family_refuses_requests_where_a_higher_diffraction_order_propagates():
    name = "ris_cp2_octagon_madp_rt5880_15"
    assert reg.check_request(name, f_ghz=11.0, theta_deg=45.0).state == reg.AVAILABLE  # below onset at every azimuth
    out = reg.check_request(name, f_ghz=11.9, theta_deg=45.0)
    assert out.state == reg.NOT_AVAILABLE and "diffraction order" in out.reasons[0]
    # with an explicit azimuth the verdict is for that azimuth (at theta 45 and 12 GHz the order propagates at azimuth 0 deg, not at 45 deg)
    assert reg.check_request(name, f_ghz=12.0, theta_deg=45.0, phi_deg=0.0).state == reg.NOT_AVAILABLE
    assert reg.check_request(name, f_ghz=12.0, theta_deg=45.0, phi_deg=45.0).state == reg.AVAILABLE
