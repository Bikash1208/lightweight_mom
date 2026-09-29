"""Sec. 9-10: control-response derivations for candidates A-D."""

import numpy as np
import pytest

from src.phase8b.holographic_modulation import phase_difference_hologram
from src.phase8e.control_response import (
    generalized_transfer_operator,
    required_phase_b,
    required_phase_c,
    required_reactance_a,
)


def test_required_phase_b_matches_phase8b_kernel_directly():
    e_ref = np.array([1 + 0j, 1j, -1 + 0j])
    e_obj = np.array([1j, 1 + 0j, 1j])
    result = required_phase_b(e_ref, e_obj)
    assert np.allclose(result.values, phase_difference_hologram(e_ref, e_obj))
    assert result.value_kind == "phase_deg"


def test_required_phase_b_and_c_share_the_identical_numeric_kernel():
    e_ref = np.array([1 + 0j, 0.5j])
    e_obj = np.array([0.3j, 2 + 1j])
    b = required_phase_b(e_ref, e_obj)
    c = required_phase_c(e_ref, e_obj)
    assert np.allclose(b.values, c.values)
    assert b.required_quantity != c.required_quantity   # different physical claim, same numbers


def test_required_reactance_a_rejects_out_of_range_modulation_index():
    e_ref, e_obj = np.array([1 + 0j]), np.array([1 + 0j])
    with pytest.raises(ValueError):
        required_reactance_a(e_ref, e_obj, x_avg_ohm=300.0, modulation_index=1.5)


def test_generalized_transfer_operator_exactly_solves_t_dot_e_ref_equals_e_obj():
    rng = np.random.default_rng(0)
    ex_ref = rng.normal(size=5) + 1j * rng.normal(size=5)
    ey_ref = rng.normal(size=5) + 1j * rng.normal(size=5)
    ex_obj = rng.normal(size=5) + 1j * rng.normal(size=5)
    ey_obj = rng.normal(size=5) + 1j * rng.normal(size=5)
    result = generalized_transfer_operator(ex_ref, ey_ref, ex_obj, ey_obj)
    v = result.values
    recon_x = v["Txx"] * ex_ref + v["Txy"] * ey_ref
    recon_y = v["Tyx"] * ex_ref + v["Tyy"] * ey_ref
    assert np.allclose(recon_x, ex_obj, atol=1e-9)
    assert np.allclose(recon_y, ey_obj, atol=1e-9)


def test_generalized_transfer_operator_reduces_to_scalar_ratio_when_reference_has_no_y_component():
    ex_ref = np.array([2.0 + 0j])
    ey_ref = np.array([0.0 + 0j])
    ex_obj = np.array([1.0 + 1j])
    ey_obj = np.array([0.0 + 0j])
    result = generalized_transfer_operator(ex_ref, ey_ref, ex_obj, ey_obj)
    assert np.allclose(result.values["Txx"], ex_obj / ex_ref)
    assert np.allclose(result.values["Txy"], 0.0)
    assert np.allclose(result.values["Tyx"], 0.0)


def test_generalized_transfer_operator_flags_undefined_where_reference_vanishes():
    ex_ref = np.array([0.0 + 0j, 1.0 + 0j])
    ey_ref = np.array([0.0 + 0j, 0.0 + 0j])
    ex_obj = np.array([1.0 + 0j, 1.0 + 0j])
    ey_obj = np.array([0.0 + 0j, 0.0 + 0j])
    result = generalized_transfer_operator(ex_ref, ey_ref, ex_obj, ey_obj)
    assert result.values["undefined_mask"][0] == True    # noqa: E712
    assert result.values["undefined_mask"][1] == False   # noqa: E712
