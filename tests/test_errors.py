"""Phase 13 Part 2: the typed error-code taxonomy."""

import pytest

from src import errors


def test_every_error_code_is_unique():
    codes = [e.code for e in errors.ALL_ERROR_CODES]
    assert len(codes) == len(set(codes))


@pytest.mark.parametrize(
    "factory, args",
    [
        (errors.config_invalid, ()),
        (errors.geometry_out_of_cell, ()),
        (errors.mesh_unaligned, ()),
        (errors.floquet_rank_deficient, ()),
        (errors.solve_failed, ()),
        (errors.not_converged, ()),
        (errors.phase_unreachable, ()),
        (errors.cad_unavailable, ()),
        (errors.impedance_unreachable, ()),
    ],
)
def test_factory_functions_raise_lwmom_error_with_stable_code(factory, args):
    err = factory(*args)
    assert isinstance(err, errors.LwmomError)
    assert err.code in {e.code for e in errors.ALL_ERROR_CODES}
    assert err.message and err.remedy


def test_wrap_preserves_original_exception_text():
    original = ValueError("period_x: metres expected")
    err = errors.config_invalid(original)
    assert err.code == "CONFIG_INVALID"
    assert "period_x" in err.debug_detail


def test_field_expression_unsupported_names_the_bad_expression():
    err = errors.field_expression_unsupported("Foo(Bar)")
    assert "Foo(Bar)" in err.message


def test_to_dict_has_every_field():
    err = errors.solve_failed(RuntimeError("singular matrix"))
    d = err.to_dict()
    assert set(d) == {"code", "message", "remedy", "debug_detail"}
    assert d["debug_detail"] == "singular matrix"


def test_no_silent_fallback_error_is_actually_raised():
    with pytest.raises(errors.LwmomError) as exc_info:
        raise errors.not_converged("fit failed, non-monotonic")
    assert exc_info.value.code == "NOT_CONVERGED"
