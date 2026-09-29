"""Sec. 22-23: provenance/validation-level taxonomy."""

import pytest

from src.phase8d_rev.validation_levels import (
    VALIDATION_LEVELS,
    Provenance,
    library_validation_status,
    record_validation_status,
)


def test_validation_levels_are_exactly_five_in_order():
    assert VALIDATION_LEVELS == (
        "DATA_PRESENT", "NUMERICALLY_VALIDATED", "CROSS_SOLVER_VALIDATED",
        "HFSS_VALIDATED", "EXPERIMENTALLY_VALIDATED",
    )


def test_provenance_rejects_unknown_validation_status():
    with pytest.raises(ValueError):
        Provenance("MoM", None, None, None, None, "SOMETHING_ELSE")


def test_provenance_accepts_every_known_level():
    for level in VALIDATION_LEVELS:
        Provenance("MoM", None, None, None, None, level)


def test_record_validation_status_numerically_validated_iff_checks_passed():
    assert record_validation_status({"checks_passed": True}) == "NUMERICALLY_VALIDATED"
    assert record_validation_status({"checks_passed": False}) == "DATA_PRESENT"
    assert record_validation_status({}) == "DATA_PRESENT"


class _FakeRecord:
    def __init__(self, checks_passed):
        self.diagnostics = {"checks_passed": checks_passed}


def test_library_validation_status_requires_every_record_to_pass():
    assert library_validation_status([_FakeRecord(True), _FakeRecord(True)]) == "NUMERICALLY_VALIDATED"
    assert library_validation_status([_FakeRecord(True), _FakeRecord(False)]) == "DATA_PRESENT"


def test_library_validation_status_rejects_empty_input():
    with pytest.raises(ValueError):
        library_validation_status([])
