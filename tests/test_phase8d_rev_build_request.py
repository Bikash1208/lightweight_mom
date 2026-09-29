"""Sec. 18-19, 21: LibraryBuildRequest generation from a missed query."""

import pytest

from src.phase8d_rev.build_request import generate_build_request
from src.phase8d_rev.query import QueryResult, UnitCellQuery


def _missing_result(**query_kwargs):
    query = UnitCellQuery(**query_kwargs)
    return QueryResult(
        status="DATABASE_COVERAGE_MISSING", query=query, matching_entries=(),
        filter_trail=(("initial", 1), ("topology", 1), ("period", 0)), eliminated_at="period",
    )


def test_generate_build_request_raises_on_a_matched_result():
    query = UnitCellQuery(topology="rectangular_patch")
    matched = QueryResult(status="MATCH", query=query, matching_entries=(object(),),
                           filter_trail=(("initial", 1),), eliminated_at=None)
    with pytest.raises(ValueError):
        generate_build_request(matched)


def test_generate_build_request_carries_the_query_fields_through():
    result = _missing_result(topology="rectangular_patch", p_x_mm=3.2, p_y_mm=3.2, frequency_hz=12.5e9,
                              polarization="x", incidence="normal")
    req = generate_build_request(result)
    assert req.topology == "rectangular_patch"
    assert req.p_x_mm == 3.2 and req.p_y_mm == 3.2
    assert req.polarization == "x"
    assert req.incidence == "normal"
    assert req.status == "BUILD_REQUIRED"
    assert "period" in req.reason


def test_generate_build_request_never_silently_fills_material_stack_unless_queried():
    result = _missing_result(topology="rectangular_patch")
    req = generate_build_request(result)
    assert req.material_stack == {}


def test_generate_build_request_carries_material_stack_when_queried():
    result = _missing_result(topology="rectangular_patch", substrate_name="Foo", eps_r=3.0)
    req = generate_build_request(result)
    assert req.material_stack["name"] == "Foo"
    assert req.material_stack["eps_r"] == 3.0
