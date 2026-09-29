"""Sec. 26-27: schema files are generated from the actual dataclasses, not hand-written."""

import json

from src.phase8d_rev.record_schema import UnitCellRecord
from src.phase8d_rev.schema_export import (
    build_database_schema_dict,
    build_library_build_request_schema_dict,
    render_database_schema_md,
    render_query_spec_md,
    write_database_schema_json,
    write_library_build_request_schema_json,
)


def test_database_schema_dict_lists_every_unit_cell_record_field():
    doc = build_database_schema_dict()
    import dataclasses

    field_names = {f.name for f in dataclasses.fields(UnitCellRecord)}
    assert set(doc["record_schema"]) == field_names


def test_database_schema_dict_lists_all_five_validation_levels():
    doc = build_database_schema_dict()
    assert len(doc["validation_levels"]) == 5


def test_write_database_schema_json_round_trips(tmp_path):
    path = write_database_schema_json(tmp_path / "database_schema.json")
    doc = json.loads(path.read_text())
    assert doc["format"] == "lightweight_mom.universal_unit_cell_database"


def test_render_database_schema_md_mentions_period_is_per_entry():
    text = render_database_schema_md()
    assert "per-LibraryEntry" in text or "per-entry" in text or "never global constants" in text


def test_render_query_spec_md_lists_the_filter_order():
    text = render_query_spec_md()
    for stage in ("topology", "material_stack", "frequency", "period", "polarization", "incidence"):
        assert stage in text


def test_library_build_request_schema_dict_notes_never_auto_submitted():
    doc = build_library_build_request_schema_dict()
    assert "never" in doc["note"].lower()


def test_write_library_build_request_schema_json_round_trips(tmp_path):
    path = write_library_build_request_schema_json(tmp_path / "lbr.json")
    doc = json.loads(path.read_text())
    assert doc["format"] == "lightweight_mom.library_build_request"
