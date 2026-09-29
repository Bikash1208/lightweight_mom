"""Phase 8B Sec. 33-36, 45-46: machine-readable output schema and source traceability."""

import json

from src.phase8b.specification import build_holographic_specification, build_source_registry


def test_specification_has_required_top_level_keys():
    doc = build_holographic_specification()
    for key in (
        "status",
        "time_convention",
        "coordinate_system",
        "reference_wave",
        "object_wave",
        "multi_beam",
        "holographic_modulation",
        "control_variable",
        "library_compatibility",
        "hfss_validation",
    ):
        assert key in doc, key


def test_specification_status_is_never_a_bare_pass():
    doc = build_holographic_specification()
    assert doc["status"] != "PASS"
    from src.phase8b.library_compatibility import FINAL_STATES

    assert doc["status"] in FINAL_STATES


def test_specification_hfss_validation_is_never_claimed_performed():
    doc = build_holographic_specification()
    assert doc["hfss_validation"] == {"performed": False}


def test_specification_time_convention_is_preserved():
    doc = build_holographic_specification()
    assert doc["time_convention"] == "exp(+j*omega*t)"


def test_specification_lists_both_candidate_formulations():
    doc = build_holographic_specification()
    names = [f["name"] for f in doc["holographic_modulation"]["formulations"]]
    assert any(n.startswith("A:") for n in names)
    assert any(n.startswith("B:") for n in names)


def test_specification_is_json_serializable():
    doc = build_holographic_specification()
    json.dumps(doc)  # must not raise


def test_source_registry_has_citation_for_every_claim():
    registry = build_source_registry()
    assert len(registry) >= 2
    for key, entry in registry.items():
        for field in (
            "claim",
            "source",
            "equation_section",
            "project_convention_compatibility",
            "implementation_consequence",
            "validation_status",
        ):
            assert field in entry, f"{key} missing {field}"
            assert entry[field]  # non-empty


def test_source_registry_fong2010_citation_is_present():
    registry = build_source_registry()
    assert "Fong" in registry["fong2010"]["source"]
    assert "2010" in registry["fong2010"]["source"]


def test_source_registry_is_json_serializable():
    json.dumps(build_source_registry())
