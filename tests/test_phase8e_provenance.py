"""Sec. 27: formulation provenance registry."""

from src.phase8e.provenance import build_provenance_registry


def test_registry_has_all_four_candidates():
    reg = build_provenance_registry()
    assert set(reg) == {"A", "B", "C", "D"}


def test_a_and_d_document_incompatibilities():
    reg = build_provenance_registry()
    assert reg["A"].incompatibilities
    assert reg["D"].incompatibilities


def test_b_and_c_share_the_identical_mathematical_kernel_reference():
    reg = build_provenance_registry()
    assert "phase_difference_hologram" in reg["B"].mathematical_reference
    assert "phase_difference_hologram" in reg["C"].mathematical_reference


def test_every_entry_has_a_non_empty_implementation_status():
    reg = build_provenance_registry()
    for letter, prov in reg.items():
        assert prov.implementation_status, letter
