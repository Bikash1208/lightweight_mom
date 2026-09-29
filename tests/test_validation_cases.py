"""Phase 6 validation-case definitions: every case builds a valid config, the requested
category checklist is covered, and construction is deterministic (`config.config_from_dict`
and `sweep.FrequencyGrid` are Phase 1/5, unmodified -- these tests are about the case *catalog*,
not about re-testing that machinery)."""

import pytest

from src.config import SimulationConfig
from src.sweep import FrequencyGrid
from src.validation_cases import CASES, CASES_BY_ID, CATEGORIES, ValidationCase, _case

REQUIRED_CATEGORIES = {
    "non-resonant", "resonant", "narrow", "wide", "different lengths/widths", "rotated", "slot",
}


def test_ten_cases_with_unique_ids():
    assert len(CASES) == 10
    ids = [c.case_id for c in CASES]
    assert len(set(ids)) == 10
    assert CASES_BY_ID == {c.case_id: c for c in CASES}


def test_every_case_is_well_typed():
    for c in CASES:
        assert isinstance(c, ValidationCase)
        assert isinstance(c.config, SimulationConfig)
        assert isinstance(c.frequency, FrequencyGrid)
        assert c.case_id and isinstance(c.case_id, str)
        assert c.category and isinstance(c.category, str)
        assert c.description and isinstance(c.description, str)


def test_all_required_categories_present():
    """Every category the Phase 6 brief asks for (non-resonant, resonant, narrow, wide, several
    length/width pairs, rotated, slot/complement) has at least one case."""
    present = set(CATEGORIES)
    assert REQUIRED_CATEGORIES <= present, f"missing categories: {REQUIRED_CATEGORIES - present}"


def test_rotated_cases_have_nonzero_rotation():
    rotated = [c for c in CASES if c.category == "rotated"]
    assert len(rotated) >= 2
    for c in rotated:
        assert c.config.unit_cell.geometry.rotation != 0.0


def test_slot_cases_are_aperture_polarity():
    slots = [c for c in CASES if c.category == "slot"]
    assert len(slots) >= 2
    for c in slots:
        assert c.config.unit_cell.geometry.polarity == "aperture"


def test_narrow_and_wide_have_opposite_aspect_ratios():
    narrow = CASES_BY_ID["case03_narrow_dipole"].config.unit_cell.geometry
    wide = CASES_BY_ID["case04_wide_patch"].config.unit_cell.geometry
    assert narrow.length > narrow.width
    assert wide.width > wide.length


def test_different_length_width_pairs_are_actually_different():
    """The module docstring's claim: cases 1, 4, 5, 6 span four unrelated (length, width) pairs
    (other cases may coincidentally share a pair with one of these while testing an unrelated
    property, e.g. case09's slot happens to reuse case08's (5, 2) mm -- that overlap is fine, it
    is these four "different lengths/widths" cases that must be pairwise distinct)."""
    ids = ("case01_non_resonant_patch", "case04_wide_patch", "case05_small_square_patch", "case06_large_rect_patch")
    pairs = {(CASES_BY_ID[i].config.unit_cell.geometry.length, CASES_BY_ID[i].config.unit_cell.geometry.width) for i in ids}
    assert len(pairs) == len(ids)


def test_offset_case_has_nonzero_center():
    offset = CASES_BY_ID["case08_rotated_offset_patch"].config.unit_cell.geometry
    assert (offset.center_x, offset.center_y) != (0.0, 0.0)


def test_all_cases_share_lattice_and_substrate():
    """The module docstring's claim: one lattice/substrate/reference-plane config shared by every
    case, so only geometry/mesh/frequency vary between them."""
    first = CASES[0].config
    for c in CASES:
        uc = c.config.unit_cell
        assert uc.lattice.period_x == first.unit_cell.lattice.period_x
        assert uc.lattice.period_y == first.unit_cell.lattice.period_y
        assert uc.substrate == first.unit_cell.substrate
        assert c.config.reference_planes == first.reference_planes


def test_case_construction_is_deterministic():
    """Calling `_case` twice with identical arguments yields configs with identical parameters
    (no hidden randomness or mutable shared state between cases)."""
    a = _case("dup", "test", "duplicate", {"type": "rectangle", "length": 3.0e-3, "width": 1.0e-3}, 8.0e9, 12.0e9)
    b = _case("dup", "test", "duplicate", {"type": "rectangle", "length": 3.0e-3, "width": 1.0e-3}, 8.0e9, 12.0e9)
    assert a.config.unit_cell == b.config.unit_cell
    assert a.frequency == b.frequency


def test_frequency_grids_are_within_a_sane_microwave_band():
    for c in CASES:
        assert 1.0e9 <= c.frequency.start < c.frequency.stop <= 40.0e9


@pytest.mark.parametrize("case", CASES, ids=[c.case_id for c in CASES])
def test_each_case_builds_a_valid_unit_cell(case):
    """UnitCell.__post_init__ (Phase 1) already enforces "feature inside the cell" and "nonzero
    metal area"; constructing every case without error is itself a check that the catalog is
    internally consistent (this reuses, not re-implements, that validation)."""
    uc = case.config.unit_cell
    assert uc.metal_area > 0.0
    assert uc.metal_fill_fraction < 1.0
