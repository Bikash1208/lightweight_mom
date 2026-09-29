"""Sec. 16-17: period-independent analysis, reusing phase8d_rev."""

import pytest

from src.phase8c import data_source
from src.phase8d_rev import library_entry as le
from src.phase8d_rev.universal_database import UniversalDatabase
from src.phase8e.candidate_library_period import assess_period_requirement

pytestmark = pytest.mark.skipif(not data_source.DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present")


@pytest.fixture
def udb(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    u = UniversalDatabase()
    u.register(entry)
    return u


def test_registered_topology_at_a_reasonable_k_required_is_current_dataset_limitation(udb):
    # k_required_max=200 rad/m brackets the registered 4.8 mm period within the Nyquist
    # candidate-period range (order-of-magnitude check); see candidate_library_period.py.
    result = assess_period_requirement(udb, "rectangular_patch", 200.0, 12.5e9, "Rogers RT/Duroid 5880")
    assert result.status == "CURRENT_DATASET_LIMITATION"
    assert result.required_library_specification is None


def test_unregistered_topology_produces_a_required_library_specification(udb):
    result = assess_period_requirement(udb, "slot_loaded_patch", 50.0, 12.5e9)
    assert result.status == "REQUIRED_LIBRARY_MISSING"
    assert result.required_library_specification is not None
    assert result.required_library_specification["topology"] == "slot_loaded_patch"


def test_candidate_periods_are_never_labelled_validated(udb):
    result = assess_period_requirement(udb, "rectangular_patch", 50.0, 12.5e9)
    # the candidate list itself only carries mm values here (labels live in candidate_period.py,
    # tested there); this test only checks the assessment does not silently call them VALIDATED.
    assert result.status != "VALIDATED"
