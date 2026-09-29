"""Shared, session-scoped fixtures for tests that need the full Phase 7 checkpoint database
(`results/unit_cell_library/_solve/database.json`, ~470 MB, ~20 s to parse). Phase 8C is the
first phase whose tests need the *full* linear-basis S-matrix (Phase 7/8A/8B tests only ever
built small, fast, synthetic-scale databases in memory) -- loading it once per test session,
rather than once per test file, keeps the added cost to a single load instead of N.
"""

import pytest

from src.phase8c import data_source


@pytest.fixture(scope="session")
def phase8c_full_database():
    if not data_source.DEFAULT_CHECKPOINT_JSON.exists():
        pytest.skip("Phase 7 production library checkpoint not present (results/unit_cell_library/_solve/)")
    return data_source.load_full_database()
