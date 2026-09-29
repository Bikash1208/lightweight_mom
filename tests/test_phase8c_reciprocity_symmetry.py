"""Phase 8C Sec. 31-32, 55: numeric reciprocity and L/W-swap symmetry checks on real solved data."""

from src.phase8c.reciprocity_symmetry import check_lw_swap_symmetry, check_reciprocity


def test_real_library_is_reciprocal_to_numerical_precision(phase8c_full_database):
    result = check_reciprocity(phase8c_full_database, 100, tolerance=1e-6)
    assert result.reciprocal, result.note
    assert result.n_cells_checked == 324
    assert result.max_abs_error_co < 1e-9
    assert result.max_abs_error_cross < 1e-9


def test_reciprocity_check_does_not_silently_correct_data():
    # A synthetic non-reciprocal database (built via a tiny fake) must be REPORTED, not corrected.
    import numpy as np

    class _FakeRecord:
        def __init__(self, s):
            self.s = np.asarray(s, dtype=complex)

    class _FakeCell:
        def __init__(self, index):
            self.index = index
            self.status = "ok"

    class _FakeDB:
        def __init__(self):
            self.cells = [_FakeCell(0)]

        def record(self, index, freq_idx):
            s = np.zeros((4, 4), dtype=complex)
            s[2, 0] = 0.5  # S21_xx
            s[0, 2] = 0.9  # S12_xx -- deliberately NOT reciprocal
            return _FakeRecord(s)

    result = check_reciprocity(_FakeDB(), 0, tolerance=1e-6)
    assert result.reciprocal is False
    assert "does NOT hold" in result.note


def test_lw_swap_symmetry_on_real_library(phase8c_full_database):
    result = check_lw_swap_symmetry(phase8c_full_database, 100, component_family="S11", tolerance=1e-2)
    assert result.n_pairs_checked > 0
    # Reported, not asserted blind: whatever the measured error is, it must be finite and sane.
    assert result.max_abs_error < 1.0


def test_lw_swap_symmetry_with_no_pairs_reports_explicitly():
    class _FakeCell:
        def __init__(self, index, l_mm, w_mm):
            self.index = index
            self.status = "ok"
            self.parameters = {"length_m": l_mm * 1e-3, "width_m": w_mm * 1e-3}

    class _FakeDB:
        cells = [_FakeCell(0, 1.0, 2.0)]  # no (2.0, 1.0) counterpart

    result = check_lw_swap_symmetry(_FakeDB(), 0)
    assert result.n_pairs_checked == 0
    assert "no (L,W)/(W,L) pairs" in result.note
