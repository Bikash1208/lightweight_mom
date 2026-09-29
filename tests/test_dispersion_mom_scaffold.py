"""Part 3: the MoM surface-wave-eigenmode scaffold -- confirms it fails loudly and specifically,
never silently, and never returns a fabricated dispersion point."""

import pytest

from src import errors
from src.dispersion import mom_scaffold


def test_find_surface_wave_mode_raises_the_documented_code():
    with pytest.raises(errors.LwmomError) as exc_info:
        mom_scaffold.find_surface_wave_mode(unit_cell=None, frequency_hz=10e9)
    assert exc_info.value.code == "DISPERSION_MOM_NOT_IMPLEMENTED"
    assert "hfss_table" in exc_info.value.remedy


def test_scaffold_labels_include_experimental():
    assert "EXPERIMENTAL" in mom_scaffold.SCAFFOLD_LABELS
    assert "DISPERSION_FROM_MOM" in mom_scaffold.SCAFFOLD_LABELS
