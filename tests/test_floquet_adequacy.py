import warnings

import pytest

from src.basis_functions import MeshSettings, floquet_adequacy


def test_adequacy_flags_minimum_cutoff():
    ok, msg = floquet_adequacy(MeshSettings(60, 60, 31))
    assert not ok and "120" in msg
    assert floquet_adequacy(MeshSettings(60, 60, 120)) == (True, "")


def test_default_is_silent(monkeypatch):
    monkeypatch.delenv("LWMOM_FLOQUET_CHECK", raising=False)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        MeshSettings(60, 60, 31)


def test_opt_in_warns_and_errors(monkeypatch):
    monkeypatch.setenv("LWMOM_FLOQUET_CHECK", "warn")
    with pytest.warns(RuntimeWarning, match="floquet_order=31"):
        MeshSettings(60, 60, 31)
    monkeypatch.setenv("LWMOM_FLOQUET_CHECK", "error")
    with pytest.raises(ValueError, match="floquet_order=31"):
        MeshSettings(60, 60, 31)
    MeshSettings(60, 60, 120)


def test_supercell_default_mesh_is_unchanged_and_opt_in_factor():
    from src.ris.supercell import default_mesh

    m = default_mesh(48, 4, 1)
    assert (m.cells_x, m.cells_y, m.floquet_order) == (192, 48, 97)
    m2 = default_mesh(48, 4, 1, floquet_factor=2)
    assert (m2.cells_x, m2.cells_y, m2.floquet_order) == (192, 48, 384)
