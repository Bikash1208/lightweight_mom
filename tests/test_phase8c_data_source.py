"""Phase 8C data source: the full linear-basis S-matrix, cross-checked against the delivered
library."""

from src.phase8c import data_source


def test_load_full_database_has_324_cells_and_full_linear_basis(phase8c_full_database):
    db = phase8c_full_database
    assert len(db.cells) == 324
    ok = [c for c in db.cells if c.status == "ok"]
    assert len(ok) == 324
    record = db.record(ok[0].index, 0)
    from src.database import family_block

    block = family_block(record.s, "S11")
    assert block.shape == (2, 2)   # xx, xy, yx, yy all present


def test_cross_check_against_delivered_library_passes(phase8c_full_database):
    result = data_source.cross_check_against_delivered_library(phase8c_full_database)
    assert result["ok"], result["mismatches"]
    assert result["checked"] > 0


def test_load_full_database_missing_checkpoint_raises(tmp_path):
    import pytest

    with pytest.raises(FileNotFoundError, match="internal Phase 7 solve checkpoint"):
        data_source.load_full_database(tmp_path / "nonexistent.json")
