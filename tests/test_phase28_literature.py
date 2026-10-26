"""Phase 28 item 1: the literature table is built from stored Crossref records; every quoted claim occurs in the record; DOIs resolve to the expected paper."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase28_literature as lit  # noqa: E402


def test_every_quote_occurs_in_the_stored_crossref_title_or_abstract():
    for doi, rec in lit.load().items():
        text = rec["title"] + " " + rec["abstract"]
        for what, quote in lit.CLAIMS[doi]:
            assert quote in text, (doi, what, quote)


def test_stored_records_are_for_the_requested_doi_and_the_corrected_dois_are_present():
    recs = lit.load()
    for doi, rec in recs.items():
        assert rec["doi_crossref"].lower() == doi.lower()
    assert recs["10.1109/TAP.2014.2377718"]["authors"][0].endswith("Minatti")
    assert "Sinusoidally-Modulated" in recs["10.1109/TAP.2011.2143668"]["title"]
    assert "10.1109/TAP.2019.2936306" not in recs  # Chen 2020: not found in Crossref, dropped


def test_literature_md_is_in_sync_with_the_records_and_every_doi_appears():
    md = (lit.OUT / "LITERATURE.md").read_text(encoding="utf-8")
    for doi, rec in lit.load().items():
        assert f"[{doi}]" in md and rec["title"] in md
    assert lit.NEEDS in md
    assert "meta-analysis" not in md.lower()


def test_build_is_deterministic(tmp_path, monkeypatch):
    before = (lit.OUT / "LITERATURE.md").read_bytes()
    monkeypatch.setattr(lit, "OUT", tmp_path)
    lit.build()
    assert (tmp_path / "LITERATURE.md").read_bytes() == before
