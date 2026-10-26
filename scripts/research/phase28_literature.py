"""Phase 28 item 1: the literature table, built ONLY from Crossref metadata (title, authors, venue, year, abstract where Crossref has one).

    python scripts/research/phase28_literature.py fetch    # network: writes results/phase28/crossref/<doi>.json (trimmed Crossref records)
    python scripts/research/phase28_literature.py build    # offline: results/phase28/crossref/*.json -> results/phase28/LITERATURE.md

Content columns the brief asks for (pitch, diodes per cell, bit depth, switching contrast, aperture, gain, scan range, bias approach, measured vs simulated) are filled only
with a QUOTE from the stored title/abstract that supports them; otherwise the cell says NEEDS_FULL_TEXT. `CLAIMS` quotes are checked against the stored record by
tests/test_phase28_literature.py. No numbers are computed from the table.
"""

import json
import re
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "phase28"
CR = OUT / "crossref"
RETRIEVED = "2026-10-26"
NEEDS = "NEEDS_FULL_TEXT (user has IEEE access)"

# DOI -> list of (what the quote supports, verbatim quote from the Crossref title or abstract).
CLAIMS = {
    "10.1109/TAP.2003.817558": [
        (
            "tunable surface, 2-D steering",
            "Two-dimensional beam steering using an electrically tunable impedance surface",
        )
    ],
    "10.1109/TAP.2010.2055812": [
        (
            "scalar and tensor (polarisation) holographic surfaces",
            "Scalar and Tensor Holographic Artificial Impedance Surfaces",
        )
    ],
    "10.1109/TAP.2011.2143668": [
        ("sinusoidally modulated reactance, printed leaky wave", "Sinusoidally-Modulated Reactance Surface")
    ],
    "10.1109/TAP.2014.2377718": [
        ("modulated metasurface antennas (passive)", "Modulated Metasurface Antennas for Space")
    ],
    "10.1109/TAP.2015.2399937": [
        ("reconfigurable holographic antenna", "Reconfigurable Holographic Metamaterial Antenna")
    ],
    "10.1109/ACCESS.2017.2712659": [
        ("reconfigurable holographic aperture", "Reconfigurable Holographic Metasurface Aperture")
    ],
    "10.1109/TAP.2018.2845746": [
        ("binary (1-bit) programmable surface, leaky wave, fixed frequency", "Binary Programmable Metasurface")
    ],
    "10.1038/s41598-021-83377-9": [
        ("varactor-tuned elements", "varactor-tuned elements"),
        ("half-wavelength spatial sampling", "spatial sampling at the Nyquist limit of half a wavelength"),
        (
            "beamsteering demonstrated experimentally in two directions",
            "experimentally demonstrate electronic beamsteering in two directions",
        ),
    ],
    "10.1109/LAWP.2023.3239992": [
        (
            "circular-polarisation reconfigurable holographic leaky wave",
            "Circular-Polarization-Reconfigurable Holographic Leaky-Wave Antenna",
        )
    ],
}
# Why each paper is in the table (role for Phase 28), from its title only.
ROLE = {
    "10.1109/TAP.2003.817558": "tunable-impedance beam steering (early)",
    "10.1109/TAP.2010.2055812": "holographic artificial impedance surface, tensor = polarisation control (passive)",
    "10.1109/TAP.2011.2143668": "sinusoidally modulated reactance surface leaky wave (passive)",
    "10.1109/TAP.2014.2377718": "modulated metasurface antennas, synthesis and realisation (passive)",
    "10.1109/TAP.2015.2399937": "reconfigurable holographic antenna, sidelobe control",
    "10.1109/ACCESS.2017.2712659": "reconfigurable holographic aperture",
    "10.1109/TAP.2018.2845746": "binary programmable leaky-wave metasurface",
    "10.1038/s41598-021-83377-9": "electronically steered waveguide-fed metasurface antenna",
    "10.1109/LAWP.2023.3239992": "CP-reconfigurable holographic leaky-wave antenna",
}
DOIS = list(CLAIMS)


def _clean(t):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", t or "")).strip()


def fetch():
    CR.mkdir(parents=True, exist_ok=True)
    for doi in DOIS:
        req = urllib.request.Request(
            f"https://api.crossref.org/works/{doi}",
            headers={"User-Agent": "lightweight_mom/phase28 (mailto:bikashkumarshanti@gmail.com)"},
        )
        m = json.load(urllib.request.urlopen(req, timeout=30))["message"]
        rec = {
            "doi_requested": doi,
            "doi_crossref": m["DOI"],
            "title": _clean(m["title"][0]),
            "authors": [f"{a.get('given', '')} {a.get('family', '')}".strip() for a in m.get("author", [])],
            "venue": (m.get("container-title") or [""])[0],
            "year": m["issued"]["date-parts"][0][0],
            "volume": m.get("volume"),
            "abstract": _clean(m.get("abstract")),
            "retrieved": RETRIEVED,
        }
        text = json.dumps(rec, indent=1, ensure_ascii=False) + "\n"
        (CR / (doi.replace("/", "_") + ".json")).write_text(text, encoding="utf-8", newline="\n")
        print(doi, "->", rec["title"][:70])


def load():
    return {doi: json.loads((CR / (doi.replace("/", "_") + ".json")).read_text(encoding="utf-8")) for doi in DOIS}


def build():
    recs = load()
    n_abs = sum(1 for r in recs.values() if r["abstract"])
    lines = [
        "# Phase 28 item 1 - literature (verified DOIs only)",
        "",
        f"Every row is built from the Crossref record of its DOI (retrieved {RETRIEVED}, stored in `results/phase28/crossref/`; title, authors, venue and year are Crossref's, not typed by hand).",
        f"Crossref holds an abstract for only {n_abs} of the {len(recs)} papers. The brief's content fields (cell pitch vs lambda, diodes per cell, bit depth, switching contrast, aperture, gain,",
        f"scan range, polarisation, bias approach, measured vs simulated) are therefore given only where the title or abstract states them (quoted); every other field is `{NEEDS}`.",
        "No PDFs are committed; no statistics are computed over the table. Status label: `LITERATURE_METADATA_VERIFIED`; content columns `NEEDS_FULL_TEXT`.",
        "",
        "Dropped from the first attempt: Chen et al. 2020 (`10.1109/TAP.2019.2936306`): Crossref returns 'Resource not found' and no real DOI for it was identified. Minatti 2015 and Patel & Grbic 2011 are included with the corrected DOIs below.",
        "",
        "| # | Authors | Title | Venue, year | DOI | Role in Phase 28 (from the title) | Stated in title / abstract (quoted) |",
        "|---|---|---|---|---|---|---|",
    ]
    for i, doi in enumerate(DOIS, 1):
        r = recs[doi]
        au = r["authors"]
        authors = au[0].split()[-1] + " et al." if len(au) > 2 else " & ".join(a.split()[-1] for a in au)
        quotes = "; ".join(f'{what}: "{q}"' for what, q in CLAIMS[doi])
        lines.append(
            f"| {i} | {authors} | {r['title']} | {r['venue']}, {r['year']} | [{doi}](https://doi.org/{doi}) | {ROLE[doi]} | {quotes} |"
        )
    lines += [
        "",
        "## Fields not available without the full text",
        "",
        f"For all {len(recs)} papers: cell pitch vs lambda, diodes per cell, bit depth, switching contrast (surface-reactance or beta/alpha change), aperture size, gain, scan range, polarisation (beyond the title), bias approach, measured vs simulated: `{NEEDS}`.",
        "The abstract of Boyarsky 2021 (the only abstract in Crossref) says: varactor-tuned elements, half-wavelength sampling, beam steering demonstrated experimentally in two directions; nothing on diodes, bias or numbers.",
        "",
        "## What this table supports (and what it does not)",
        "",
        "- The papers split into passive holographic / modulated-reactance antennas (Fong 2010, Patel & Grbic 2011, Minatti 2015) and electronically reconfigurable ones (Sievenpiper 2003, Johnson 2015, Yurduseven 2017, Wan 2018, Boyarsky 2021, Wang 2023), judged from the titles.",
        "- Tensor surfaces for polarisation control (Fong 2010) and CP-reconfigurable holographic leaky-wave antennas (Wang 2023) exist in the literature, judged from the titles.",
        "- It does NOT support any statement on a required switching contrast, a recommended pitch, bit-depth loss, diode counts or bias architecture: those need the full texts. Item 2 derives the contrast requirement from the Phase 26 model instead of from this table.",
    ]
    (OUT / "LITERATURE.md").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print("wrote", OUT / "LITERATURE.md")


if __name__ == "__main__":
    {"fetch": fetch, "build": build}[sys.argv[1]]()
