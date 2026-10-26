# References

Reference papers are kept here locally for methodology lookup but **never committed**
(`docs/references/*.pdf` is gitignored — copyrighted third-party journal articles). Cite by DOI;
if a file listed below is missing from your checkout, ask whoever set up the project for a copy.

| File (gitignored) | Citation | Used for |
|---|---|---|
| `Mishra_TVT2026_multiport_1bit_OCMA.pdf` | P. K. Mishra, D. Barad, D. Sarkar, "A Multi-Port Multi-Beam Binary Coded Millimeter-Wave Metasurface Antenna for Full-Duplex Vehicular Communications," *IEEE Transactions on Vehicular Technology*, 2026. DOI: [10.1109/TVT.2026.3733618](https://doi.org/10.1109/TVT.2026.3733618). | Phase 15 Part C's multi-port sub-aperture layout and Phase 15 Part E's 1-bit coded design method (`docs/history/PHASE15_POLARIZATION_VORTEX.md`); reproduced at 26 GHz in `docs/submission/paper26/comparison.md`. |
| `Ghosh_OJAP2026_multipolarized_HMS.pdf` | Ghosh, Ghosh, Saha, *IEEE Open Journal of Antennas and Propagation*, 2026 (full title/DOI not independently confirmed by this project — the file itself is the source of record). | Phase 15 Part A/B's per-beam polarization compensation and vector far-field model (Table 1 and eqs. 1–17 of the paper); see `docs/history/PHASE15_POLARIZATION_VORTEX.md`. |

Also cited, not held as a local PDF:

- G. Minatti et al., "Metasurface Antennas: New Models, Applications and Realizations," surveyed
  in multiple IEEE TAP papers (2011–2016) — general polarization-control framing for Phase 15 Part A.
- Fong, Huang, Itoh et al. (2010) — the classic monopole-fed, grounded holographic-impedance-
  surface technique; found structurally incompatible with this project's ungrounded-slab solver
  core in Phase 8B (`docs/history/PHASE8B_FORMULATION.md`), motivating the reflectarray-style
  (reflection-phase) formulation this project actually uses for Formulation B.

Phase 20 R3 (2-bit PIN-diode RIS cells), cited by DOI only (none held locally): `results/phase20_r3/LITERATURE.md` lists them with what was verified
(10.3390/electronics12204220; 10.1038/srep42302; 10.1109/LAWP.2022.3188491; 10.1109/TAP.2024.3511945; 10.1364/OE.415730; 10.1109/TAP.2023.3291773;
10.1109/LAWP.2025.3606542; 10.1109/LAWP.2025.3646699; 10.36227/techrxiv.175019406.64903296; 10.1109/APUSNCURSINRSM.2017.8073084; 10.1155/2022/6011557).

If you add a new reference PDF, add its row here (citation + DOI if available + what it's used
for) in the same commit that cites it from a phase doc or module docstring.
