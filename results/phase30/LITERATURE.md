# Phase 30 item 1 - literature table (journals only)

All DOIs verified against Crossref on 2026-11-02 (title/journal/year match). Cells marked
`NOT_VERIFIED` could not be read from an accessible full text in this session (IEEE Xplore pages are
JS-only here; the Ulm GeMiC conference PDF is excluded by the journals-only rule). No PDFs committed.

| Ref | DOI (verified) | Journal |
|---|---|---|
| Faenzi, Minatti, Gonzalez-Ovejero, Caminita, Martini, Della Giovampaola, Maci, "Metasurface Antennas: New Models, Applications and Realizations" | 10.1038/s41598-019-46522-z | Sci. Rep. 9, 2019 (open access, PMC6629853) |
| Minatti, Faenzi, Martini, Caminita, De Vita, Gonzalez-Ovejero, Maci, "Modulated Metasurface Antennas for Space: Synthesis, Analysis and Realizations" | 10.1109/TAP.2014.2377718 | IEEE TAP 63(4), 1288-1300, 2015 |
| Fong et al., "Scalar and Tensor Holographic Artificial Impedance Surfaces" | 10.1109/TAP.2010.2055812 | IEEE TAP, 2010 |
| Teniou et al., "Implementation of Radiating Aperture Field Distribution Using Tensorial Metasurfaces" | 10.1109/TAP.2017.2751652 | IEEE TAP, 2017 |
| Frey, Doering, Riese, Waldschmidt, Chaloun, "Angle-Dependent Synthesis Method for Holographic Multi-Feed Antennas" | 10.1109/OJAP.2023.3261903 | IEEE OJAP 4, 392-407, 2023 |

Note: the brief lists Frey OJAP 2023 as "elliptical patch with centre slot"; Crossref confirms only the title, the cell geometry in that paper is UNCONFIRMED (see status section below).

## Table (Faenzi 2019, read from open full text; other rows not read)

| Design (Faenzi 2019) | f (GHz) | Substrate | Cell | Aperture | Measured |
|---|---|---|---|---|---|
| High-efficiency broadside | 29.75 | RO3003 0.762 mm, er 3, tand 0.001 | grain of rice | radius 13.5 lambda | 37 dBi, ~70 % overall (feed x ohmic ~80 %) |
| 30 deg tilted, RHCP | 20 | not stated in text read | grain of rice | radius 10 lambda | 33 dBi, ~58 % aperture eff., AR < 1 dB in HPBW, XP ~ -30 dB |
| Multibeam (2 superimposed modulations, tapered) | 17 | RO3010 0.635 mm, er 10.2 | elliptical patch | radius 6 lambda | 23.5 dBi per beam, 3 dB BW 5.9 % |
| Wideband (simulated only) | 20-30 / 24-30 | RO3006 0.635 mm, er 6.15 | continuous-impedance MoM model | 16.7 / 11.1 cm | ~29 dB avg gain; mean X -450 / -260 ohm |
| Isoflux | 26.7 | er 9.8, 0.5 mm | double anchor (prototype families) | radius 9 lambda | gain not read |

Other facts stated there: pixel families = coffee bean, patch with slot, grain of rice, patch with cross
slot, double pi, double anchor; cell size lambda/10 to lambda/5; each pixel parameterised by size and
rotation psi; impedance is a tensor (diagonal + off-diagonal) in cylindrical components; "high
anisotropy control" credited to grain of rice; amplitude (leakage) set by modulation index M, phase by
the modulation function; feed = coax monopole at centre with annular patch/slot matching. M range,
pitch in mm, and anisotropy ratio are not given in the part of the text read (`NOT_AVAILABLE_FROM_CURRENT_DATA`).
Edge absorber/rim: not mentioned in the part read; published boards are 6-13.5 lambda in radius vs
our 150 mm hexagon (~2.6 lambda radius), i.e. 2-5x larger, so they reach low edge power by length
rather than by a rim.

## Status of the four papers not read (2026-11-02, item 2 of the follow-up)
Bibliographic data verified against Crossref (title, journal, year, volume, pages; authors for Fong and Teniou):
- Fong, Colburn, Ottusch, Visher, Sievenpiper, IEEE TAP 58(10), 3212-3221, 2010, 10.1109/TAP.2010.2055812 - DOI VERIFIED, content NOT read.
- Minatti, Faenzi, Martini, Caminita, De Vita, Gonzalez-Ovejero, Maci, IEEE TAP 63(4), 1288-1300, 2015, 10.1109/TAP.2014.2377718 - DOI VERIFIED, content NOT read.
- Teniou, Roussel, Capet, Piau, Casaletti, IEEE TAP 65(11), 5895-5907, 2017, 10.1109/TAP.2017.2751652 - DOI VERIFIED, content NOT read.
- Frey, Doering, Riese, Waldschmidt, Chaloun, IEEE OJAP 4, 392-407, 2023, 10.1109/OJAP.2023.3261903 - DOI VERIFIED, content NOT read.
Full text could not be retrieved by the tools available (IEEE Xplore needs a browser session; an automated download was
refused by the session's safety layer and is not repeated). **No claim in this repository is based on these four papers' contents**:
every table row above is Faenzi 2019 (open access, read in full). The cell-geometry attribution of the brief ("Frey: elliptical
patch with centre slot") is therefore UNCONFIRMED. If the PDFs/text are provided, the Fong/Minatti/Teniou/Frey rows
(cell, substrate, pitch, reactance range, anisotropy, efficiency, edge treatment) can be filled in.

## Conclusions (from verified rows + physics; labelled as inference where so)
1. Every published efficient design uses substrates thinner/higher-er than ours (er 3-10, 0.5-0.76 mm vs
   RT5880 er 2.2, 1.575 mm); high-er thin slabs give a more strongly bound TM wave and a wider X range
   per cell size (inference; to be quantified in item 2).
2. Cells are lambda/10-lambda/5 pitch with size + rotation as the two tensor knobs; our 6.4 mm
   (~lambda/4.4 free-space at 10.54 GHz) is at the coarse end of that range.
3. Grain of rice / elliptical patch are the demonstrated CP/tensor workhorses; cross-slot and slotted
   patches are listed families but their quantitative swing was not read.
4. Aperture electrical size is the main lever the literature uses to get efficient leakage: with radius
   6-13.5 lambda, a modest alpha suffices. At ~2.6 lambda radius our board needs alpha/k several times higher
   (consistent with the brief's alpha/k >= 0.03 target) - confirm with item 2.
5. Multibeam (Faenzi) superimposes modulations with a taper; this is the reference for the 7-port aim.
