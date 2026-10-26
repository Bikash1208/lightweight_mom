# Phase 28 item 1 - literature (verified DOIs only)

Every row is built from the Crossref record of its DOI (retrieved 2026-10-26, stored in `results/phase28/crossref/`; title, authors, venue and year are Crossref's, not typed by hand).
Crossref holds an abstract for only 1 of the 9 papers. The brief's content fields (cell pitch vs lambda, diodes per cell, bit depth, switching contrast, aperture, gain,
scan range, polarisation, bias approach, measured vs simulated) are therefore given only where the title or abstract states them (quoted); every other field is `NEEDS_FULL_TEXT (user has IEEE access)`.
No PDFs are committed; no statistics are computed over the table. Status label: `LITERATURE_METADATA_VERIFIED`; content columns `NEEDS_FULL_TEXT`.

Dropped from the first attempt: Chen et al. 2020 (`10.1109/TAP.2019.2936306`): Crossref returns 'Resource not found' and no real DOI for it was identified. Minatti 2015 and Patel & Grbic 2011 are included with the corrected DOIs below.

| # | Authors | Title | Venue, year | DOI | Role in Phase 28 (from the title) | Stated in title / abstract (quoted) |
|---|---|---|---|---|---|---|
| 1 | Sievenpiper et al. | Two-dimensional beam steering using an electrically tunable impedance surface | IEEE Transactions on Antennas and Propagation, 2003 | [10.1109/TAP.2003.817558](https://doi.org/10.1109/TAP.2003.817558) | tunable-impedance beam steering (early) | tunable surface, 2-D steering: "Two-dimensional beam steering using an electrically tunable impedance surface" |
| 2 | Fong et al. | Scalar and Tensor Holographic Artificial Impedance Surfaces | IEEE Transactions on Antennas and Propagation, 2010 | [10.1109/TAP.2010.2055812](https://doi.org/10.1109/TAP.2010.2055812) | holographic artificial impedance surface, tensor = polarisation control (passive) | scalar and tensor (polarisation) holographic surfaces: "Scalar and Tensor Holographic Artificial Impedance Surfaces" |
| 3 | Patel & Grbic | A Printed Leaky-Wave Antenna Based on a Sinusoidally-Modulated Reactance Surface | IEEE Transactions on Antennas and Propagation, 2011 | [10.1109/TAP.2011.2143668](https://doi.org/10.1109/TAP.2011.2143668) | sinusoidally modulated reactance surface leaky wave (passive) | sinusoidally modulated reactance, printed leaky wave: "Sinusoidally-Modulated Reactance Surface" |
| 4 | Minatti et al. | Modulated Metasurface Antennas for Space: Synthesis, Analysis and Realizations | IEEE Transactions on Antennas and Propagation, 2015 | [10.1109/TAP.2014.2377718](https://doi.org/10.1109/TAP.2014.2377718) | modulated metasurface antennas, synthesis and realisation (passive) | modulated metasurface antennas (passive): "Modulated Metasurface Antennas for Space" |
| 5 | Johnson et al. | Sidelobe Canceling for Reconfigurable Holographic Metamaterial Antenna | IEEE Transactions on Antennas and Propagation, 2015 | [10.1109/TAP.2015.2399937](https://doi.org/10.1109/TAP.2015.2399937) | reconfigurable holographic antenna, sidelobe control | reconfigurable holographic antenna: "Reconfigurable Holographic Metamaterial Antenna" |
| 6 | Yurduseven et al. | Design and Analysis of a Reconfigurable Holographic Metasurface Aperture for Dynamic Focusing in the Fresnel Zone | IEEE Access, 2017 | [10.1109/ACCESS.2017.2712659](https://doi.org/10.1109/ACCESS.2017.2712659) | reconfigurable holographic aperture | reconfigurable holographic aperture: "Reconfigurable Holographic Metasurface Aperture" |
| 7 | Wan et al. | Beam Forming of Leaky Waves at Fixed Frequency Using Binary Programmable Metasurface | IEEE Transactions on Antennas and Propagation, 2018 | [10.1109/TAP.2018.2845746](https://doi.org/10.1109/TAP.2018.2845746) | binary programmable leaky-wave metasurface | binary (1-bit) programmable surface, leaky wave, fixed frequency: "Binary Programmable Metasurface" |
| 8 | Boyarsky et al. | Electronically steered metasurface antenna | Scientific Reports, 2021 | [10.1038/s41598-021-83377-9](https://doi.org/10.1038/s41598-021-83377-9) | electronically steered waveguide-fed metasurface antenna | varactor-tuned elements: "varactor-tuned elements"; half-wavelength spatial sampling: "spatial sampling at the Nyquist limit of half a wavelength"; beamsteering demonstrated experimentally in two directions: "experimentally demonstrate electronic beamsteering in two directions" |
| 9 | Wang et al. | A Circular-Polarization-Reconfigurable Holographic Leaky-Wave Antenna With Fixed-Frequency Beam-Scanning Property | IEEE Antennas and Wireless Propagation Letters, 2023 | [10.1109/LAWP.2023.3239992](https://doi.org/10.1109/LAWP.2023.3239992) | CP-reconfigurable holographic leaky-wave antenna | circular-polarisation reconfigurable holographic leaky wave: "Circular-Polarization-Reconfigurable Holographic Leaky-Wave Antenna" |

## Fields not available without the full text

For all 9 papers: cell pitch vs lambda, diodes per cell, bit depth, switching contrast (surface-reactance or beta/alpha change), aperture size, gain, scan range, polarisation (beyond the title), bias approach, measured vs simulated: `NEEDS_FULL_TEXT (user has IEEE access)`.
The abstract of Boyarsky 2021 (the only abstract in Crossref) says: varactor-tuned elements, half-wavelength sampling, beam steering demonstrated experimentally in two directions; nothing on diodes, bias or numbers.

## What this table supports (and what it does not)

- The papers split into passive holographic / modulated-reactance antennas (Fong 2010, Patel & Grbic 2011, Minatti 2015) and electronically reconfigurable ones (Sievenpiper 2003, Johnson 2015, Yurduseven 2017, Wan 2018, Boyarsky 2021, Wang 2023), judged from the titles.
- Tensor surfaces for polarisation control (Fong 2010) and CP-reconfigurable holographic leaky-wave antennas (Wang 2023) exist in the literature, judged from the titles.
- It does NOT support any statement on a required switching contrast, a recommended pitch, bit-depth loss, diode counts or bias architecture: those need the full texts. Item 2 derives the contrast requirement from the Phase 26 model instead of from this table.
