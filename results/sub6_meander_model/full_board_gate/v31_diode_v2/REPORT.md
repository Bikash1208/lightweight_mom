# Corrected B-ALL-DIODE-V2 build

**STATUS:** BUILD COMPLETE, NO SOLVE. The source was the converged B-ALL-SHORT project; the invalid B-ALL-DIODE archive was not used.

**GATES:** Corrected local-frame cut x −0.30..+0.30, y −0.15..+0.15 mm and sheet x −0.25..+0.25, y −0.15..+0.15 mm were applied at all 250 slotted B sites. At every site, the three mid-copper gap probes at u = −0.2, 0, +0.2 mm found no metal. Each DIODE sheet reads 0.15 mm²; one `DIODE` impedance boundary contains 250 sheets and reads `Resistance=Zr`, `Reactance=Zx`. Metal objects stayed 518. All 266 untouched sites, including 16 gate discs, retained position and volume. Design validation passed before and after CAD history purge. Read-only source/V2 comparisons of the port, absorber, substrate, bias layer, AIR, old boundaries, sphere, Setup1 and mesh operations all passed.

**DATA:** Default `Zr=3.333ohm`, `Zx=-232.667ohm` (ON-A); recorded alternatives ON-B 37.05 / −559.0 and OFF-D 2.667 / +48.667 ohm/square. Gap-cut copper loss ranged 0.00524075–0.00526095 mm³ around the analytic 0.00525 mm³; tolerance was 0.00005 mm³ to allow imported CAD facets, while all 750 gap probes provided the geometry gate. Setup1 remains Direct at 5.68 GHz, min 6/max 12 passes, max ΔS 0.02. No mesh or solve was run.

**FILES:** `BOARD_R150_V31_B_ALL_DIODE_V2.aedtz` is 1,351,047 bytes; SHA-256 `54203F4835F49CA5F14CE0B799B03A96EC5D5119AA7C35103C6A9CF2E83363C3`. The local project is `C:\lwmom_meander_model\hfss_work\full_board_gate_v31_diode_v2\BOARD_R150_V31_B_ALL_DIODE_V2.aedt` (9,814,986 bytes). The archive contains embedded geometry, no STEP references and no results. `audit.json` has all 250 cells and `unchanged_readback.json` records the source/V2 comparison.
