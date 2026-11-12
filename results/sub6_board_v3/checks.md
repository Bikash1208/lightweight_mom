# Board V3 exporter checks (numeric, from manifest.json and the re-imported STEP)

Board: circular, R = 180 mm (D = 360 mm), rot180 frame, RT5880 3.175 / ground / RT5880 bias layer 0.254 mm, no FR4, vias only (no bias network). 696 sites: 636 STATIC (A) + 60 SWITCHED (B), 60 diodes. Exporter code did not launch AEDT.

| Check | Result | Numbers |
|---|---|---|
| No overlapping copper (top layer, 756 copper bodies; bottom pads; barrels vs ground) | **PASS** | overlapping pairs: 0; pad overlaps: 0; barrels touching ground: 0 |
| Copper gaps >= 0.20 mm | **PASS** | smallest gap 0.2988 mm = the 0.30 mm diode gap (['disc_294', 'island_294'], 0.001 mm coordinate grid); other kinds: disc_static|disc_static 0.800; disc_static|disc_switched 0.765; disc_switched|disc_switched 0.730; disc_switched|island 0.299; bottom pads: none within 2 mm of another; anti-pads: none within 2 mm of another; anti-pad to feed opening 22.38 mm |
| Via-pad margins >= 0.10 mm | **PASS** | island via pad margin min 0.1092 mm; disc via pad edge inside disc edge min 0.2013 mm (>= 0.20 required, 0.202 nominal); barrel margin min 0.2592 mm; every pad inside its host: True |
| Every anti-pad clears its via barrel >= 0.35 mm radial | **PASS** | anti-pad r 0.50 - barrel r 0.15 = 0.3500 mm; concentric for all 120 vias |
| No switched node (island or disc) DC-connected to the ground | **PASS** | ground component = ['GROUND']; barrels touching the ground: 0; pin-to-ground distance 1.415 mm |
| Each island is connected to its disc only through its diode sheet | **PASS** | islands in a different copper component from their own disc: True; island-to-disc min gap 0.2988 mm, overlap area 0.0; sheet corners on island/slot polygon vertices: 60/60; sheet sizes (mm, 0.001 grid): [[0.299, 0.599], [0.299, 0.6], [0.3, 0.599], [0.3, 0.6], [0.3, 0.601], [0.301, 0.6]] |
| Diode count = 60 | **PASS** | sheets in manifest 60; switched cells 60; diode products in the re-imported STEP 60 |
| Nothing in the 12 x 12 mm SMA window | **PASS** | items in window: []; nearest top copper 12.65 mm from the window edge; nearest anti-pad 16.46 mm |
| Nothing across the board edge (R = 180 mm) | **PASS** | largest copper radius 179.400 mm, margin 0.600 mm |
| Bottom layer: only via pads (no bias lines / resistors / connector pads) | **PASS** | 120 pads, features ['via_pad'] |
| STEP re-import | **PASS** | products 66, solids 761 (expected 761: 756 copper lumps + ground + pin + substrate + bias + air), 0 inexact literals, bbox [-195.0, -195.0, -18.499, 195.0, 195.0, 40.035] |

## Layout and cell provenance
- Partition: `jobs1_iter1_study.py` (k = 1, N = 60, R = 180 mm, ring exclusion 25 mm around the feed), placement chosen by the study's own key: **inner**, static-cell cap 0.5 Np/m, rot180 frame (beam A at phi 180, B at phi 0). Study read-out for this layout (model only): A/B OFF -9.9 / -17.3 dB (balance A-B +7.37 dB), B suppression ON 8.8 dB, converted power OFF 0.203 / ON 0.118. The study could NOT meet the 3 dB balance with 60 diodes at this radius (best case shown).
- Switched cell: `island`, LS5p5 + WS2p4 (slot 5.5 x 2.4 mm, island 3.5 x 1.5 mm, disc r 5.435 mm, diode gap 0.30 mm, opposite clearance 0.60 mm), disc via pad edge 0.202 mm inside the disc edge, island via at the island centre-line, 1.34 mm from the island centre; rotated by (slot_axis - 90) deg about each cell centre (canonical slot along +y).
- Static A cells: the existing NOCUT table rows (slot L x W per cell, axis from the design). **FLAG: FR4-calibrated**, not re-derived for the RT5880 bias-layer stack.
- Diode sheet 0.30 x 0.60 mm at z = 0.0175 mm; anode edge on the island, cathode edge on the disc slot wall; the sheet impedances are set per design in the HFSS project (ELS values); `diode_sheets.json` lists the old OFF 22.23 - j335.4 / ON 1.60 + j29.2 Ohm/sq.
- Feed: centre pin r 0.635 mm (z -3.175 ... +8.925), ground opening r 2.05 mm, port annulus at z -3.175 mm (feed block from `results/sub6_feed/fix1/feed_block.json`); the 12 x 12 mm SMA window is cut through the bias layer only.
- Coordinates are rounded to the 0.001 mm grid after rotation, so nominal 0.300 / 0.600 sheet sizes appear as 0.299-0.301 / 0.599-0.601 mm.
