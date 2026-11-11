# Option 1 lab B: top-layer bias cell

AEDT 2025.2 / PyAEDT 1.9.0. BAR90-02LRH sheet 0.30 x 0.60 mm; OFF `Rs=22.23ohm`, `Xs=-335.4ohm`; ON `Rs=1.60ohm`, `Xs=29.2ohm`.

## Geometry and local copper DC check

All ten geometries passed the analytic pre-build check and AEDT design validation. The island and disc remain separate top-copper solids, with no local copper contact except across the diode sheet's 0.30 mm span. The line leaving the island occupies the opened channel; a copper receiver tab at the opposite periodic boundary connects the incoming line to the local disc rim. The disc rim ends 0.365 mm inside the 11.6 mm cell boundary, so that tab is required for the stated series string.

| Variant | Line width | Channel width | Clearance on each side |
|---|---:|---:|---:|
| REF | none | none | n/a |
| Y/X narrow | 0.10 mm | 0.40 mm | 0.15 mm |
| Y/X wide | 0.20 mm | 0.50 mm | 0.15 mm |

The 0.20 mm variant changes **both** the line and channel width, per the latest instruction. Each variant has island lengths 4.348 and 2.500 mm. The stack has top copper, 3.175 mm RT5880, and a solid PEC ground at the model bottom; it has no vias, antipads, FR4, stubs, or air below ground. Mesh operations: none. See `geometry_precheck.json` and `builds/`.

## Fresh REF determinism VERIFY

Two independently rebuilt `REF_L4.348` projects solved OFF / along / 90 degrees using new local AEDT instances. Each yielded eight modes. All eight frequencies and Q values matched exactly in the raw PyAEDT outputs; both variable readbacks passed. The gate uses 0.5% frequency and 5% Q tolerances and no external reference. See `verify/determinism.json`, the two row JSONs, and HFSS profiles.

QUEUE_B Option 1: 120 planned rows (10 geometries x OFF/ON x along/across x 80/90/100 degrees). The runner measures free CPU before each geometry, uses up to five new local AEDT processes at 8 cores each, retries a failing row once after a reconnect, and pushes results after about 24 rows. Analysis will exclude Q >= 1000 and Q <= 60 modes.
