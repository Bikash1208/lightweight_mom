# HFSS export packages

Each subfolder is a self-contained HFSS import package (STEP + layout CSV/JSON/DXF + an
`HFSS_README.md` with Floquet-port/PEC/thickness setup notes) for one design, plus whatever
real HFSS results exist for it. Raw HFSS Touchstone exports (`*.s36p`, `*.s4p`) are gitignored —
`docs/DEVELOPING.md` has the full account of what's excluded and why.

| Folder | What | Status |
|---|---|---|
| `phase12_grounded_8mm/` | Single-cell (G1–G6) and supercell (S1–S3) grounded 8 mm unit cells | **Validated**: G1 off-resonance (mean phase diff -0.29°); G4 at-resonance (13.51° diff, decomposed into mesh/thickness/residual — `docs/validation/README.md`). G2/G3/G5 and all supercells await an HFSS run. |
| `two_beam_20deg_phi0_180/` | Two-beam (20°, phi=0/180°) single-aperture board, 150 mm and 300 mm | Awaiting HFSS. |
| `three_port_25deg/` | 1×3 port board, 25° RHCP beams per port (`scripts/adhoc/three_port.py`) | **Awaiting HFSS** (`docs/STATUS.md` "Open next steps"). `*.step` gitignored; everything else committed. |
| `four_port_1bit_gap0/`, `four_port_1bit_gap28/` | 4-port 1-bit coded boards, two guard-gap variants (`scripts/adhoc/coded4.py`) | **Awaiting HFSS** (`docs/STATUS.md` "Open next steps"). `*.step` gitignored; everything else committed. |
| `fss_square_loop_10mm/` | Square-loop FSS unit cell, periodic MoM + finite-11×11-panel plane-S (`scripts/adhoc/fss_loop.py`, `scripts/adhoc/fss/plane_s.py`) | **Self-consistency reference, not HFSS-validated** — no HFSS data exists for this geometry (`docs/validation/README.md` item 5). No `*.step` (CSV/PNG only). |

## `*.step` files and git

`hfss_exports/**/*.step` is gitignored **except** the STEP files already tracked before this rule
was added (`phase12_grounded_8mm/`, `two_beam_20deg_phi0_180/`) — those stay tracked. Any new
package's own `.step` (large, binary, fully regenerable from its `layout.json` via the app's
Export page or the generating script above) is not committed; everything else in the package
(README, CSV, JSON, PNG) is.
