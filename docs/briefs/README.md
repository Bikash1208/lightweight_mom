# Briefs — the working method

A **brief** is a single markdown file describing one phase or task: what to build, what not to
touch, how to verify it, and (for anything that could lose data) explicit safety rules. The user
writes or pastes a brief; the implementer (human or AI agent) reads it, works through it in order,
and reports back against its own acceptance checks.

## Conventions

- **One commit per logical item** in the brief, fast suite (`make test-fast`) green before each.
- **No silent physics/solver changes** unless the brief explicitly says so for this phase.
- **Never move or delete a release tag.** A brief that needs a safety tag says so explicitly
  (e.g. "tag `pre-cleanup-<date>` first") — follow that before touching anything else.
- When a brief is finished (its own acceptance checks pass, merged to the target branch), move it
  to `docs/history/briefs/` (`git mv`, keep the name) so `docs/briefs/` only ever shows what's
  still active.
- An idea that isn't ready to become a brief, or that was explicitly deferred, lives in
  `docs/briefs/_parked/` instead — not deleted, not acted on until someone writes a real brief for it.

## What's here now

- **Primary active track (Phase 28 — Monopole-Fed Coded Holographic Metasurface Antenna):**
  - `PHASE28_CODED_HOLOGRAPHIC_ANTENNA.md`: Monopole-fed coded (PIN-diode reconfigurable) holographic metasurface antenna. Direction change 2026-10-24: reflectarray / plane-wave RIS designs with external feed dropped; focus on integrated monopole feed launching surface wave across PIN-diode switched holographic cells. Current status: items 0–2 done, at the item-2 stop point (awaiting the user's decision).
  - `INFRA_COMPUTE_SERVER_SYNC.md`: Architecture and sync rules for remote Windows compute servers.
- **Active on hold:**
  - `PHASE26_LEAKY_TENSOR_CP_HEXAGON.md`: Leakage-calibrated tensor-cell pure-CP hexagon (synthesis and Board A export done; HFSS simulation run parked per `docs/STATUS.md`).
- **Parked:**
  - `_parked/PHASE20_R5_RIS_DESIGN_MODE.md`: 2-bit X-band PIN-diode RIS design mode (parked 2026-10-24; R5-1/2/3/5 and R5-6 infra done, R5-4 stopped; direction change to monopole-fed holographic antenna).
  - `PHASE20_RIS_ROADMAP.md`: High-level roadmap for the RIS track (R1–R4 completed and archived in `docs/history/briefs/`; R5 parked).
  - `_parked/PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md`: 28 GHz nylon lens (Parts A–D done; Parts E/F parked).
  - `_parked/PHASE25_BUTLER_SECTOR_SCANNER.md`: Phase 25 Butler sector scanner (paused).
  - `_parked/PHASE16_UNIT_CELL_SPARAM_TOOL.md`: Phase 16 unit-cell/FSS S-parameter tool.

See `docs/history/briefs/` for every completed brief (phases 12A onward; phases 1–11's own task
descriptions are narrated inline in `docs/history/PROJECT_HISTORY.md` instead, predating this
one-file-per-brief convention).
