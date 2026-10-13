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

- Active: `PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md` (Parts A–D done; Parts E/F wait on the user's decision whether
  ±30.00° is a hard specification and on an HFSS run of the lens — see `docs/STATUS.md` "Open next steps") and
  `PHASE24_CURATION_AND_CLEANUP.md` (the consolidation pass in progress).
- `_parked/` — deferred ideas, each its own file, not current work: the unit-cell/FSS S-parameter tool (Phase 16) and
  the RIS roadmap + R1 brief (Phase 20, paused, not started).

See `docs/history/briefs/` for every completed brief (phases 12A onward; phases 1–11's own task
descriptions are narrated inline in `docs/history/PROJECT_HISTORY.md` instead, predating this
one-file-per-brief convention).
