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

- Active: whatever hasn't been merged/archived yet (check `docs/STATUS.md` for the current one).
- `_parked/` — deferred ideas, each its own file, not current work.

See `docs/history/briefs/` for every completed brief (phases 12A onward; phases 1–11's own task
descriptions are narrated inline in `docs/history/PROJECT_HISTORY.md` instead, predating this
one-file-per-brief convention).
