# lightweight_mom (lwmom)

A lightweight Python Method-of-Moments solver for **periodic metasurface unit cells** (not a
general-purpose full-wave solver), plus a local web app for reflectarray / holographic-metasurface
(HMSA) design built on top of it.

```
Unit-cell geometry -> Periodic MoM -> Surface current -> Reflected/transmitted fields
                    -> S11/S21 (mag, phase) -> Unit-cell database -> Holographic/reflectarray synthesis
```

**Status:** Phases 1-12D (solver core, validation, grounded supercells, real-HFSS ingestion) and
Phase 13 (production repo + local web app) done. The full phase-by-phase narrative, module table
and the project's assumptions/resolved/still-open lists live in
**[docs/PROJECT_HISTORY.md](docs/PROJECT_HISTORY.md)** -- read that for the science; this file is
just the front door. Guidance for an AI assistant (or a new contributor) working in this repo:
**[CLAUDE.md](CLAUDE.md)**.

## Install

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,app]"       # add ",cad" for STEP export / plotting (cadquery, matplotlib)
```

Python >= 3.11 required (developed and tested on a 3.14 venv). A pinned, known-working environment
is recorded in `requirements-lock.txt` (`pip install -r requirements-lock.txt` reproduces it
exactly). See `docs/DEVELOPING.md` for the fresh-clone verification steps.

## Run the app

```bash
make app          # or: lwmom-app   (installed console script, same thing)
```

Starts the local FastAPI backend and opens the browser UI at `http://127.0.0.1:8000`. Runs
entirely on `localhost`; no external services; data never leaves the machine (see
`docs/ARCHITECTURE.md`). Every design it produces is explicitly a lookup-table /
local-periodicity design to be verified in full-wave -- the app labels outputs accordingly, it
does not hide the project's open scientific questions (below).

## Run the tests

```bash
make test-fast     # or: pytest -m "not slow"   -- excludes real-MoM-solve tests, < 5 min
make test-full      # or: pytest                 -- everything, tens of minutes
```

Tests needing the optional `cad` extra, the ~470 MB Phase 7 checkpoint, or real HFSS export files
not tracked in git (`hfss_exports/**/*.s36p`, `*.s4p`) skip cleanly when those are absent
(`cad`/`checkpoint`/`hfss_data` pytest markers; see `docs/DEVELOPING.md` for what's excluded from
git and how to regenerate or obtain it).

## Before Phase 13's app is trusted for a real design -> still open

These are scientific findings, not bugs, and the app surfaces them as validation labels rather
than hiding them (`docs/ARCHITECTURE.md` "Validation labels"):

- **Local-periodicity step limit** (design-relevant metric, Phase 12D): holds to 5.8 deg/cell
  off-resonance; unresolved near/at resonance (genuinely infeasible to sample finely there --
  see `docs/phases/PHASE12D_LOCAL_PERIODICITY.md`).
- **Lookup-table forward model at the aperture level** (Phase 12D Task B + Phase 13 Part 0):
  no tested gradient size (N=6/8/12) reaches 90% of reflected power in the intended order, but
  Phase 13's own idealized (no-coupling) prediction from the SAME lookup table reproduces the
  measured full-MoM result to within ~2 percentage points -- i.e. this design's 186 deg (of 360
  deg) achievable phase span, not the lookup-table/local-periodicity forward model itself, is the
  limit (`docs/phases/PHASE13_WEBAPP.md`).
- **HFSS single-cell validation**: off-resonance (G1) matches to ~0.3 deg mean phase; at resonance
  (G4) a ~3-4 deg residual remains unexplained after accounting for mesh convergence (~5.9 deg) and
  patch thickness (~3.7 deg) (`docs/phases/PHASE13_WEBAPP.md`, carrying `PHASE12D_ADDENDUM.md`
  forward). G2/G3/G5 still have no real HFSS run at all.
- **Geometry inversion and full-array STEP export remain `GEOMETRY_INVERSION_LOCKED` /
  `STEP_EXPORT_LOCKED`** as scientific-finding records (Phase 12D), even though Phase 13's app now
  *produces* designs and exports anyway, per the user's explicit Phase 13 decision -- every such
  output carries these labels, it does not claim the open questions are resolved.
- **Surface-wave-dispersion-based HMSA control** (the physically correct per-cell quantity for a
  guided-wave hologram, as opposed to normal-incidence S11 phase) has no MoM capability yet; Phase
  13 ships HFSS-eigenmode-table import only (`DISPERSION_FROM_HFSS`) and leaves a MoM eigenmode
  solver as a documented, tested scaffold (`DISPERSION_FROM_MOM`, not implemented).

## More

- `docs/ARCHITECTURE.md` -- backend/job-queue/UI/solver layers, data flow, where logs and run
  manifests go, and why the app has no frontend build step.
- `docs/USER_GUIDE.md` -- walkthrough of every app page.
- `docs/DEVELOPING.md` -- setup, test markers, adding a solver capability, reading a diagnostic
  bundle, what's excluded from git and how to regenerate/obtain it.
- `docs/theory/` -- `CONVENTIONS.md` (single source of truth for signs/units), `GREEN_FUNCTION.md`,
  `MOM_FORMULATION.md`, `S_PARAMETERS.md`, `DATABASE.md`.
- `docs/phases/` -- one write-up per phase (the actual science). `docs/briefs/` -- the original
  task briefs each phase was executed from.
- `CHANGELOG.md` -- per-phase summary of what changed. `CONTRIBUTING.md` -- how to contribute.
