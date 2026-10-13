# Architecture

## Module map

Four layers, outermost first:

**(d) App** — `src/app/`: FastAPI backend + static HTML/JS frontend. No physics here; every
router calls an already-tested function below it. See "The app" below for its own internals.

**(c) HMSA design layer** — built on the core solver, used by the app:

| Package | What it does |
|---|---|
| `src/dispersion/` | Surface-wave dispersion (analytic TM0, HFSS-table import, sheet-admittance model; a MoM eigenmode scaffold, not implemented) |
| `src/fields/` | HFSS field-calculator `.fld` import (the monopole reference-wave upload) |
| `src/hologram/` | Hologram synthesis (both formulations), two-state/coded cells, multi-port layouts (`ports.py`: rectangles and convex polygons/hexagons), automatic coded-cell states (`auto_states.py`, Phase 22), polarization compensation, pointing correction/level matching, and the sidelobe optimiser (`sll_optimizer.py` greedy q-norm cell flips; `sll_board.py` job result <-> optimiser, held-out grid; Phase 23) |
| `src/pattern/` | Far-field pattern (scalar legacy + vector polarization-aware model), SLL/XPD/vortex measurement; `cellwise_pattern.py` (Phase 23) the same vector far field as a linear function of each cell's contribution |
| `src/reflectarray/` | Standalone reflectarray designer (phase-gradient lookup, not HMSA) |
| `src/export/` | CSV/JSON/DXF/STEP layout export; `exact_board.py` / `aedt_script.py` / `board_readme.py` (Phase 21-23): exact-coordinate multi-port board STEP, AEDT IronPython script and per-export README |
| `src/dielectric_cell/` | Phase 19: RCWA solver for dielectric (nylon) unit cells and the cell library; **not imported by the app** (model-only, used by `src/transmitarray/` and `examples/`/`scripts/research/phase19_*`) |
| `src/transmitarray/` | Phase 19: horn-fed transmitarray lens design, local-periodicity coupling models, Monte Carlo error bar; imports `src/dielectric_cell/` only, not wired into the app |
| `src/verify/` | Real full-supercell MoM check of the local-periodicity forward model |
| `src/hfss_compare/` | HFSS-vs-MoM single-cell comparison ledger |
| `src/library_sweep.py` | Generalizes `library_builder.py`'s sweep to an arbitrary (incl. grounded) config — not new physics, see its own docstring |

**(b) Research-phase packages** (`src/phase8b/` … `src/phase13/`) — frozen, one package per
phase, **never renamed or edited to suit a new caller** (`CLAUDE.md`). Still imported directly by
the layers above; import-count (modules outside the package itself that import it, as of Phase 17):

| Package | What it is | Imported by |
|---|---|---|
| `phase8b` | Object-wave/multi-beam formulation, reference-wave models | 15 |
| `phase8c` | Unit-cell scattering architecture validation | 9 |
| `phase8d_rev` | Universal unit-cell database architecture | 6 |
| `phase8e` | Control-variable/holographic-scattering bridge, beam reconstruction | 8 |
| `phase9` | Control-variable & local-periodicity validation, HFSS Touchstone reference reader | 23 |
| `phase10` | Environment-aware neighborhood/surrogate/inverse/beam validation | 9 |
| `phase11` | Heterogeneous-MoM nonlocal-coupling study (physics-only, no ML) | 4 |
| `phase12a` | Optional PEC ground plane | 3 |
| `phase12b` | Mesh alignment / duplicate-geometry audit / grounded period study | 0 (analysis-only; its own tests/examples) |
| `phase12c` | Grounded supercells, HFSS export | 7 |
| `phase12d` | General N-port Touchstone reader, aperture Floquet-power test | 5 |
| `phase13` | Richardson mesh-convergence extrapolation, resonance restatement | 1 |

**(a) Core periodic MoM solver** — flat modules directly under `src/`: `config`, `units`,
`constants`, `materials`, `geometry`, `periodic_bc`, `excitation`, `green_function`,
`basis_functions`, `mom_matrix`, `solver`, `sparameters`, `sweep`, `runner`, `database`,
`polarization`, `phase`, plus the library/feasibility modules. Everything above is built on this,
unmodified (the documented, additive exceptions in `CLAUDE.md` "Code conventions" aside).

## Data flow — the app's 3 steps

```
Step 1  upload .fld  --fields.fld_parser-->  reference wave (phase, optional magnitude)
Step 2  beams + board + design method
           --hologram.hmsa_synthesis / two_state_cells / multi_port_synthesis-->  hologram / cells
Step 3  --pattern.vector_aperture_pattern--> far field (SLL/XPD/vortex, pointing/level reports)
        --export.array_export--> CSV/JSON/DXF (+ optional STEP) ready for HFSS
```

Each step is a `POST` that the app either answers directly or runs as a background job
(`src/app/jobs.py`); the browser polls until done. See "The job queue" below.

## The app

`src/app/`: a FastAPI backend plus a static HTML + vanilla-JS frontend, wrapping the `src/*`
service layer above. **No physics lives in `src/app/`** — every router calls an already-tested
`src/*` function and translates its result to JSON.

## Layers

```
Browser (static HTML/JS, src/app/static/)
   | fetch() JSON, multipart uploads
FastAPI routers (src/app/routers/*.py) -- thin HTTP glue, one file per page
   | calls directly, or via the job queue for anything slower than ~1s
   |   (`routers/hmsa.py`: synthesis incl. `/synthesize_ports`, `/auto_states`; `routers/export.py`: board export;
   |    `routers/sidelobes.py`: `/api/sidelobes/optimize`, Phase 23; no router touches `dielectric_cell`/`transmitarray`)
Job queue (src/app/jobs.py) -- ThreadPoolExecutor, in-process, config-hash cache
   | calls
Service layer (src/*, src/dispersion, src/fields, src/hologram, src/reflectarray,
                src/pattern, src/verify, src/export, src/hfss_compare, src/library_sweep.py)
   | calls
Solver core (src/config.py, sparameters.py, mom_matrix.py, ... Phases 1-12D, unmodified)
```

## The job queue (`src/app/jobs.py`)

One process, one `ThreadPoolExecutor`, no Celery/Redis — "keep it simple" per the brief; this is a
single local machine. A router calls `request.app.state.jobs.submit(kind, fn, *args, cache_key=...)`;
`fn(job, *args)` runs in a worker thread and returns a JSON-serializable-after-`to_jsonable` value
(or raises `errors.LwmomError`, whose `.to_dict()` becomes `job.error`; any other exception is
caught too, with a generic `SOLVE_FAILED` code, so a bug in a router never crashes the server).
The browser polls `GET /api/jobs/{id}` (`src/app/static/app.js`'s `LWMOM.pollJob`, every 500 ms)
until the status leaves `queued`/`running`.

**Caching**: `config_hash()` hashes the JSON-serialized request payload; re-submitting an
identical request returns the already-computed (or in-flight) job instead of re-running it. There
is no cache *eviction* — jobs and their results live for the process's lifetime, held in memory
(`JobQueue._jobs`, a plain dict). For a long-running local session with many large sweeps this is
a real, currently-unaddressed memory-growth issue; restarting the app clears it.

**Job chaining without re-uploading data**: `Job.result` holds the *raw* Python/numpy return value
of `fn`, not the JSON-flattened version — only `Job.to_dict(include_result=True)` (the HTTP
response) calls `to_jsonable`. Routers that build on a previous job (Pattern and Export chaining
off a Reflectarray or HMSA design) fetch the *upstream* `Job` object directly from
`request.app.state.jobs` and read its `.result` in-process, instead of asking the browser to
re-upload a potentially large array. This is why the frontend passes a `source_job_id` rather than
the data itself.

**Cancellation** is cooperative (`job.cancelled`, checked between whole solves in a loop — e.g.
`unit_cell.py`'s frequency sweep). No in-flight single MoM solve is preemptible mid-factorization.
The Library sweep's own progress callback does not currently check `job.cancelled` — a known gap,
not wired up (progress is reported, cancellation is not honored mid-sweep there).

## Storage (`src/app/storage.py`, `src/app/library_registry.py`)

Everything the app writes lives under `<app data dir>/app_outputs/` and `<app data dir>/uploads/`,
where `<app data dir>` defaults to `~/.lwmom` (matching `src/logging_setup.py`'s own
`~/.lwmom/logs/`) and is overridable via `LWMOM_APP_DATA_DIR` — every app test sets this to an
isolated `tmp_path` so test runs never touch the real user's home directory. `app_outputs/<job_id>/`
holds whatever files a job produced (STEP exports, layout CSV/JSON, an app-built unit-cell library
JSON); `uploads/<upload_id>/` holds a user-uploaded `.fld`/dispersion-CSV/Touchstone file, read
back by a later job via its ID (`GET /api/files/{job_id}/{filename}` serves the former back to the
browser, with a path-traversal check against the resolved output root).

**Library registry** (`app/library_registry.py`) is a small JSON file (`library_registry.json`
under `app_outputs/`) mapping a `library_id` to its JSON path and a metadata summary (period,
backing, substrate, n_geometries). It auto-seeds the read-only Phase 7 benchmark
(`results/unit_cell_library/database/unit_cell_library.json`) if that file is present on disk
(it's excluded from git, ~49 MB — see `docs/DEVELOPING.md`). This registry is deliberately lighter
than `phase8d_rev.universal_database.UniversalDatabase`: that class answers "which
topologies/periods exist across every registered library" and has no persistence of its own; the
app's actual moment-to-moment need — "give me phase candidates from library X at frequency f" — is
answered directly by `unit_cell_library.UnitCellLibrary`, one library's own query API. **Wiring the
app's libraries into a live `UniversalDatabase` is a documented follow-up, not done here.**

**The HFSS-compare ledger** (`src/hfss_compare/service.py`'s `DEFAULT_LEDGER_PATH`) is the one
exception to "everything lives under the app data dir": it defaults to
`results/hfss_compare_ledger.json`, *inside* the repo tree, by design — it's meant to be a durable
project record of HFSS-vs-MoM comparisons (the brief's own "feed the HFSS gate records"), not
per-user scratch state. Running a real comparison through the app therefore modifies a file `git
status` will show; that's intentional, review it like any other change before committing.

## Frontend (`src/app/static/`)

Plain HTML + vanilla JS, no framework, no build step, one file per page plus shared
`style.css`/`app.js` (`LWMOM.postJson`/`postForm`/`pollJob`/`renderLabels` helpers). **Why no
build step**: this is a local, single-user, offline tool — a bundler/framework buys nothing here
and would add a Node toolchain dependency to a Python-only project whose one deployment target is
"the developer's own machine." Pages that chain off a prior job (Pattern, Export) accept a
`?job=<id>` query-string parameter that a completed design page links to directly.

**Plots** (`src/app/plots.py`, matplotlib, `Agg` backend): the original brief called for
server-rendered PNGs; Part 3 shipped without them (every page rendered an HTML `<table>` instead).
Part 3C built this for the home-page wizard's own flow (reference phase map, layout/hologram
preview, pattern u-v map + per-beam cuts) since that flow explicitly required it — the "More
tools" pages (S-parameter sweeps, phase-vs-geometry tables, HFSS-vs-MoM overlays) still render
tables, not revisited. Two delivery styles: a synchronous endpoint with no job id yet
(`upload_field`) inlines the PNG as base64 in the JSON response; a job-backed endpoint writes the
PNG into `job_output_dir(job.id)` and returns a normal `/api/files/{job_id}/{filename}` URL, the
same mechanism STEP export already used. Every plot call is wrapped best-effort — a matplotlib
failure never fails the real numeric result the user is waiting on.

**The 3-step wizard** (`index.html`, the home page as of Part 3C) is a single-page flow across
three `<section>` panels toggled by JS, not three separate HTML files — it calls the same
`/api/hmsa/*`, `/api/pattern/*`, `/api/export/*` endpoints the old standalone `hmsa.html` (removed)
and the "More tools" pages use. `POST /api/hmsa/ensure_library` hides the Library page from this
flow entirely by building (or reusing, via both the job cache and the persisted `library_registry`)
a small square-patch library for the chosen board/frequency automatically (the reflection-phase
method's own unit-cell data). `POST /api/hmsa/ensure_dispersion` (Phase 14) is its counterpart for
the surface-wave (eigenmode) method: builds a `dispersion.sheet_model` table for the same board and
writes it as a CSV at a deterministic `upload_dir` path (cache hit = file already exists), in the
exact schema an uploaded HFSS eigenmode table uses -- so Formulation A's synthesis code reads either
source identically, with no branch on where the table came from.

## Diagnostics, errors, labels, manifests (Part 2, `src/logging_setup.py`/`errors.py`/`labels.py`/`manifest.py`)

- **Logging**: `create_app()` calls `logging_setup.configure_logging()` once (writes to
  `~/.lwmom/logs/lwmom.log`, JSON-formatted, `LWMOM_LOG_LEVEL`-controlled). `diagnostics.py`'s
  `instrumented_compute_s_parameters` is given `run_id=job.id` by `unit_cell.py`'s solve endpoint,
  so a unit-cell solve's log lines carry the job ID; other routers do not currently pass a
  `run_id` through to the modules they call.
- **Errors**: every router either lets a solver-layer `errors.LwmomError` propagate (the job queue
  translates it to `job.error` with a stable `code`/`message`/`remedy`) or wraps a caught
  lower-level exception (`config.ConfigError`, `ValueError`) into one explicitly
  (`errors.config_invalid(exc, ...)` is the most common). A few endpoints raise a plain
  `HTTPException` instead, for checks that happen *before* a job is even queued (an unknown
  `library_id`, an unfinished source job) — these are HTTP-layer 404/409s, not `LwmomError`s, since
  no job exists yet to attach the error to.
- **Labels**: every design-producing job attaches the fixed vocabulary from `labels.py`
  (`LOCAL_PERIODICITY_APPROXIMATION`, `MESH_CONVERGED`/`MESH_NOT_CONVERGED`,
  `ZERO_THICKNESS_METAL_MODEL`, `STEEP_GRADIENT(...)`, `FEED_REGION_MASKED`, `EXPERIMENTAL`, ...)
  and the frontend renders them as pills (`LWMOM.renderLabels`).
- **Run manifests and the "download diagnostic bundle" action are NOT wired into the app.**
  `src/manifest.py`'s `RunManifest`/`build_diagnostic_bundle` exist and are tested
  (`tests/test_manifest.py`) but no router currently constructs one or exposes a download
  endpoint for it. This is a real gap against the brief ("a run manifest saved with every job...
  a 'download diagnostic bundle' action") — a documented follow-up, not a hidden one.

## Config-hash caching in practice

`unit_cell.py`'s solve endpoint and `hmsa.py`'s `ensure_library`/`ensure_dispersion` endpoints use
`cache_key=config_hash(...)`; other routers (library sweep, reflectarray design, HMSA synthesis,
export) do not cache by config —
re-submitting an identical reflectarray design, for instance, re-runs it. Caching wasn't extended
to every endpoint because most of them are already fast relative to a unit-cell solve (a
reflectarray design over a small library is sub-second) or explicitly write files as a side effect
each time (export, library sweep) where silently returning a stale cached job could be surprising.
