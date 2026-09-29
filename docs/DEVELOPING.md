# Developing

## Setup

```bash
python3 -m venv .venv && .venv/bin/pip install -e ".[dev,app,cad]"
# or, to reproduce the exact pinned environment this project was developed against:
.venv/bin/pip install -r requirements-lock.txt
```

`make setup` does the venv + install + `pre-commit install` in one step.

## Running tests

```bash
make test-fast   # pytest -m "not slow"  -- 2090+ tests, ~3.5 min on this machine
make test-full   # pytest                -- everything, ~35-40 min (real MoM solves)
.venv/bin/python -m pytest tests/test_mom_matrix.py -q   # a single file
```

Run from the project root (`pytest.ini` sets `pythonpath = .`).

### Markers (`pytest.ini`)

- `slow` — a real MoM solve (or several) expensive enough to exclude from the fast suite.
  Retrofitted onto ~89 pre-existing tests in this phase (see the "Retrofit @pytest.mark.slow"
  commit) as well as Phase 13's own new expensive tests. When adding a new test that does a real
  solve taking more than a couple of seconds: if it shares a module-scoped (or session-scoped)
  fixture with other tests in the same file, mark the **whole file**
  (`pytestmark = pytest.mark.slow` near the top) — excluding only the one test that happened to
  trigger the fixture doesn't save anything, since another unmarked test in the file would just
  trigger it instead. If the expense is local to one test function, mark only that function.
- `cad` — needs the optional `cadquery`/`matplotlib` extras (STEP export, plots).
- `checkpoint` — needs the ~470 MB Phase 7 unit-cell-library checkpoint on disk (see below).
- `hfss_data` — needs real HFSS exports under `hfss_exports/` (not regenerable without HFSS).

All three self-skip cleanly (`pytest.mark.skipif(not <path>.exists(), ...)`) when their data isn't
present, rather than failing — a fresh clone's fast suite genuinely excludes almost nothing by
skip count; it's the `slow` marker that does the actual exclusion.

### App tests specifically

Every `tests/test_app_*.py` file does `pytest.importorskip("fastapi")` at the top (skips cleanly
without the `app` extra) and uses a `client` fixture that sets `LWMOM_APP_DATA_DIR` to an isolated
`tmp_path` via `monkeypatch` — **always use this pattern for a new app test**; without it, the
test would write real files (job outputs, uploads, the library registry) into the actual user's
`~/.lwmom`, and repeated runs would accumulate state across test sessions (a real bug this phase
hit once — see the "unit-cell and library pages" commit). The one deliberate exception is
`src/hfss_compare/service.py`'s ledger, which defaults to `results/hfss_compare_ledger.json`
*inside* the repo by design (`docs/ARCHITECTURE.md`); a test that exercises it should still pass
an explicit `path=` (or `monkeypatch.setattr` the module's `DEFAULT_LEDGER_PATH` attribute, which
works because that default is resolved at call time — see the "HFSS compare and verify pages"
commit for why it wasn't always safe to monkeypatch).

Tests that start a real live server (a smoke test, not just `TestClient`) are not part of the
suite — they were run manually during development (`python -c "from src.app.main import main;
main(port=..., open_browser=False)"` with `LWMOM_APP_DATA_DIR` set, then plain `curl`/`urllib`) to
catch frontend/backend field-name mismatches a `TestClient`-only test can't, since nothing actually
parses the JS. If you change a JSON response shape a page's JS reads, re-run that page's flow this
way before trusting it — this phase caught a real bug exactly this way (`to_record()`'s magnitude
is linear, not dB; the unit-cell page was reading a `_mag_db` key that doesn't exist).

## Linting

```bash
make lint     # ruff check + ruff format --check
make format   # ruff check --fix + ruff format
```

`pyproject.toml`'s `[tool.ruff.lint]` has a pragmatic ignore list (documented there) plus
`extend-immutable-calls` for `fastapi.File`/`Form`/`Depends`/`Query` — FastAPI's own required
dependency-injection idiom (`param: T = File(...)`) otherwise trips `B008` ("no function calls in
argument defaults"), which is a real bug pattern everywhere *except* this one framework
convention.

## What's excluded from git and how to regenerate or obtain it

- `results/unit_cell_library/{_solve,database,geometry}/` — the Phase 7 checkpoint and its
  derived exports (~554 MB). Regenerate: `examples/build_unit_cell_library.py` (resumable,
  appends to `checkpoint.jsonl`; takes a long time — 324 geometries × 201 frequencies).
- `hfss_exports/**/*.s36p`, `*.s4p` — real HFSS Touchstone exports (~24 MB each). Not regenerable
  from this repo's own code; obtain from wherever the original HFSS runs are archived. `hfss_data`
  -marked tests skip cleanly without them.
- `.venv/`, `__pycache__/`, `.pytest_cache/`, `.ruff_cache/`, OS/IDE files — standard, regenerate
  via `make setup` / normal tool usage.
- The app's own runtime state (`~/.lwmom/`: logs, job outputs, uploads, the library registry) is
  never part of the repo at all (`docs/ARCHITECTURE.md`).

## Adding a solver capability

Follow the established pattern (CLAUDE.md "Code conventions"): new phases/capabilities reuse
earlier, frozen modules **unmodified** — wrap or duck-type, never edit a frozen module's existing
behavior to suit a new caller. Phase 13's own two exceptions, both additive and both already
covered above:
- `src/phase13/convergence.py`'s Richardson-extrapolation option (the one permitted additive
  numerics change, default only inside the app).
- `src/library_sweep.py` generalizes `library_builder.py`'s sweep orchestration to an arbitrary
  config (in particular, a grounded one) rather than editing `library_schema.LibraryConfig` to
  grow a `backing` field it was deliberately never given — read that module's own docstring for
  the reasoning before adding a third way to build a library sweep.

To add a new app page: one router file under `src/app/routers/`, included in
`src/app/main.py::create_app()`; if it wraps a real solve, submit it as a job
(`request.app.state.jobs.submit(...)`) rather than blocking the request; return
`to_jsonable(...)`-safe values (or rely on `Job.to_dict`'s own call to it) — watch specifically for
non-finite floats (`inf`/`nan`), which `to_jsonable` maps to `null` (a bare `Infinity` token is
invalid JSON that breaks browser `JSON.parse` — this phase hit that with `verify.supercell_check`'s
refused-size cost estimate). Add a `src/app/static/<page>.html` following the existing pages'
structure (shared nav, `style.css`, `app.js`'s `postJson`/`postForm`/`pollJob`/`renderLabels`).

## Reading a diagnostic bundle

**Not currently produced by the app** (`docs/ARCHITECTURE.md` "Diagnostics, errors, labels,
manifests" — a documented gap). `src/manifest.py`'s `build_diagnostic_bundle` exists and is
tested (`tests/test_manifest.py`): it zips a `RunManifest`'s JSON, the slice of the log file
matching that manifest's `run_id`, and any input files, into `diagnostic_bundle_<run_id>.zip`. To
wire it into a page: construct a `RunManifest` in that router's job function (inputs, resolved
config, labels/warnings, timings), call `manifest.add_output_file(...)` for each file the job
produces, and expose a download endpoint the same way `routers/files.py` already does for job
outputs.
