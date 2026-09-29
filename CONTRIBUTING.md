# Contributing

This is an ISRO-internal research project (not currently accepting external contributions), but
the workflow below is the one to follow for any change, internal or not.

## Setup

```bash
git clone <repo> && cd lightweight_mom
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,app,cad]"
pre-commit install
make test-fast
```

See `docs/DEVELOPING.md` for the full walkthrough (markers, adding a solver capability, reading a
diagnostic bundle) and what's excluded from git / how to regenerate or obtain it.

## Rules that apply to every change (see `CLAUDE.md` for the full list)

1. **Read `docs/theory/CONVENTIONS.md` first.** If code and it disagree, that's a bug: fix one and
   add a test.
2. **Never fabricate data.** No invented HFSS/measured results -- report
   `NOT_AVAILABLE_FROM_CURRENT_DATA` instead. Report measured numbers, don't fit decay laws or
   models before checking the data.
3. **Reuse earlier phases/modules unmodified** (wrap or duck-type). Editing a frozen core module
   to suit a new feature needs an explicit, documented reason (the project's own "deliberate
   exceptions" list in `CLAUDE.md`) and must stay bit-identical for every pre-existing
   config/sweep/result.
4. **Gates return explicit named states**, never a bare pass/fail; don't upgrade a status without
   evidence.
5. **Tests check mathematical/numerical properties** (reciprocity, energy conservation, symmetry,
   convergence, known-answer cases), not just "it runs." Add tests with every change.
6. **No ML / fitted surrogates** in the solver core or `src/phase11/` (a test greps for this).
7. Run `make lint` and `make test-fast` before opening a PR; `make test-full` before a release.
   Logical, reviewable commits with clear messages.

## Reporting a bug

Use the app's "download diagnostic bundle" action (manifest + logs + inputs) if the bug happened
in the web app -- it captures everything needed to reproduce without you needing to describe the
exact steps from memory.
