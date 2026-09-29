"""Phase 13 Part 3: where the app's own generated files live -- STEP exports, layout CSV/JSON,
app-built unit-cell libraries, diagnostic bundles, and the library registry
(`app/library_registry.py`). One directory per job under ``<app data dir>/app_outputs/`` (matches
`src/logging_setup.py`'s ``~/.lwmom/logs/`` convention): never inside the repo tree, so
``results/`` stays exactly what the test suite and build scripts expect (CLAUDE.md).

The app data dir defaults to ``~/.lwmom`` but is resolved fresh on every call (not cached at
import time) so tests can point it at an isolated ``tmp_path`` via the ``LWMOM_APP_DATA_DIR``
env var instead of writing real files into the user's home directory on every test run.
"""

import os
from pathlib import Path


def app_data_dir():
    override = os.environ.get("LWMOM_APP_DATA_DIR")
    return Path(override) if override else Path.home() / ".lwmom"


def output_root():
    return app_data_dir() / "app_outputs"


def job_output_dir(job_id):
    d = output_root() / job_id
    d.mkdir(parents=True, exist_ok=True)
    return d
