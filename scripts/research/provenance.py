"""Provenance record for every result produced on any machine (infra brief rule 3): machine, OS, git commit and dirty flag, Python / numpy / scipy versions, BLAS vendor, thread settings,
worker count, parameters, start time. `provenance(extra)` returns a JSON-serialisable dict; `finish(rec)` adds the end time.
"""

import datetime
import os
import platform
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _git(*args):
    try:
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception:  # noqa: BLE001
        return "NOT_AVAILABLE"


def _blas():
    try:
        import numpy as np

        cfg = np.show_config(mode="dicts")
        b = cfg.get("Build Dependencies", {}).get("blas", {})
        return {"name": b.get("name"), "version": b.get("version")}
    except Exception:  # noqa: BLE001
        return "NOT_AVAILABLE"


def provenance(extra=None, workers=None):
    import numpy as np
    import scipy

    rec = {
        "machine": platform.node(),
        "os": f"{platform.system()} {platform.release()} {platform.machine()}",
        "cpu_count": os.cpu_count(),
        "git_commit": _git("rev-parse", "HEAD"),
        "git_branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "git_dirty": bool(_git("status", "--porcelain")),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "blas": _blas(),
        "threads": {
            k: os.environ.get(k)
            for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
        },
        "workers": workers,
        "started": datetime.datetime.now(datetime.UTC).isoformat(),
    }
    rec.update(extra or {})
    return rec


def finish(rec):
    rec["finished"] = datetime.datetime.now(datetime.UTC).isoformat()
    return rec
