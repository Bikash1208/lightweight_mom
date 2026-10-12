"""Phase 13 Part 3: the local web app's entry point. ``python -m src.app.main`` / ``lwmom-app`` /
``make app`` all resolve here -- one FastAPI app, a job queue on ``app.state.jobs``, server-rendered
plots (matplotlib, already a dependency of the ``app`` extra), and a static HTML + vanilla-JS
frontend with no build step (see docs/ARCHITECTURE.md "why no frontend build step"). Runs on
localhost only; nothing here calls out to the network.
"""

import subprocess
import webbrowser
from importlib.metadata import PackageNotFoundError, version as _package_version
from pathlib import Path

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from ..logging_setup import configure_logging
from .jobs import JobQueue

STATIC_DIR = Path(__file__).resolve().parent / "static"
REPO_ROOT = Path(__file__).resolve().parent.parent.parent


def app_version() -> str:
    try:
        return _package_version("lwmom")
    except PackageNotFoundError:
        return "unknown"


def git_short_sha() -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout.strip() or None if result.returncode == 0 else None


def create_app():
    configure_logging()
    app = FastAPI(
        title="lwmom",
        description="Reflectarray / holographic metasurface antenna (HMSA) design tool -- local, offline.",
        version=app_version(),
    )
    app.state.jobs = JobQueue()

    from .routers import export as export_router
    from .routers import files as files_router
    from .routers import hfss_compare as hfss_compare_router
    from .routers import hmsa as hmsa_router
    from .routers import jobs as jobs_router
    from .routers import library as library_router
    from .routers import pattern as pattern_router
    from .routers import reflectarray as reflectarray_router
    from .routers import sidelobes as sidelobes_router
    from .routers import unit_cell as unit_cell_router
    from .routers import verify as verify_router

    app.include_router(jobs_router.router)
    app.include_router(files_router.router)
    app.include_router(unit_cell_router.router)
    app.include_router(library_router.router)
    app.include_router(reflectarray_router.router)
    app.include_router(pattern_router.router)
    app.include_router(export_router.router)
    app.include_router(hmsa_router.router)
    app.include_router(hfss_compare_router.router)
    app.include_router(verify_router.router)
    app.include_router(sidelobes_router.router)

    @app.get("/api/health")
    def health():
        return {"status": "ok"}

    @app.get("/api/version")
    def version_info():
        return {"version": app_version(), "commit": git_short_sha()}

    # Mounted last: a "/" StaticFiles mount would otherwise shadow every API route registered
    # after it (Starlette matches routes in registration order).
    app.mount("/", StaticFiles(directory=STATIC_DIR, html=True), name="static")
    return app


def main(host="127.0.0.1", port=8000, open_browser=True):
    import uvicorn

    app = create_app()
    if open_browser:
        webbrowser.open(f"http://{host}:{port}")
    uvicorn.run(app, host=host, port=port)


if __name__ == "__main__":
    main()
