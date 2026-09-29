"""Phase 13 Part 3: the local web app's entry point. ``python -m src.app.main`` / ``lwmom-app`` /
``make app`` all resolve here -- one FastAPI app, a job queue on ``app.state.jobs``, server-rendered
plots (matplotlib, already a dependency of the ``app`` extra), and a static HTML + vanilla-JS
frontend with no build step (see docs/ARCHITECTURE.md "why no frontend build step"). Runs on
localhost only; nothing here calls out to the network.
"""

import webbrowser
from pathlib import Path

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from ..logging_setup import configure_logging
from .jobs import JobQueue

STATIC_DIR = Path(__file__).resolve().parent / "static"


def create_app():
    configure_logging()
    app = FastAPI(
        title="lwmom",
        description="Reflectarray / holographic metasurface antenna (HMSA) design tool -- local, offline.",
        version="0.13.0",
    )
    app.state.jobs = JobQueue()

    from .routers import files as files_router
    from .routers import jobs as jobs_router
    from .routers import library as library_router
    from .routers import unit_cell as unit_cell_router

    app.include_router(jobs_router.router)
    app.include_router(files_router.router)
    app.include_router(unit_cell_router.router)
    app.include_router(library_router.router)

    @app.get("/api/health")
    def health():
        return {"status": "ok"}

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
