"""Phase 13 Part 3: an in-process job queue for the app backend. Solves take seconds to minutes
(CLAUDE.md), so every solve-backed endpoint submits a job and the frontend polls its status
instead of blocking an HTTP request -- "keep it simple" per the brief: one process, one
`ThreadPoolExecutor`, no Celery/Redis (this is a single local machine). Results are cached by a
config hash so re-submitting an identical request returns the already-computed job instead of
re-solving.
"""

import hashlib
import json
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from enum import StrEnum

from ..errors import LwmomError
from .jsonutil import to_jsonable


class JobStatus(StrEnum):
    QUEUED = "queued"
    RUNNING = "running"
    DONE = "done"
    FAILED = "failed"
    CANCELLED = "cancelled"


def config_hash(payload):
    """A stable cache key for a JSON-serializable request payload (order-independent: dict keys
    are sorted before hashing)."""
    blob = json.dumps(to_jsonable(payload), sort_keys=True)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


@dataclass
class Job:
    id: str
    kind: str
    cache_key: str | None
    status: JobStatus = JobStatus.QUEUED
    progress: float = 0.0
    progress_message: str = ""
    result: object = None
    error: dict | None = None
    created_at: float = field(default_factory=time.time)
    started_at: float | None = None
    finished_at: float | None = None
    _cancel_event: threading.Event = field(default_factory=threading.Event, repr=False, compare=False)

    @property
    def cancelled(self):
        """Cooperative cancellation flag: a running job's own function is expected to check this
        between solves in a sweep/loop and stop early. No in-flight single MoM solve is
        preemptible mid-matrix-factorization, so cancellation only takes effect at those
        checkpoints -- documented at each router that loops over a sweep."""
        return self._cancel_event.is_set()

    def request_cancel(self):
        self._cancel_event.set()

    def to_dict(self, include_result=False):
        body = {
            "id": self.id,
            "kind": self.kind,
            "status": self.status.value,
            "progress": self.progress,
            "progress_message": self.progress_message,
            "error": self.error,
            "created_at": self.created_at,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
        }
        if include_result and self.status == JobStatus.DONE:
            body["result"] = to_jsonable(self.result)
        return body


class JobQueue:
    """One instance lives on ``app.state.jobs``, shared by every router. ``submit`` is the only
    entry point routers use: ``fn`` is called as ``fn(job, *args, **kwargs)`` (the running `Job`
    itself, so it can report progress via ``job.progress = ...`` and poll ``job.cancelled``) and
    must return a JSON-serializable-after-``to_jsonable`` value.
    """

    def __init__(self, max_workers=2):
        self._executor = ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="lwmom-job")
        self._jobs: dict[str, Job] = {}
        self._by_cache_key: dict[str, str] = {}
        self._lock = threading.Lock()

    def submit(self, kind, fn, *args, cache_key=None, **kwargs):
        with self._lock:
            if cache_key is not None:
                cached_id = self._by_cache_key.get(cache_key)
                cached = self._jobs.get(cached_id) if cached_id else None
                if cached is not None and cached.status in (JobStatus.DONE, JobStatus.QUEUED, JobStatus.RUNNING):
                    return cached
            job = Job(id=str(uuid.uuid4()), kind=kind, cache_key=cache_key)
            self._jobs[job.id] = job
            if cache_key is not None:
                self._by_cache_key[cache_key] = job.id
        self._executor.submit(self._run, job, fn, args, kwargs)
        return job

    def _run(self, job, fn, args, kwargs):
        job.status = JobStatus.RUNNING
        job.started_at = time.time()
        try:
            job.result = fn(job, *args, **kwargs)
            job.status = JobStatus.CANCELLED if job.cancelled else JobStatus.DONE
        except LwmomError as exc:
            job.status = JobStatus.FAILED
            job.error = exc.to_dict()
        except Exception as exc:  # noqa: BLE001 -- job-boundary translation, mirrors errors.py's own rationale
            job.status = JobStatus.FAILED
            job.error = {"code": "SOLVE_FAILED", "message": str(exc), "remedy": None, "debug_detail": str(exc)}
        finally:
            job.finished_at = time.time()

    def get(self, job_id):
        return self._jobs.get(job_id)

    def cancel(self, job_id):
        job = self._jobs.get(job_id)
        if job is not None:
            job.request_cancel()
        return job

    def list(self):
        return list(self._jobs.values())
