"""Phase 13 Part 3: the app backend's in-process job queue -- lifecycle (queued -> running ->
done/failed/cancelled), config-hash cache hits, and typed-error propagation into `job.error`."""

import time

from src.app.jobs import JobQueue, JobStatus, config_hash
from src.errors import solve_failed


def _wait(job, timeout=5.0):
    deadline = time.monotonic() + timeout
    while job.status in (JobStatus.QUEUED, JobStatus.RUNNING):
        if time.monotonic() > deadline:
            raise TimeoutError(f"job {job.id} did not finish within {timeout}s")
        time.sleep(0.01)
    return job


def test_submit_runs_the_function_and_reports_the_result():
    q = JobQueue()
    job = q.submit("echo", lambda _job, x: x * 2, 21)
    _wait(job)
    assert job.status == JobStatus.DONE
    assert job.result == 42


def test_failing_job_reports_a_stable_error_code():
    q = JobQueue()

    def _boom(_job):
        raise solve_failed(message="deliberate test failure")

    job = q.submit("boom", _boom)
    _wait(job)
    assert job.status == JobStatus.FAILED
    assert job.error["code"] == "SOLVE_FAILED"
    assert job.error["message"] == "deliberate test failure"


def test_unexpected_exception_still_produces_an_error_dict_not_a_crash():
    q = JobQueue()
    job = q.submit("boom", lambda _job: 1 / 0)
    _wait(job)
    assert job.status == JobStatus.FAILED
    assert job.error is not None and "code" in job.error


def test_cache_key_returns_the_same_job_without_rerunning():
    q = JobQueue()
    calls = []

    def _fn(_job):
        calls.append(1)
        return "result"

    key = config_hash({"a": 1, "b": 2})
    job1 = q.submit("kind", _fn, cache_key=key)
    _wait(job1)
    job2 = q.submit("kind", _fn, cache_key=key)
    assert job2.id == job1.id
    assert len(calls) == 1  # not re-run


def test_config_hash_is_order_independent_and_stable():
    assert config_hash({"a": 1, "b": 2}) == config_hash({"b": 2, "a": 1})
    assert config_hash({"a": 1}) != config_hash({"a": 2})


def test_cancel_sets_the_cooperative_cancel_flag():
    q = JobQueue()
    seen_cancelled = []

    def _fn(job):
        # A real sweep checks job.cancelled between solves; simulate one checkpoint.
        time.sleep(0.05)
        seen_cancelled.append(job.cancelled)
        return "done"

    job = q.submit("kind", _fn)
    q.cancel(job.id)
    _wait(job)
    assert seen_cancelled == [True]
    assert job.status == JobStatus.CANCELLED


def test_get_and_list_and_cancel_of_unknown_id():
    q = JobQueue()
    assert q.get("does-not-exist") is None
    assert q.cancel("does-not-exist") is None
    job = q.submit("kind", lambda _job: None)
    _wait(job)
    assert job in q.list()


def test_job_to_dict_omits_result_unless_asked_and_only_when_done():
    q = JobQueue()
    job = q.submit("kind", lambda _job: {"value": 1})
    _wait(job)
    assert "result" not in job.to_dict()
    assert job.to_dict(include_result=True)["result"] == {"value": 1}
