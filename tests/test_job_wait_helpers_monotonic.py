"""Phase 20 R4 item 0: job-waiting test helpers must use the monotonic clock. The machine's wall clock jumps (observed: `time.time()` advanced 903 s during a
2 s sleep), which tripped 120-300 s wait timeouts after a few seconds."""

import re
import time
from pathlib import Path

from src.app.jobs import JobQueue, JobStatus

TESTS = Path(__file__).resolve().parent


def test_no_job_wait_helper_reads_the_wall_clock():
    bad = []
    for p in sorted(TESTS.glob("test_app_*.py")) + [TESTS / "test_ui_coded_panel.py"]:
        for i, line in enumerate(p.read_text().splitlines(), 1):
            if re.search(r"time\.time\(\)", line):
                bad.append(f"{p.name}:{i}: {line.strip()}")
    assert not bad, "wall-clock reads in job-wait helpers:\n" + "\n".join(bad)


def test_a_wall_clock_jump_does_not_trip_the_wait_timeout(monkeypatch):
    from tests.test_app_jobs import _wait

    jumps = iter(range(0, 10**9, 10**6))
    monkeypatch.setattr(time, "time", lambda: next(jumps))  # every read jumps by 1e6 s
    q = JobQueue()
    job = q.submit("slow", lambda _job: (time.sleep(0.2), 7)[1])
    assert _wait(job, timeout=5.0).status == JobStatus.DONE
