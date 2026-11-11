"""Validate a built Stage 1 project in its own AEDT 2025.2 instance."""
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

from ansys.aedt.core import Hfss

project = Path(sys.argv[1]).resolve()
out = Path(sys.argv[2]).resolve()
out.mkdir(parents=True, exist_ok=True)
report = out / (project.stem + "_validation.json")
d = {"project": str(project), "own_python_pid": os.getpid(),
     "started_utc": datetime.now(timezone.utc).isoformat()}
report.write_text(json.dumps(d, indent=2))
h = None
try:
    h = Hfss(project=str(project), design="AnchorPreflight", solution_type="Modal",
             version="2025.2", non_graphical=True, new_desktop=True,
             close_on_exit=False)
    d["aedt_version"] = str(h.desktop_class.aedt_version_id)
    d["own_aedt_pid"] = h.desktop_class.aedt_process_id
    report.write_text(json.dumps(d, indent=2))
    messages, ok = h.validate_full_design(output_dir=str(out), ports=1)
    d["messages"] = messages
    d["ok"] = bool(ok)
except Exception as exc:
    d["error"] = repr(exc)
    d["ok"] = False
finally:
    d["ended_utc"] = datetime.now(timezone.utc).isoformat()
    report.write_text(json.dumps(d, indent=2, default=str))
    if h is not None:
        h.release_desktop(close_projects=True, close_desktop=True)
