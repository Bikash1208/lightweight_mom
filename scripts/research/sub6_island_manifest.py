"""Finalize model references, provenance, validation and delivery manifests."""
import contextlib
import hashlib
import io
import json
import platform
import subprocess
import sys
from datetime import datetime,timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).resolve().parent))
import sub6_island_finalize as F
import sub6_iter1_island_retune as R
import cadquery
import pytest
import scipy

RES,PKG=F.RES,F.PKG
validation=dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED",HFSS_solutions_run=False,
    geometry_packages=json.loads((PKG/"geometry_summary.json").read_text(encoding="utf-8")),
    CAD_to_MoM_mapping_all_pass=json.loads((PKG/"cad_to_mom_mapping.json").read_text(encoding="utf-8"))["all_pass"],
    full_fast_suite=dict(passed=2953,skipped=129,failed=8,seconds=945.08,
        accepted_ENV_ONLY_hash_failures=6,poll_timeout_failures=2,
        timeout_rerun_passed=2,timeout_rerun_seconds=41.42,
        exception_authority="user revised ENV_ONLY norm/hash criterion, parent job; runtime timeouts rerun without changing tests"),
    targeted_checks=dict(passed=5,seconds=4.06),src_changes=False,existing_tests_or_fixtures_changed=False,
    optional_mesh192=json.loads((RES/"mesh192_status.json").read_text(encoding="utf-8")))
assert all(validation["geometry_packages"].values()) and validation["CAD_to_MoM_mapping_all_pass"]
(RES/"validation.json").write_text(json.dumps(validation,indent=2),encoding="utf-8",newline="\n")
(PKG/"validation.json").write_text(json.dumps(validation,indent=2),encoding="utf-8",newline="\n")
buf=io.StringIO()
with contextlib.redirect_stdout(buf):
    R.np.show_config()
source_names=["sub6_iter1_island_retune.py","sub6_island_meshtrend.py","sub6_island_cell_packages.py",
    "sub6_island_package_predictions.py","sub6_island_validate_packages.py","sub6_island_finalize.py","sub6_island_manifest.py"]
provenance=dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED",utc_end=datetime.now(timezone.utc).isoformat(),
    machine=platform.node(),OS=platform.platform(),executable=sys.executable,python=sys.version,
    numpy=R.np.__version__,scipy=scipy.__version__,pytest=pytest.__version__,cadquery=cadquery.__version__,
    blas=buf.getvalue(),workers=64,BLAS_threads=1,worker_cap=None,FFT_row_chunk=256,
    git_parent_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
    git_branch=subprocess.check_output(["git","branch","--show-current"],cwd=ROOT,text=True).strip(),
    dirty_flag=bool(subprocess.check_output(["git","status","--porcelain"],cwd=ROOT,text=True)),
    source_sha256={name:F.digest(ROOT/"scripts/research"/name) for name in source_names},
    model="fixed geometry; no tuning; existing FFT sum + exact load update; first-order tensor leakage and two-node sheet beta",
    raw_node_file="nodes.jsonl",expected_completed_nominal_nodes=24,
    actual_completed_nodes=len((RES/"nodes.jsonl").read_text(encoding="utf-8").splitlines()),
    optional_mesh192_status="ALLOCATION_FAILED_NO_NUMERIC_RESULT",HFSS_solutions_run=False,
    inputs_from="initial.log, initial_provenance.json, 192_provenance.json; raw per-node dimensions included")
assert provenance["actual_completed_nodes"]==24
(RES/"provenance.json").write_text(json.dumps(provenance,indent=2),encoding="utf-8",newline="\n")
F.report()
for name,source in [("MESH_TREND.md",RES/"REPORT.md"),("mesh_trend.csv",RES/"mesh_trend.csv"),
                    ("mesh_trend.json",RES/"mesh_trend.json"),("mesh192_status.json",RES/"mesh192_status.json")]:
    (PKG/name).write_bytes(source.read_bytes())
# Normalize this job's text before hashing, matching git eol=lf. STEP files
# are already exactified LF and must retain their recorded SHA256.
for directory in (RES,PKG):
    for path in directory.rglob("*"):
        if not path.is_file() or path.suffix.lower()==".step":
            continue
        data=path.read_bytes().replace(b"\r\n",b"\n")
        if path.suffix in (".log",".err"):
            data=b"\n".join(line.rstrip(b" \t") for line in data.split(b"\n"))
        if data!=path.read_bytes():
            path.write_bytes(data)
for directory in (RES,PKG):
    files=[dict(path=p.relative_to(directory).as_posix(),size_bytes=p.stat().st_size,sha256=F.digest(p))
           for p in sorted(directory.rglob("*")) if p.is_file() and p.name!="manifest.json"]
    (directory/"manifest.json").write_text(json.dumps(dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED",files=files),indent=2),encoding="utf-8",newline="\n")
for folder in PKG.iterdir():
    if folder.is_dir():
        checks=json.loads((folder/"geometry_checks.json").read_text(encoding="utf-8"))
        step=next(folder.glob("*.step"))
        assert F.digest(step)==checks["step_sha256"],"STEP checksum changed"
F.deliver()
print("Finalized: six numerical mesh rows, two mesh192 allocation-failure rows, four checked HFSS packages")
