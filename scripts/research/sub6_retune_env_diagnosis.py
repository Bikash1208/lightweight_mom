"""Capture S11 arrays without changing src or tests; diagnose env differences.

Stored fingerprints cannot supply elementwise reference errors. This script
reports norms vs that fixture and preserves raw arrays for an honest future
comparison when the original reference values become available.
"""
import importlib.util
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pytest
import scipy

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location(
    "regression", ROOT/"scripts/research/rotated_load_regression_cases.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
original_hx = module.hx


def capture(a):
    result = original_hx(a)
    arr = np.asarray(a.toarray() if hasattr(a,"toarray") else a,complex)
    if arr.shape == (2,2):
        result["values"] = np.stack([arr.real,arr.imag],axis=-1).tolist()
    return result


module.hx = capture
rows = module.compute()
reference = json.loads((ROOT/"tests/data/rotated_load_main_reference.json").read_text())
for key,row in rows.items():
    if key.endswith("_s11"):
        row["stored_reference_norm"] = reference[key]["norm"]
        row["absolute_norm_difference"] = abs(row["norm"]-reference[key]["norm"])
        row["relative_norm_difference"] = row["absolute_norm_difference"]/reference[key]["norm"]
        row["hash_matches"] = row["sha256"] == reference[key]["sha256"]
        row["max_abs_dS11_vs_stored_reference"] = None
        row["max_relative_difference_vs_stored_reference"] = None
output = dict(machine=platform.node(),executable=sys.executable,python=sys.version,
              numpy=np.__version__,scipy=scipy.__version__,pytest=pytest.__version__,
              threads={k:os.environ.get(k) for k in
                       ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS")},
              stored_reference="SHA256, shape, norm only; complex entries unavailable",
              rows={k:v for k,v in rows.items() if k.endswith("_s11")})
Path(sys.argv[1]).write_text(json.dumps(output,indent=2),encoding="utf-8")
print(json.dumps({"env":output["executable"],"rows":output["rows"]},indent=2))
