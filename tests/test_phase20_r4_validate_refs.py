"""The R4 validation script can re-evaluate a library on SAVED direct-solve references (no solves): the server's v2 draws are reproduced from their reference files."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
R4 = ROOT / "results" / "phase20_r4"


@pytest.mark.parametrize("name", ["v2_a", "v2_a_phi90", "v1_a_phi90"])
def test_use_refs_reproduces_the_recorded_summary(name, tmp_path):
    lib = "cp2_library_v1.npz" if name.startswith("v1") else "cp2_library_v2_fullgrid.npz"
    out = tmp_path / "o.json"
    env = {**os.environ, "PHASE20_R4_LIB": str(R4 / lib)}
    subprocess.run([sys.executable, str(ROOT / "scripts/research/phase20_r4_validate.py"), "a", "0", str(out), "--use-refs", str(R4 / f"validation_{name}_refs.npz")], check=True, env=env, capture_output=True)
    new, old = json.loads(out.read_text()), json.loads((R4 / f"validation_{name}.json").read_text())
    assert new["gate_theta_le_30"]["result"] == old["gate_theta_le_30"]["result"]
    assert new["gate_theta_le_30"]["p95_phase_deg"] == pytest.approx(old["gate_theta_le_30"]["p95_phase_deg"], abs=1e-6)
    assert "in_band_10p4_10p9.theta_le_30.phase_err_deg" in new["summary"]
