"""Cross-machine equivalence tools (infra brief rule 5): comparison arithmetic and the provenance record."""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))
import equivalence_compare as ec  # noqa: E402
import provenance  # noqa: E402


def _write(path, s11, machine):
    cases = [{"state": 1, "theta_deg": 0.0, "phi_deg": 0.0, "s11": [[[z.real, z.imag] for z in row] for row in s11]}]
    path.write_text(json.dumps({"provenance": {"machine": machine}, "cases": cases}))


def test_identical_and_perturbed_cases(tmp_path):
    s = np.array([[0.9 * np.exp(1j * 0.3), 0.01 + 0.02j], [0.01 + 0.02j, 0.8 * np.exp(-1j * 1.1)]])
    _write(tmp_path / "a.json", s, "A")
    _write(tmp_path / "b.json", s, "B")
    out = ec.main(tmp_path / "a.json", tmp_path / "b.json")
    assert out["verdict"] == "INTERCHANGEABLE" and out["max_rel_complex_diff"] == 0.0
    s2 = s.copy()
    s2[0, 0] *= np.exp(1j * np.radians(0.5)) * 10 ** (0.1 / 20)  # 0.5 deg, 0.1 dB
    _write(tmp_path / "c.json", s2, "C")
    out = ec.main(tmp_path / "a.json", tmp_path / "c.json")
    assert abs(out["max_abs_dphase_deg"] - 0.5) < 1e-9 and abs(out["max_abs_dlevel_db"] - 0.1) < 1e-9
    assert out["verdict"].startswith("DIFFERENCE_ABOVE_TOLERANCE")


def test_provenance_has_the_required_fields():
    rec = provenance.finish(provenance.provenance({"x": 1}, workers=3))
    for k in (
        "machine",
        "os",
        "git_commit",
        "git_dirty",
        "python",
        "numpy",
        "blas",
        "threads",
        "workers",
        "started",
        "finished",
        "x",
    ):
        assert k in rec
    assert rec["workers"] == 3
