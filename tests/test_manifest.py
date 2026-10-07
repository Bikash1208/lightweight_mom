"""Phase 13 Part 2: the run manifest and diagnostic bundle."""

import json
import zipfile

from src import manifest as mf


def test_run_manifest_fills_provenance_fields_automatically():
    m = mf.RunManifest(run_id="abc123", job_kind="unit_cell_solve", inputs={"length_mm": 2.0})
    assert m.python_version
    assert m.platform
    assert isinstance(m.dependency_versions, dict) and "numpy" in m.dependency_versions
    assert m.dependency_versions["numpy"] is not None  # numpy is always installed here


def test_git_hash_is_a_real_commit_in_this_checkout():
    h = mf.git_hash()
    assert h is None or (isinstance(h, str) and len(h) == 40)


def test_add_output_file_hashes_it(tmp_path):
    f = tmp_path / "result.json"
    f.write_text('{"a": 1}')
    m = mf.RunManifest(run_id="r1", job_kind="test")
    m.add_output_file("result.json", f)
    assert len(m.output_files["result.json"]) == 64  # sha256 hex digest length
    # deterministic: same content -> same hash
    f2 = tmp_path / "result2.json"
    f2.write_text('{"a": 1}')
    m.add_output_file("result2.json", f2)
    assert m.output_files["result.json"] == m.output_files["result2.json"]


def test_write_and_reload_roundtrip(tmp_path):
    m = mf.RunManifest(run_id="r2", job_kind="reflectarray_design", labels=("EXPERIMENTAL",))
    path = m.write(tmp_path / "manifest.json")
    reloaded = json.loads(path.read_text())
    assert reloaded["run_id"] == "r2"
    assert reloaded["labels"] == ["EXPERIMENTAL"]


def test_build_diagnostic_bundle_contains_manifest_and_inputs(tmp_path):
    log_file = tmp_path / "lwmom.log"
    log_file.write_text('{"run_id": "r3", "message": "solve started"}\n{"run_id": "other", "message": "noise"}\n')
    input_file = tmp_path / "input.yaml"
    input_file.write_text("period_x: 8mm\n")
    m = mf.RunManifest(run_id="r3", job_kind="unit_cell_solve")

    bundle_path = mf.build_diagnostic_bundle(m, log_file, [input_file], tmp_path / "bundles")
    assert bundle_path.exists()
    with zipfile.ZipFile(bundle_path) as zf:
        names = zf.namelist()
        assert "manifest.json" in names
        assert "log_excerpt.jsonl" in names
        assert "inputs/input.yaml" in names
        excerpt = zf.read("log_excerpt.jsonl").decode()
        assert "solve started" in excerpt
        assert "noise" not in excerpt
