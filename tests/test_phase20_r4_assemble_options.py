"""The R4 assemble script's grid / chunk-directory options: defaults unchanged, manifest check refuses a modified chunk."""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import phase20_r4_assemble as asm  # noqa: E402


def test_default_chunk_path_is_the_build_scripts_and_override_works(tmp_path, monkeypatch):
    assert asm.chunk_path("grid", 8.0) == asm.B.chunk_path("grid", 8.0)
    monkeypatch.setitem(asm.CHUNK_DIRS, "grid", str(tmp_path))
    monkeypatch.setitem(asm.CHUNK_DIRS, "corr", str(tmp_path / "c"))
    assert asm.chunk_path("grid", 8.0) == tmp_path / "grid_f8.00.npz"
    assert asm.chunk_path("corr90", 8.0) == tmp_path / "c" / "corr90_f8.00.npz"


def test_manifest_sha256_check(tmp_path):
    f = tmp_path / "grid_f8.00.npz"
    np.savez(f, a=np.arange(3))
    man = tmp_path / "MANIFEST.json"
    man.write_text(
        json.dumps(
            {
                "files": [
                    {"file": "chunks\\grid_f8.00.npz", "sha256": hashlib.sha256(f.read_bytes()).hexdigest()},
                    {"file": "run.log", "sha256": "x"},
                ]
            }
        )
    )
    assert asm.verify_manifest(man, [tmp_path]) == 1
    np.savez(f, a=np.arange(4))
    with pytest.raises(ValueError):
        asm.verify_manifest(man, [tmp_path])
