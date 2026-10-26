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


def test_corr_f_list_option_changes_only_the_correction_frequencies():
    B = asm.B
    old_corr, old_grid = list(B.CORR_F), list(B.GRID_F)
    try:
        pos, kw, dry = B.parse_grid_options(
            [
                "x",
                "corr90",
                "40",
                "--corr-f-list",
                "10.45,10.55",
                "--f-list",
                "10.45,10.55",
                "--theta-step",
                "5",
                "--phi-step",
                "15",
            ]
        )
        assert pos == ["x", "corr90", "40"] and kw["corr_f_list"] == [10.45, 10.55]
        B.configure_grid(**kw)
        assert B.CORR_F == [10.45, 10.55] and B.GRID_F == [10.45, 10.55]
        assert B.CORR_POINTS[:2] == [(0.0, 0.0), (30.0, 0.0)]
    finally:
        B.CORR_F, B.GRID_F = old_corr, old_grid
        B.configure_grid(theta_step=7.5, phi_step=30.0)


def test_missing_corr60_chunk_is_cut_from_the_grid_chunk(tmp_path, monkeypatch):
    pts = [(0.0, 0.0), (5.0, 15.0), (30.0, 0.0), (30.0, 90.0), (45.0, 0.0), (45.0, 90.0)]
    s = np.arange(len(pts) * 16, dtype=float).reshape(len(pts), 4, 2, 2).astype(complex)
    np.savez(tmp_path / "grid_f10.40.npz", f_ghz=10.4, points=np.array(pts), s=s, mesh=60)
    monkeypatch.setitem(asm.CHUNK_DIRS, "grid", str(tmp_path))
    monkeypatch.setitem(asm.CHUNK_DIRS, "corr", str(tmp_path))
    z = asm._corr_chunk("corr60", 10.4)
    assert np.allclose(z["points"], np.array(asm.B.CORR_POINTS))
    idx = [pts.index(tuple(p)) for p in asm.B.CORR_POINTS]
    assert np.array_equal(z["s"], s[idx])


def test_pchip_correction_interpolation_is_shape_preserving_and_cubic_default_unchanged():
    from src.ris.state_library import interpolate_correction

    f_nodes = np.array([10.0, 10.4, 10.45, 10.5, 11.0])
    th = np.array([0.0, 30.0, 45.0])
    step = np.array([0.0, 0.0, 0.0, 40.0, 0.0])  # a narrow peak between uneven nodes
    vals = np.broadcast_to(step[:, None, None, None, None], (5, 3, 1, 1, 1)).astype(complex)
    f_grid = np.linspace(10.0, 11.0, 101)
    cub = interpolate_correction(f_nodes, th, vals, f_grid, np.array([0.0, 15.0]))
    pch = interpolate_correction(f_nodes, th, vals, f_grid, np.array([0.0, 15.0]), "pchip")
    assert np.array_equal(cub, interpolate_correction(f_nodes, th, vals, f_grid, np.array([0.0, 15.0]), "cubic"))
    assert pch.real.max() <= 40.0 + 1e-9 and pch.real.min() >= -1e-9  # no overshoot
    assert np.allclose(interpolate_correction(f_nodes, th, vals, f_nodes, th), vals)  # the default is exact at the nodes too
    assert np.allclose(interpolate_correction(f_nodes, th, vals, f_nodes, th, "pchip"), vals)  # exact at the nodes
    with pytest.raises(ValueError):
        interpolate_correction(f_nodes, th, vals, f_grid, th, "bogus")
