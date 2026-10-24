"""Phase 20 R4 item 3: assemble the CP-2 state library from the chunks of `phase20_r4_build_library.py`.

    python scripts/research/phase20_r4_assemble.py [out.npz]      # default results/phase20_r4/cp2_library_v1.npz (+ .manifest.json)
    python scripts/research/phase20_r4_assemble.py out.npz --f-step 0.1 --theta-step 5 --phi-step 15 --grid-chunks <dir> [--corr-chunks <dir>] [--server-manifest MANIFEST.json] [--name vN]

The grid options are those of `phase20_r4_build_library.py` (defaults unchanged: with no option this script behaves exactly as before). `--grid-chunks` / `--corr-chunks` give the
directories of the `grid_f*.npz` and `corr*_f*.npz` chunks (default: the build script's chunk directory); `--server-manifest` verifies every chunk's sha256 against a server MANIFEST.json before use.

Reads the grid chunks (mesh 60) and the correction chunks (mesh 60 / 90 / 120), fills theta = 0 by the exact frame rotation, refers everything to the fixed structure frame,
builds the (f, theta) mesh-correction map and the uncertainty (src/ris/state_library.py), writes one compressed npz plus a JSON manifest. MODEL_ONLY.
"""

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import phase20_r4_build_library as B

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.ris.polarization import circular_reflection, to_structure_frame  # noqa: E402
from src.ris.state_library import THETA_LABEL_LIMIT_DEG, assemble, fill_normal_incidence, save_library  # noqa: E402

CHUNK_DIRS = {}  # kind -> directory override ("grid", "corr")


def chunk_path(kind, f):
    d = CHUNK_DIRS.get("grid" if kind == "grid" else "corr")
    return (Path(d) / f"{kind}_f{f:.2f}.npz") if d else B.chunk_path(kind, f)


def verify_manifest(manifest_path, directories):
    """sha256 of every chunk named in a server MANIFEST.json (paths with backslashes) against the file found in `directories`; raises on any mismatch or missing file."""
    man = json.loads(Path(manifest_path).read_text())
    checked = 0
    for ent in man["files"]:
        name = ent["file"].replace("\\", "/").split("/")[-1]
        found = next((Path(d) / name for d in directories if (Path(d) / name).exists()), None)
        if found is None:
            if name.endswith(".npz"):
                raise FileNotFoundError(name)
            continue  # logs / provenance files live next to the chunks directory, not in it
        if hashlib.sha256(found.read_bytes()).hexdigest() != ent["sha256"]:
            raise ValueError(f"sha256 mismatch: {found}")
        checked += 1
    return checked


def load_grid():
    pts = B.grid_points()
    s = np.zeros((len(B.GRID_F), len(B.GRID_THETA), len(B.GRID_PHI), 4, 2, 2), dtype=complex)
    for i, f in enumerate(B.GRID_F):
        z = np.load(chunk_path("grid", f))
        assert np.allclose(z["points"], np.array(pts))
        for (th, ph), val in zip(pts, z["s"], strict=True):
            s[i, B.GRID_THETA.index(th), B.GRID_PHI.index(ph)] = val
    return fill_normal_incidence(s, B.GRID_THETA, B.GRID_PHI)


def _corr_chunk(kind, f):
    """The correction chunk; a missing mesh-60 one is cut out of the grid chunk of that frequency (the correction points are grid points, same solves)."""
    path = chunk_path(kind, f)
    if kind == "corr60" and not path.exists():
        z = np.load(chunk_path("grid", f))
        idx = [next(k for k, pt in enumerate(z["points"]) if np.allclose(pt, p)) for p in B.CORR_POINTS]
        return {"points": np.array(B.CORR_POINTS), "s": z["s"][idx]}
    return np.load(path)


def load_corr(kind):
    CORR_F, CORR_POINTS = B.CORR_F, B.CORR_POINTS
    out = np.zeros((len(CORR_F), len(CORR_POINTS), 4, 2, 2), dtype=complex)
    for i, f in enumerate(CORR_F):
        z = _corr_chunk(kind, f)
        assert np.allclose(z["points"], np.array(CORR_POINTS))
        for j, ((_, ph), sk) in enumerate(zip(CORR_POINTS, z["s"], strict=True)):
            for k in range(4):
                out[i, j, k] = to_structure_frame(circular_reflection(sk[k]), np.radians(ph))
    return out


def main(out, name="ris_cp2_octagon_madp_rt5880_15", version=1, note=None):
    GRID_F, GRID_THETA, GRID_PHI = B.GRID_F, B.GRID_THETA, B.GRID_PHI
    CORR_F, CORR_POINTS = B.CORR_F, B.CORR_POINTS
    s_raw = load_grid()
    c60, c90, c120 = (load_corr(k) for k in ("corr60", "corr90", "corr120"))
    arrays = assemble(GRID_F, GRID_THETA, GRID_PHI, s_raw, CORR_F, CORR_POINTS, c60, c90, c120)
    dc = arrays["correction"][:, :, 0]  # (f, theta, state, 2, 2)
    rr = arrays["c_raw"][:, :, 0, :, 0, 0]
    phase_corr = np.degrees(np.angle((rr + dc[..., 0, 0]) / rr))
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=ROOT).stdout.strip()
    manifest = {
        "name": name,
        "version": version,
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "PROVISIONAL_DATASHEET", "VIA_OMITTED_APPROXIMATION"],
        "axes": {"f_ghz": GRID_F, "theta_deg": GRID_THETA, "phi_deg": GRID_PHI, "state": [1, 2, 3, 4]},
        "theta_label": {
            f"{t:g}": ("EXTRAPOLATED_BEYOND_R3_TARGETS" if t > THETA_LABEL_LIMIT_DEG else "within_R3_targets")
            for t in GRID_THETA
        },
        "frame": "circular reflection C[reflected, incident] (R, L; propagation-relative) referred to the fixed structure axes; s_corr in the incidence-plane (TE, TM) frame",
        "grid_mesh": 60,
        "floquet_cutoff": "2 n (n = cells per period) at every mesh",
        "grid_reduction_note": note
        or "reduced from the approved 41 x 10 x 24 grid after the Floquet-cutoff finding (results/phase20_r4/FLOQUET_CONVERGENCE.md)",
        "correction": {
            "meshes": [60, 90, 120],
            "f_ghz": CORR_F,
            "points": CORR_POINTS,
            "definition": "dC = C_inf(p=2) - C(60) with C_inf = C(120) + R2 (C(120) - C(90)), phi-mean, interpolated in (f, theta), added to every phi; uncertainty = |C_inf(p=1) - C_inf(p=2)| + phi spread; all meshes with Floquet cutoff F = 2n",
            "r_RR_phase_correction_deg": {
                "min": float(phase_corr.min()),
                "max": float(phase_corr.max()),
                "median": float(np.median(phase_corr)),
            },
            "max_uncertainty_abs": float(arrays["uncertainty"].max()),
        },
        "cell": {
            "geometry_json": "results/phase20_r3/opt/cp2_madp_final.json",
            "diode": "MACOM MADP-000907-14020 (L_s = 0 assumed)",
            "period_mm": 15.0,
        },
        "built_at_commit": commit,
    }
    save_library(out, arrays, manifest)
    manifest["sha256"] = hashlib.sha256(Path(out).read_bytes()).hexdigest()
    Path(out).with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=1))
    print("wrote", out, f"{Path(out).stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    argv, grid_kw, dry = B.parse_grid_options(sys.argv)
    B.configure_grid(**grid_kw)
    opts = {}
    for key in ("--grid-chunks", "--corr-chunks", "--server-manifest", "--name"):
        if key in argv:
            i = argv.index(key)
            opts[key] = argv[i + 1]
            del argv[i : i + 2]
    if "--grid-chunks" in opts:
        CHUNK_DIRS["grid"] = opts["--grid-chunks"]
    if "--corr-chunks" in opts:
        CHUNK_DIRS["corr"] = opts["--corr-chunks"]
    if "--server-manifest" in opts:
        print(
            "sha256 verified:",
            verify_manifest(
                opts["--server-manifest"], [d for d in (opts.get("--grid-chunks"), opts.get("--corr-chunks")) if d]
            ),
            "files",
        )
    full = bool(grid_kw)
    main(
        argv[1] if len(argv) > 1 else ROOT / "results/phase20_r4/cp2_library_v1.npz",
        name=opts.get("--name", "ris_cp2_octagon_madp_rt5880_15"),
        version=2 if full else 1,
        note="full approved grid 41 x 10 x 24 (server run r4-full-grid, mesh 60, F = 2n; corrections from the R4 mesh 60/90/120 set, mesh 90 re-run on the server agrees to 1e-13)"
        if full
        else None,
    )
