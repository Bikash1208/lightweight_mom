"""Phase 30 cell LIBRARY (MODEL_ONLY, NOT_HFSS_VALIDATED): E1 family at eps_r 11.2, 10.05 GHz -- RO3010 0.635 mm, 3.2 mm lattice, full-fill super-ellipse n = 4
(major = minor = 3.05 mm, gap 0.15 mm), centre slot L x W, OUTLINE FIXED, slot turns.  Slot sizes are multiples of two mesh cells (mesh 48: 0.0667 mm) so the
rasterised slot is symmetric and not quantisation-noisy (odd cell counts collapse onto the even one below).  Three products:

  uniform  : per-slot uniform-cell data at angle 0 (beta/k0, X_par/X_perp, M_t, first-order alpha, Floquet alpha) -> `uniform_table.json` (feeds the ray-march design)
  lut      : full 2 x 2 Bloch sheet admittance Y(slot, a, psi) at kt = 1.30 k0, a, psi = 0..165 deg step 15 -> `angle_lut_v1.npz/json`
  heldout  : N random off-grid (L, a, psi) points solved by MoM and compared with the interpolated LUT -> `heldout.json`

    P30_EPS=11.2 P30_F0=10.05e9 python scripts/research/phase30_library.py {uniform|lut|heldout} results/phase30/library [--workers 44]
Heavy: run through jobs/phase30-library.md.  Single-threaded BLAS per worker."""

import json
import os
import sys
import time
from multiprocessing import get_context
from pathlib import Path

os.environ.setdefault("P30_EPS", "11.2")
os.environ.setdefault("P30_F0", "10.05e9")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

import phase30_refine as R  # noqa: E402

MESH = 48
HC = R.PITCH / MESH  # mesh-cell size
MAJOR = R.PITCH - R.GAP
KT = 1.30
ANG = np.arange(0.0, 180.0, 15.0)
NL = [0, 12, 18, 22, 26, 30, 34]  # slot length in mesh cells (L = n h: 0 .. 2.27 mm; bridge >= 0.39 mm)
NW = [2, 4, 6, 8, 10]  # slot width in mesh cells (0.13 .. 0.67 mm)
HELD_NL = [14, 16, 20, 24, 28, 32]


def slots():
    out = [(0, 0)]
    out += [(ln, w) for ln in NL[1:] for w in NW]
    return out


def dims(nl, nw):
    return nl * HC, nw * HC


def _spec(nl, nw, angle=0.0):
    ln, w = dims(nl, nw)
    return dict(major=MAJOR, ratio=1.0, slot_l=ln, slot_w=w, angle=angle)


def _cell(nl, nw, slot_angle_deg):
    from src.geometry import UnitCell
    from src.materials import Metal
    from src.periodic_bc import RectangularLattice

    ln, w = dims(nl, nw)
    geo = R.EllipticalPatch(MAJOR, MAJOR, 0.0, ln, w, float(np.radians(slot_angle_deg)), False, R.EXPO)
    return UnitCell(RectangularLattice(R.PITCH, R.PITCH), geo, R.SUB, Metal(), backing="pec")


def _uniform_job(a):
    nl, nw = a
    r = R.evaluate(_spec(nl, nw), R.F0, MESH, with_loss=False, floquet=True)
    r["nl"], r["nw"] = nl, nw
    return r


def _tensor_row(cell, angles_psi):
    from src.basis_functions import MeshSettings, RooftopBasis
    from src.dispersion import tensor_sheet as ts
    from src.frequency import free_space_wavenumber

    mesh = MeshSettings(MESH, MESH, MESH)
    basis = RooftopBasis(cell, mesh.cells_x, mesh.cells_y)
    kt = KT * free_space_wavenumber(R.F0)
    return np.array([ts.bloch_sheet_tensor(cell, R.F0, mesh, kt * np.cos(p), kt * np.sin(p), basis=basis) for p in np.radians(angles_psi)])


def _lut_job(a):
    nl, nw, ia = a
    return nl, nw, ia, _tensor_row(_cell(nl, nw, ANG[ia]), ANG)


def _held_job(a):
    nl, nw, adeg, pdeg = a
    return a, _tensor_row(_cell(nl, nw, adeg), [pdeg])[0]


def _pool(n):
    return get_context("spawn").Pool(n)


def stage_uniform(out, workers):
    t0 = time.time()
    with _pool(workers) as p:
        rows = p.map(_uniform_job, slots(), chunksize=1)
    (out / "uniform_table.json").write_text(
        json.dumps({"rows": rows, "eps_r": R.EPS_R, "f_hz": R.F0, "mesh": MESH, "seconds": time.time() - t0, "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"]}, indent=1, default=float)
    )
    print("uniform", len(rows), "rows,", sum(1 for r in rows if r["ok"]), "ok,", round(time.time() - t0), "s", flush=True)


def stage_lut(out, workers):
    t0 = time.time()
    sl = slots()
    jobs = [(ln, w, ia) for (ln, w) in sl for ia in range(len(ANG))]
    Y = np.zeros((len(sl), len(ANG), len(ANG), 2, 2), dtype=complex)
    idx = {s: i for i, s in enumerate(sl)}
    with _pool(workers) as p:
        for nl, nw, ia, y in p.imap_unordered(_lut_job, jobs, chunksize=1):
            Y[idx[(nl, nw)], ia] = y
    np.savez_compressed(out / "angle_lut_v1.npz", Y=Y, slots_cells=np.array(sl), a_deg=ANG, psi_deg=ANG)
    R90 = np.array([[0.0, -1.0], [1.0, 0.0]])
    e = [
        np.linalg.norm(Y[s, (ia + 6) % 12, (ip + 6) % 12] - R90 @ Y[s, ia, ip] @ R90.T) / max(np.linalg.norm(Y[s, ia, ip]), 1e-30)
        for s in range(len(sl))
        for ia in range(12)
        for ip in range(12)
    ]
    asym = np.abs(Y[..., 0, 1] - Y[..., 1, 0]) / np.maximum(np.abs(Y[..., 0, 1]), 1e-30)
    meta = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "BLOCH_SHEET_TENSOR_FROM_PERIODIC_MOM"],
        "eps_r": R.EPS_R,
        "f_hz": R.F0,
        "kt_over_k0": KT,
        "mesh": MESH,
        "mesh_cell_mm": HC * 1e3,
        "slots_cells_nl_nw": sl,
        "slots_mm": [list(np.round(np.array(dims(*s)) * 1e3, 4)) for s in sl],
        "a_deg": ANG.tolist(),
        "psi_deg": ANG.tolist(),
        "axes": "Y[slot, a_i, psi_j, 2, 2] siemens; slot turned by a (outline fixed), kt along psi",
        "c4_relative_error_median": float(np.median(e)),
        "c4_relative_error_max": float(np.max(e)),
        "reciprocity_max_asym": float(np.nanmax(asym)),
        "seconds": time.time() - t0,
    }
    (out / "angle_lut_v1.json").write_text(json.dumps(meta, indent=1))
    print("lut done", round(time.time() - t0), "s", meta["c4_relative_error_median"], meta["c4_relative_error_max"], flush=True)


def trig_interp(v, a_deg, axis):
    """Periodic (period 180 deg) trigonometric interpolation of `v` along `axis` (12 samples) at the single angle `a_deg`."""
    n = v.shape[axis]
    c = np.fft.fft(v, axis=axis)
    k = np.fft.fftfreq(n, d=1.0 / n)
    ph = np.exp(2j * np.pi * (a_deg / 180.0) * k)
    ph = np.where(np.abs(k) == n / 2, np.cos(2.0 * np.pi * (a_deg / 180.0) * k), ph)  # Nyquist term: cosine (real, symmetric)
    shape = [1] * v.ndim
    shape[axis] = n
    return np.sum(c * ph.reshape(shape), axis=axis) / n


def lut_eval(Y, sl, nl, nw, a_deg, psi_deg):
    """Interpolate Y at (nl [off-grid], nw [grid], a, psi): trig in a then psi, PCHIP in L (through L = 0), linear in W."""
    from scipy.interpolate import PchipInterpolator

    yp = trig_interp(trig_interp(Y, a_deg, 1), psi_deg, 1)  # [slot, 2, 2]
    d = {s: yp[i] for i, s in enumerate(sl)}
    nls = sorted({s[0] for s in sl if s[0] > 0})
    nws = sorted({s[1] for s in sl if s[1] > 0})
    x = np.array([0] + nls, float)

    def at_w(w):
        v = np.array([d[(0, 0)]] + [d[(ln, w)] for ln in nls])
        return np.array([[PchipInterpolator(x, v[:, i, j].real)(nl) + 1j * PchipInterpolator(x, v[:, i, j].imag)(nl) for j in range(2)] for i in range(2)])

    per_w = {w: at_w(w) for w in nws}
    out = np.zeros((2, 2), dtype=complex)
    for i in range(2):
        for j in range(2):
            out[i, j] = np.interp(nw, nws, [per_w[w][i, j].real for w in nws]) + 1j * np.interp(nw, nws, [per_w[w][i, j].imag for w in nws])
    return out


def stage_heldout(out, workers, n=60, seed=7):
    rng = np.random.default_rng(seed)
    t0 = time.time()
    pts = [(int(rng.choice(HELD_NL)), int(rng.choice(NW)), float(rng.uniform(0, 180)), float(rng.uniform(0, 180))) for _ in range(n)]
    with _pool(workers) as p:
        res = p.map(_held_job, pts, chunksize=1)
    d = np.load(out / "angle_lut_v1.npz")
    Y, sl = d["Y"], [tuple(map(int, s)) for s in d["slots_cells"]]
    rows = []
    for (nl, nw, a, p_), y in res:
        yl = lut_eval(Y, sl, nl, nw, a, p_)
        rel = float(np.linalg.norm(yl - y) / np.linalg.norm(y))
        rows.append(
            {
                "nl": nl,
                "nw": nw,
                "L_mm": nl * HC * 1e3,
                "W_mm": nw * HC * 1e3,
                "a_deg": a,
                "psi_deg": p_,
                "rel_err_Y": rel,
                "Y_mom_mS": [[complex(z).imag * 1e3 for z in r] for r in y],
                "Y_lut_mS": [[complex(z).imag * 1e3 for z in r] for r in yl],
            }
        )
    e = np.array([r["rel_err_Y"] for r in rows])
    (out / "heldout.json").write_text(
        json.dumps(
            {
                "n": n,
                "seed": seed,
                "rel_err_median": float(np.median(e)),
                "rel_err_p90": float(np.percentile(e, 90)),
                "rel_err_max": float(e.max()),
                "frac_below_5pct": float(np.mean(e < 0.05)),
                "rows": rows,
                "seconds": time.time() - t0,
            },
            indent=1,
        )
    )
    print("heldout", n, "median", np.median(e), "p90", np.percentile(e, 90), "max", e.max(), flush=True)


if __name__ == "__main__":
    st, outd = sys.argv[1], Path(sys.argv[2])
    outd.mkdir(parents=True, exist_ok=True)
    nw_ = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 44
    {"uniform": stage_uniform, "lut": stage_lut, "heldout": stage_heldout}[st](outd, nw_)
