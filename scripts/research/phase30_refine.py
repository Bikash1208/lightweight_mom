"""Phase 30 item 2b (MODEL_ONLY): refined screen of the approved family -- RO3010 0.635 mm, 3.2 mm pitch, full-fill
super-ellipse (n = 4) with a centre slot (long side along the major axis x), gap fixed at 0.15 mm along x.

    python scripts/research/phase30_refine.py grid  <out_dir> [--mesh 48] [--workers 44]
    python scripts/research/phase30_refine.py stage <out_dir> --picks picks.json [--workers 20]

grid : minor/major ratio x slot length x slot width (W = 0: no slot) at 10.54 GHz, per cell: tensor (Y_par, Y_perp), X, M_t,
       anisotropy, beta/k0, first-order alpha/k0, dielectric and conductor loss (Cu, sigma 5.8e7 S/m, smooth, 18 um thick
       >> skin depth: ground plane + patch).
stage: for picked cells, 10.3 / 10.54 / 10.8 GHz (beta, alpha, first-order and Floquet), mesh convergence (48/64/80),
       gap +-0.025 mm, rotation check (22.5/45/90 deg).

Conductor loss (ESTIMATE, a lower bound: rooftop currents at mesh resolution under-resolve edge crowding; surface
roughness ignored): with the mode normalised to unit tangential field E0 = 1 V/m at the patch plane,
  patch  : p = 0.5 Rs * integral |J|^2 dA / (Px Py), J from the Bloch MoM solve for E_tot = x_hat (the same solve as Y_par);
  ground : p = 0.5 Rs |H_g|^2 with |H_g| = |Y_slab(beta)| / |cos(kz h)| (TM slab field H ~ cos(kz (z + h)));
  alpha_c = p / (2 P_sw), P_sw = w k0^2 / (eta0 kappa)^2 per |E0|^2 (leakage_model.surface_wave_power_per_h2).
"""

import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.basis_functions import MeshSettings, RooftopBasis  # noqa: E402
from src.constants import ETA0, MU0  # noqa: E402
from src.dispersion import tensor_sheet as ts  # noqa: E402
from src.dispersion.elliptical_cells import EllipticalPatch  # noqa: E402
from src.dispersion.leakage_model import leakage_alpha_first_order_tensor, surface_wave_power_per_h2  # noqa: E402
from src.dispersion.sheet_model import _y_slab, find_beta_from_sheet_admittance  # noqa: E402
from src.dispersion.tensor_floquet import floquet_complex_beta_tensor  # noqa: E402
from src.frequency import angular_frequency, free_space_wavenumber  # noqa: E402
from src.geometry import UnitCell  # noqa: E402
from src.green_function import LayeredGreenFunction  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402
from src.solver import solve  # noqa: E402

PITCH = 3.2e-3
EPS_R = float(os.environ.get("P30_EPS", 10.2))  # RO3010 process Dk 10.2 +-0.30 @10 GHz (z axis); 'design Dk' 11.2 (sensitivity run)
TAND, H = 0.0022, float(os.environ.get("P30_H", 0.635e-3))  # Rogers RO3010 typical DF @10 GHz; 25 mil

SIGMA_CU = 5.8e7
SUB = Substrate("ro3010", EPS_R, TAND, H)
THETA = np.radians(17.0)
EXPO = 4.0
GAP = 0.15e-3
F0 = float(os.environ.get("P30_F0", 10.54e9))  # default unchanged; the eps_r 11.2 library sets 10.05e9


def make_cell(major, ratio, slot_l, slot_w, angle=0.0):
    geo = EllipticalPatch(major, major * ratio, angle, slot_l, slot_w, 0.0, False, EXPO)
    return UnitCell(RectangularLattice(PITCH, PITCH), geo, SUB, Metal(), backing="pec")


def conductor_alpha(cell, f, mesh, beta, basis):
    """alpha_c [Np/m] (ESTIMATE, see module docstring) for the x-directed Bloch wave at |kt| = beta."""
    omega = angular_frequency(f)
    k0 = free_space_wavenumber(f)
    green = LayeredGreenFunction(cell.substrate, f, backing="pec")
    z = assemble_impedance_matrix(basis, green, beta, 0.0, mesh.floquet_order)
    v = np.stack([assemble_excitation_vector(basis, e, beta, 0.0) for e in np.eye(2)], axis=1)
    cur = solve(z, v)
    j = np.stack([basis.current_coefficients(cur[:, k], beta, 0.0) for k in range(2)], axis=1)
    g00 = np.asarray(green.spectral_dyad(beta, 0.0), dtype=complex)
    e_tot = np.eye(2) + g00 @ j
    coeff = cur @ np.linalg.inv(e_tot)[:, 0]  # currents for E_tot = x_hat (unit field at the patch plane)
    n = 96
    g = (np.arange(n) + 0.5) / n * PITCH - PITCH / 2
    xx, yy = np.meshgrid(g, g, indexing="ij")
    jxy = basis.evaluate(xx, yy, coeff, beta, 0.0)
    rs = np.sqrt(omega * MU0 / (2.0 * SIGMA_CU))
    p_patch = 0.5 * rs * np.sum(np.abs(jxy) ** 2) * (PITCH / n) ** 2 / PITCH**2
    y_slab = _y_slab(beta, omega, cell.substrate.eps_complex, k0, H)
    kz = np.sqrt(cell.substrate.eps_complex * k0**2 - beta**2)
    p_gnd = 0.5 * rs * (abs(y_slab) / abs(np.cos(kz * H))) ** 2
    kappa = np.sqrt(beta**2 - k0**2)
    w = surface_wave_power_per_h2(beta, EPS_R, H, f)
    p_sw = w * k0**2 / (ETA0 * kappa) ** 2  # per |E0|^2 = 1
    return float((p_patch + p_gnd) / (2.0 * p_sw)), float(p_patch / (2 * p_sw)), float(p_gnd / (2 * p_sw))


def dielectric_alpha(y_bar, beta, f):
    omega = angular_frequency(f)
    k0 = free_space_wavenumber(f)
    from src.dispersion.sheet_model import _y_air

    def g(b):
        return 1j * y_bar + _y_air(b, omega, k0) + _y_slab(b, omega, SUB.eps_complex, k0, H)

    b0, b1 = complex(beta), complex(beta) * (1 + 1e-6)
    g0, g1 = g(b0), g(b1)
    for _ in range(30):
        if g1 == g0:
            break
        b0, b1, g0 = b1, b1 - g1 * (b1 - b0) / (g1 - g0), g1
        g1 = g(b1)
        if abs(g1) < 1e-14:
            break
    return abs(b1.imag)


def evaluate(spec, f=F0, mesh_n=48, with_loss=True, floquet=False):
    """spec: dict(major, ratio, slot_l, slot_w, angle). Returns a result dict (never raises)."""
    k0 = free_space_wavenumber(f)
    out = dict(spec)
    out.update({"f_hz": f, "mesh": mesh_n, "ok": False})
    try:
        mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
        cell = make_cell(spec["major"], spec["ratio"], spec["slot_l"], spec["slot_w"], spec.get("angle", 0.0))
        basis = RooftopBasis(cell, mesh.cells_x, mesh.cells_y)
        yn = ts.sheet_admittance_tensor_for_cell(cell, f, mesh)
        beta0 = find_beta_from_sheet_admittance(1j * 0.5 * (yn[0, 0].imag + yn[1, 1].imag), SUB, f).beta_rad_per_m
        yp = ts.bloch_sheet_tensor(cell, f, mesh, beta0, 0.0, basis=basis)
        yq = ts.bloch_sheet_tensor(cell, f, mesh, 0.0, beta0, basis=basis)
        y_par, y_perp = yp[0, 0].imag, yq[1, 1].imag
        y_bar, delta = 0.5 * (y_par + y_perp), 0.5 * abs(y_par - y_perp)
        mode = find_beta_from_sheet_admittance(1j * y_bar, SUB, f)
        beta = mode.beta_rad_per_m
        omega = angular_frequency(f)
        y_slab = _y_slab(beta, omega, SUB.eps_complex, k0, H)
        x_par = float(np.imag(1.0 / (1j * y_par + y_slab)))
        x_perp = float(np.imag(1.0 / (1j * y_perp + y_slab)))
        alpha = leakage_alpha_first_order_tensor(beta, 1j * y_bar, delta, THETA, SUB, f)
        out.update({
            "ok": True, "y_par_mS": y_par * 1e3, "y_perp_mS": y_perp * 1e3, "y_bar_mS": y_bar * 1e3, "delta_mS": delta * 1e3,
            "beta_over_k0": beta / k0, "x_bar_ohm": mode.x_s_ohm, "x_par_ohm": x_par, "x_perp_ohm": x_perp,
            "M_t": abs(x_par - x_perp) / abs(x_par + x_perp),
            "aniso_ratio": max(abs(x_par), abs(x_perp)) / max(min(abs(x_par), abs(x_perp)), 1e-9),
            "alpha_np_m": alpha, "alpha_over_k0": alpha / k0,
        })
        if with_loss:
            a_d = dielectric_alpha(y_bar, beta, f)
            a_c, a_cp, a_cg = conductor_alpha(cell, f, mesh, beta0, basis)
            out.update({"alpha_diel_np_m": a_d, "alpha_cond_np_m": a_c, "alpha_cond_patch_np_m": a_cp, "alpha_cond_gnd_np_m": a_cg,
                        "rad_over_diel": alpha / a_d if a_d else None, "rad_over_cond": alpha / a_c if a_c else None})
        if floquet:
            n_cell = max(int(round(2.0 * np.pi / beta / PITCH)), 4)
            tens = []
            for i in range(n_cell):
                a = np.pi * i / n_cell
                c, s = np.cos(a), np.sin(a)
                r = np.array([[c, -s], [s, c]])
                tens.append(1j * (y_bar * np.eye(2) + delta * r @ np.diag([1.0, -1.0]) @ r.T))
            try:
                bf, _ = floquet_complex_beta_tensor(np.array(tens), n_cell * PITCH, SUB, f, beta, n_max=4)
                out.update({"N_cells_period": n_cell, "alpha_floquet_over_k0": float(abs(bf.imag)) / k0,
                            "beta_floquet_over_k0": float(bf.real) / k0})
                out["floquet_ratio"] = out["alpha_floquet_over_k0"] / out["alpha_over_k0"] if out["alpha_over_k0"] else None
            except Exception as e:  # noqa: BLE001
                out["floquet_error"] = f"{type(e).__name__}: {e}"
    except Exception as e:  # noqa: BLE001 -- record, never fabricate
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def grid_specs():
    major = PITCH - GAP
    specs = []
    for ratio in (1.0, 0.97, 0.94):
        specs.append(dict(major=major, ratio=ratio, slot_l=0.0, slot_w=0.0, angle=0.0))
        for lf in (0.40, 0.50, 0.62, 0.75):
            for w in (0.10e-3, 0.20e-3, 0.30e-3, 0.45e-3, 0.58e-3, 0.72e-3, 0.90e-3):
                specs.append(dict(major=major, ratio=ratio, slot_l=lf * major, slot_w=w, angle=0.0))
    return specs


def _grid_job(a):
    return evaluate(a[0], F0, a[1])


def run_grid(out_dir, mesh_n=48, workers=44):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    specs = grid_specs()
    t0 = time.time()
    rows = []
    with Pool(workers) as pool:
        for i, r in enumerate(pool.imap(_grid_job, [(s, mesh_n) for s in specs], chunksize=1)):
            rows.append(r)
            with open(out / "refine_grid.jsonl", "a") as fh:
                fh.write(json.dumps(r) + "\n")
            if i % 10 == 0:
                print(f"{i + 1}/{len(specs)} [{time.time() - t0:.0f}s]", flush=True)
    (out / "refine_grid.json").write_text(json.dumps({"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"], "mesh": mesh_n, "rows": rows}, indent=1))
    print("done", len(rows), sum(r["ok"] for r in rows), f"{time.time() - t0:.0f}s")


def _stage_job(a):
    kind, spec, f, mesh_n, extra = a
    if kind == "freq":
        return kind, evaluate(spec, f, mesh_n, with_loss=True, floquet=True)
    if kind == "mesh":
        return kind, evaluate(spec, f, mesh_n, with_loss=False)
    if kind == "gap":
        s = dict(spec)
        s["major"] = spec["major"] + extra  # extra = -(delta gap)
        return kind, evaluate(s, f, mesh_n, with_loss=False)
    if kind == "rot":
        s = dict(spec)
        s["angle"] = extra
        mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
        # board construction: the OUTLINE stays fixed (4-fold super-ellipse), only the slot turns.  (Turning the whole
        # cell is NOT equivalent for a super-ellipse: rel. error 1.1-1.3, results/phase30/refine/rotation_check.json.)
        def _cell(slot_angle):
            geo = EllipticalPatch(spec["major"], spec["major"] * spec["ratio"], 0.0, spec["slot_l"], spec["slot_w"], slot_angle, False, EXPO)
            return UnitCell(RectangularLattice(PITCH, PITCH), geo, SUB, Metal(), backing="pec")

        c0 = _cell(0.0)
        ca = _cell(extra)
        ks = free_space_wavenumber(f) * 1.1
        th = extra
        ya = ts.bloch_sheet_tensor(ca, f, mesh, ks, 0.0)
        y0 = ts.bloch_sheet_tensor(c0, f, mesh, ks * np.cos(-th), ks * np.sin(-th))
        pred = ts.rotate_tensor(y0, th)
        return kind, {**spec, "angle": extra, "f_hz": f, "mesh": mesh_n, "rot_rel_error": float(np.linalg.norm(ya - pred) / np.linalg.norm(ya))}
    raise ValueError(kind)


def run_stage(out_dir, picks, workers=20):
    out = Path(out_dir)
    jobs = []
    for pi, spec in enumerate(picks):
        for f in (10.3e9, 10.54e9, 10.8e9):
            jobs.append(("freq", spec, f, 48, pi))
        for mn in (48, 64, 80):
            jobs.append(("mesh", spec, F0, mn, pi))
        for dg in (-0.025e-3, 0.025e-3):
            jobs.append(("gap", spec, F0, 48, -dg))  # gap +dg -> major -dg
        for a in (np.radians(22.5), np.radians(45.0), np.radians(90.0)):
            jobs.append(("rot", spec, F0, 48, a))
    owners = [j[4] if j[0] in ("freq", "mesh") else None for j in jobs]
    t0 = time.time()
    with Pool(workers) as pool:
        res = pool.map(_stage_job, jobs, chunksize=1)
    pk = []
    for pi, spec in enumerate(picks):
        pk.append({"spec": spec, "freq": [], "mesh": [], "gap": [], "rot": []})
    for (kind, r), j in zip(res, jobs, strict=True):
        spec = j[1]
        pi = next(i for i, s in enumerate(picks) if s == spec)
        pk[pi][kind].append({**r, "gap_delta_mm": (-j[4] * 1e3) if kind == "gap" else None})
    del owners
    (out / "refine_stage.json").write_text(json.dumps({"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"], "picks": pk, "seconds": time.time() - t0}, indent=1))
    print("done", len(res), f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["grid", "stage"])
    ap.add_argument("out_dir")
    ap.add_argument("--mesh", type=int, default=48)
    ap.add_argument("--workers", type=int, default=44)
    ap.add_argument("--picks")
    a = ap.parse_args()
    if a.cmd == "grid":
        run_grid(a.out_dir, a.mesh, a.workers)
    else:
        run_stage(a.out_dir, json.loads(Path(a.picks).read_text()), a.workers)
