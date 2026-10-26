"""Phase 20 R3 item 3: optimisation of the 2-bit RIS cells (OUTSIDE the solver core; Nelder-Mead on coarse aligned-pitch meshes, parallel over
frequency / state / corner). MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.

This module holds the shared evaluation helpers; the per-family drivers are `phase20_r3_opt_cp2.py` and `phase20_r3_opt_lp.py`.
"""

import math
import os
import sys
import warnings
from pathlib import Path

for _k in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.basis_functions import MeshSettings  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.ris.polarization import cartesian_reflection  # noqa: E402
from src.sparameters import compute_s_parameters_oblique  # noqa: E402

PLANES = ReferencePlanes(0.0, None, 1.575e-3)
PITCH = 0.25e-3  # nominal mesh pitch of the search meshes


def mesh_for(period, pitch=PITCH):
    return max(8, int(round(period / pitch)))


def solve_s11(task):
    """task = (cell, state, kit, f_hz, theta_deg, mesh_n) -> s11 (TE, TM) 2x2 complex, or None when the cell is invalid (e.g. a stub reaches
    the cell boundary)."""
    cell, state, kit, f, theta, n = task
    try:
        uc = cell.unit_cell(state, kit)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return compute_s_parameters_oblique(
                uc, f, MeshSettings(n, n, n // 2 + 1), PLANES, math.radians(theta), 0.0
            ).s11
    except ValueError:
        return None


def wrap180(d):
    return (d + 180.0) % 360.0 - 180.0


def db(x):
    return 20 * np.log10(max(abs(x), 1e-9))


def pb_single(cell, kit, freqs, pool, n=None, theta=0.0):
    """One-diode-ON PB condition of a CP cell in the stub frame (state 0: stub along +x): per frequency (|r_xx|, phase_xx, |r_yy|, phase_yy,
    |r_xy|) in the Cartesian frame at normal incidence (theta=0). Returns None if invalid."""
    n = n or mesh_for(cell.period)
    st = cell.states()[0]
    s = pool.map(solve_s11, [(cell, st, kit, f, theta, n) for f in freqs])
    if any(v is None for v in s):
        return None
    out = []
    for v in s:
        r = cartesian_reflection(v, 0.0)
        out.append(
            (abs(r[0, 0]), math.degrees(np.angle(r[0, 0])), abs(r[1, 1]), math.degrees(np.angle(r[1, 1])), abs(r[0, 1]))
        )
    return np.array(out)


def pb_errors(rows, tol_deg=10.0, floor_db=-1.0):
    """Per-frequency normalised violation of the single-state PB condition: 180 +- tol_deg phase difference and both |r| >= floor_db."""
    e = []
    for rxx, pxx, ryy, pyy, _ in rows:
        dphi = abs(wrap180(pxx - pyy))  # 0..180; ideal 180
        e.append(max(0.0, (180.0 - dphi) - tol_deg) / tol_deg + max(0.0, floor_db - min(db(rxx), db(ryy))) / 1.0)
    return np.array(e)


def pb_band(freqs, rows, f0=10e9, tol_deg=10.0, floor_db=-1.0):
    """Contiguous bandwidth (Hz) around f0 satisfying the single-state PB condition on the sampled grid (0 if f0 itself fails)."""
    ok = pb_errors(rows, tol_deg, floor_db) == 0
    freqs = np.asarray(freqs)
    i0 = int(np.argmin(abs(freqs - f0)))
    if not ok[i0]:
        return 0.0
    lo = hi = i0
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    while hi < len(freqs) - 1 and ok[hi + 1]:
        hi += 1
    return float(freqs[hi] - freqs[lo])


# ---------------------------------------------------------------------------------------------------------------------
# LP families: 4-state objective with tolerance corners
# ---------------------------------------------------------------------------------------------------------------------
import itertools  # noqa: E402

from src.ris.cells import DiodeKit  # noqa: E402


def corner_kits(kit):
    """Nominal + one-at-a-time tolerance corners. SMP1345: C_T and L_s +-20 % (both from the datasheet nominal). MADP: C_T +-20 %; its L_s is
    NOT_AVAILABLE (nominal 0), so an L_s corner is the one-sided ASSUMED offset +0.1 nH (ASSUMED_SENSITIVITY_NOT_DEVICE_DATA)."""
    out = {"nominal": kit}

    def sc(**f):
        return DiodeKit(
            kit.name,
            kit.on.scaled(**{k: v for k, v in f.items() if k != "C_T"}) if "L_s" in f else kit.on,
            kit.off.scaled(**f),
        )

    for name, f in (("C_T-20%", {"C_T": 0.8}), ("C_T+20%", {"C_T": 1.2})):
        out[name] = DiodeKit(kit.name, kit.on, kit.off.scaled(**f))
    if kit.on.L_s:
        for name, f in (("L_s-20%", 0.8), ("L_s+20%", 1.2)):
            out[name] = sc(L_s=f)
    else:
        from src.ris.cells import madp_kit

        out["L_s+0.1nH(assumed)"] = madp_kit(0.1e-9)
    return out


def lp_responses(cell, kits, thetas, pool, n=None, basis="linear"):
    """{(corner, theta): (phase_deg[4], amp[4], xpol[4])} of the TM co-polar reflection for the 4 states (TE-TM leakage as cross-pol)."""
    n = n or mesh_for(cell.period)
    tasks, keys = [], []
    for kname, kit in kits.items():
        for th in thetas if kname == "nominal" else (0.0,):  # corners at theta = 0 only (cost)
            for st in cell.states():
                tasks.append((cell, st, kit, 10e9, th, n))
                keys.append((kname, th))
    res = pool.map(solve_s11, tasks)
    out = {}
    for k, s in zip(keys, res, strict=True):
        if s is None:
            return None
        ph, am, xp = out.setdefault(k, ([], [], []))
        if basis == "circular":  # handedness-preserving r_RR; undesired = opposite-handed (plate-like) r_LR and r_RL
            from src.ris.polarization import L, R, circular_reflection

            c = circular_reflection(s)
            ph.append(math.degrees(np.angle(c[R, R])))
            am.append(abs(c[R, R]))
            xp.append(max(abs(c[L, R]), abs(c[R, L])))
        else:
            ph.append(math.degrees(np.angle(s[1, 1])))
            am.append(abs(s[1, 1]))
            xp.append(abs(s[0, 1]))
    return {k: tuple(np.array(v) for v in vals) for k, vals in out.items()}


def lattice_fit(phases):
    """Best assignment of the 4 phases to the ideal 0/90/180/270 lattice (any permutation, any offset): returns (max deviation deg, perm) with
    perm[i] = lattice index of state i."""
    best = (1e9, None)
    for perm in itertools.permutations(range(4)):
        d = np.exp(1j * np.radians(np.asarray(phases) - 90.0 * np.array(perm)))
        c = np.angle(d.sum())
        dev = float(np.max(np.abs(np.degrees(np.angle(d * np.exp(-1j * c))))))
        if dev < best[0]:
            best = (dev, perm)
    return best


def step_error(phases, perm):
    """Max |step - 90| (deg) over the 4 cyclic steps between lattice-adjacent states, with the nominal assignment `perm` held fixed."""
    order = np.argsort(perm)  # state index at lattice position 0..3
    p = np.asarray(phases)[order]
    steps = [wrap180(p[(k + 1) % 4] - p[k]) for k in range(4)]
    # cyclic steps are taken modulo 360 around the circle: each should be +90 (or all -90 for a reversed ladder)
    return float(min(max(abs(s - 90) for s in steps), max(abs(s + 90) for s in steps)))


def lp_score(resp, thetas_report=(0.0, 30.0)):
    """Scalar objective (0 = every target met) + a metrics dict from `lp_responses` output."""
    ph, am, xp = resp[("nominal", 0.0)]
    dev, perm = lattice_fit(ph)
    amp_db = 20 * np.log10(np.maximum(am, 1e-9))
    xp_db = 20 * np.log10(np.maximum(xp.max(), 1e-9))
    e = max(0.0, dev - 15) / 15 + dev / 600
    e += max(0.0, -1.5 - amp_db.min()) / 1.5 + max(0.0, (amp_db.max() - amp_db.min()) - 1.0)
    e += max(0.0, xp_db + 15) / 5
    metrics = {"dev0": dev, "amp_db": amp_db.tolist(), "xpol_db": float(xp_db)}
    if ("nominal", 30.0) in resp:
        ph3, am3, xp3 = resp[("nominal", 30.0)]
        dev3, _ = lattice_fit(ph3)
        a3 = 20 * np.log10(np.maximum(am3, 1e-9))
        e += max(0.0, dev3 - 15) / 15 + max(0.0, -1.5 - a3.min()) / 1.5 + max(0.0, a3.max() - a3.min() - 1.0)
        x3 = 20 * np.log10(max(float(xp3.max()), 1e-9))
        e += max(0.0, x3 + 15) / 5
        metrics.update(dev30=dev3, amp30_db=a3.tolist(), xpol30_db=float(x3))
    for (kname, th), (p, _, _) in resp.items():
        if kname != "nominal" and th == 0.0:
            se = step_error(p, perm)
            metrics[f"step_err[{kname}]"] = se
            e += max(0.0, se - 20) / 20 + se / 1000
    return float(e), metrics
