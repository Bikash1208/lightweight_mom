"""Corrected seven-vertical-leg 5.8 GHz MODEL_ONLY cell study.

Run one MoM extraction at a time. Every completed frequency/mesh row is saved
individually and a heartbeat is written every 60 seconds. No AEDT is used.
"""

import argparse
import json
import math
import sys
import threading
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.basis_functions import MeshSettings, RooftopBasis  # noqa: E402
from src.dispersion import tensor_sheet as ts  # noqa: E402
from src.dispersion.elliptical_cells import VerticalLegMeanderSlotPatch  # noqa: E402
from src.dispersion.leakage_model import leakage_alpha_first_order_tensor  # noqa: E402
from src.frequency import free_space_wavenumber  # noqa: E402
from src.geometry import UnitCell  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402

PITCH = 11.6e-3
SUB = Substrate("rt5880", 2.2, 0.0009, 3.175e-3)
SPECS = {
    "baseline_7leg_w0p4": dict(diameter=11.2e-3, leg_count=7, leg_pitch=1.2e-3, leg_depth=3.0e-3, slot_width=0.4e-3),
    "7leg_w0p6": dict(diameter=11.2e-3, leg_count=7, leg_pitch=1.2e-3, leg_depth=3.0e-3, slot_width=0.6e-3),
    "9leg_w0p4": dict(diameter=11.2e-3, leg_count=9, leg_pitch=1.2e-3, leg_depth=3.0e-3, slot_width=0.4e-3),
    "5leg_w0p4": dict(diameter=11.2e-3, leg_count=5, leg_pitch=1.2e-3, leg_depth=3.0e-3, slot_width=0.4e-3),
}


def _complex_pairs(y):
    return [[[float(z.real), float(z.imag)] for z in row] for row in y]


def _restore(pairs):
    return np.array([[complex(*z) for z in row] for row in pairs])


def extract(name, f_hz, mesh_n, out):
    spec = SPECS[name]
    if min(spec["slot_width"], spec["leg_pitch"] - spec["slot_width"]) / (PITCH / mesh_n) < 3:
        raise ValueError("slot/finger gap has fewer than three mesh cells")
    geo = VerticalLegMeanderSlotPatch(**spec)
    cell = UnitCell(RectangularLattice(PITCH, PITCH), geo, SUB, Metal(), backing="pec")
    mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
    basis = RooftopBasis(cell, mesh_n, mesh_n)
    row = {"label": "MODEL_ONLY / NOT_HFSS_VALIDATED", "geometry": name, "parameters_si": spec,
           "frequency_hz": f_hz, "pitch_m": PITCH, "substrate": {"name": SUB.name, "eps_r": 2.2, "tan_delta": .0009, "thickness_m": .003175},
           "mesh": {"nx": mesh_n, "ny": mesh_n, "floquet_order": mesh_n,
                    "cells_across_slot": spec["slot_width"] / (PITCH / mesh_n), "unknowns": basis.n_unknowns}}
    t0 = time.time()
    y = ts.sheet_admittance_tensor_for_cell(cell, f_hz, mesh)
    row["sheet_y_siemens"] = _complex_pairs(y)
    row["seconds_extraction"] = time.time() - t0
    k0 = free_space_wavenumber(f_hz)
    axes = {"along_x": 0.0, "across_y": math.pi / 2}
    for axis, phi in axes.items():
        tm = ts.polarized_surface_wave_roots(y, phi, SUB, f_hz, "TM")
        te = ts.polarized_surface_wave_roots(y, phi, SUB, f_hz, "TE")
        hybrid = None
        try:
            beta, character = ts.find_beta_tensor(y, phi, SUB, f_hz)
            hybrid = {"beta_over_k0": beta / k0, "te_over_tm_field": character}
        except Exception as exc:
            hybrid = {"error": f"{type(exc).__name__}: {exc}"}
        row[axis] = {"TM_projection": tm, "TE_projection": te, "hybrid_TM_like": hybrid}
    _, delta, principal_axis, offdiag = ts.decompose(y)
    row["anisotropy"] = {"delta_mS": delta * 1e3, "principal_axis_deg": math.degrees(principal_axis),
                         "reciprocity_error": ts.reciprocity_error(y), "offdiag_asym_siemens": float(offdiag)}
    betas = [row[a]["hybrid_TM_like"].get("beta_over_k0") for a in axes]
    if all(v is not None for v in betas):
        from src.dispersion.sheet_model import _y_slab
        from src.frequency import angular_frequency
        beta_mean = 0.5 * sum(betas) * k0
        y_bar = 1j * 0.5 * (y[0, 0].imag + y[1, 1].imag)
        slab = _y_slab(beta_mean, angular_frequency(f_hz), SUB.eps_complex, k0, SUB.thickness)
        eig = np.linalg.eigvalsh(0.5 * (y.imag + y.imag.T))
        x = [float(np.imag(1 / (1j * b + slab))) for b in eig]
        row["anisotropy"]["M_t_principal"] = abs(x[0] - x[1]) / max(abs(x[0] + x[1]), 1e-30)
        if abs(f_hz - 5.8e9) < 1:
            row["alpha_radiating_potential_np_per_m"] = leakage_alpha_first_order_tensor(
                beta_mean, y_bar, delta, math.radians(30), SUB, f_hz)
            row["alpha_note"] = "First-order rotating-axis modulation potential at theta=30, not uniform-cell leakage"
    row["status"] = "COMPLETE"
    out.write_text(json.dumps(row, indent=2) + "\n", encoding="utf-8")
    return row


def main(mode, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    if mode == "screen":
        jobs = [(name, 5.8e9, 88) for name in SPECS]
    elif mode == "dispersion":
        jobs = [("baseline_7leg_w0p4", f * 1e9, 88) for f in (5.8, 5.3, 5.6, 6.0, 6.3, 4.5, 5.0, 6.6, 7.0)]
    elif mode == "edge":
        jobs = [("baseline_7leg_w0p4", f * 1e9, 88) for f in (6.02, 6.04, 6.06, 6.08, 6.10)]
    else:
        jobs = [("baseline_7leg_w0p4", 5.8e9, 176)]
    state = {"status": "STARTING", "mode": mode, "total": len(jobs), "complete": 0, "active": None}
    lock, done = threading.Lock(), threading.Event()

    def heartbeat():
        while not done.is_set():
            with lock:
                state["heartbeat_utc"] = datetime.now(timezone.utc).isoformat()
                payload = json.dumps(state, indent=2) + "\n"
            temp = outdir / f"status_{mode}.tmp"
            temp.write_text(payload, encoding="utf-8")
            temp.replace(outdir / f"status_{mode}.json")
            done.wait(60)

    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    try:
        for i, (name, f_hz, mesh_n) in enumerate(jobs):
            path = outdir / f"{name}_{f_hz/1e9:.2f}GHz_n{mesh_n}.json"
            with lock:
                state.update(status="SOLVING", active=path.name, complete=i)
            if path.exists():
                old = json.loads(path.read_text(encoding="utf-8"))
                if old.get("status") == "COMPLETE" and old.get("parameters_si") == SPECS[name] and old.get("mesh", {}).get("nx") == mesh_n:
                    print(f"REUSE {path.name}", flush=True)
                    continue
            row = extract(name, f_hz, mesh_n, path)
            print(f"COMPLETE {path.name} {row['seconds_extraction']:.1f}s", flush=True)
        with lock:
            state.update(status="COMPLETE", active=None, complete=len(jobs))
    except BaseException:
        with lock:
            state.update(status="FAILED", traceback=traceback.format_exc())
        raise
    finally:
        done.set()
        thread.join(timeout=2)
        with lock:
            state["heartbeat_utc"] = datetime.now(timezone.utc).isoformat()
            payload = json.dumps(state, indent=2) + "\n"
        (outdir / f"status_{mode}.json").write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("screen", "dispersion", "edge", "mesh"))
    parser.add_argument("outdir", type=Path)
    args = parser.parse_args()
    main(args.mode, args.outdir)
