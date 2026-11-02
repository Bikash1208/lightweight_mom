"""Phase 30: tensor-rotation check, slot-only vs whole-cell rotation (see phase30_refine.py). Writes results/phase30/refine/rotation_check.json."""
import sys, json, math
sys.path.insert(0, r"C:\lightweight_mom"); sys.path.insert(0, r"C:\lightweight_mom\scripts\research")
from multiprocessing import Pool
import numpy as np
import phase30_refine as R
from src.basis_functions import MeshSettings
from src.dispersion import tensor_sheet as ts
from src.dispersion.elliptical_cells import EllipticalPatch
from src.geometry import UnitCell
from src.periodic_bc import RectangularLattice
from src.materials import Metal
from src.frequency import free_space_wavenumber

def cell(spec, slot_angle, outline_rot=0.0):
    geo = EllipticalPatch(spec["major"], spec["major"]*spec["ratio"], outline_rot, spec["slot_l"], spec["slot_w"], slot_angle, False, R.EXPO)
    return UnitCell(RectangularLattice(R.PITCH, R.PITCH), geo, R.SUB, Metal(), backing="pec")

def job(a):
    spec, th, mode, mesh_n = a
    mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
    f = R.F0; ks = 1.1 * free_space_wavenumber(f)
    if mode == "slot_only":
        ca = cell(spec, th); c0 = cell(spec, 0.0)
    else:  # whole cell rotated
        ca = cell(spec, 0.0, th); c0 = cell(spec, 0.0)
    ya = ts.bloch_sheet_tensor(ca, f, mesh, ks, 0.0)
    y0 = ts.bloch_sheet_tensor(c0, f, mesh, ks*math.cos(-th), ks*math.sin(-th))
    pred = ts.rotate_tensor(y0, th)
    return dict(mode=mode, angle_deg=math.degrees(th), rel_err=float(np.linalg.norm(ya-pred)/np.linalg.norm(ya)))

if __name__ == "__main__":
    picks = json.load(open(r"C:\lightweight_mom\results\phase30\refine\picks.json"))
    jobs = []
    for pi in (0, 5, 3):
        for th in (22.5, 45.0, 90.0):
            for mode in ("slot_only", "whole_cell"):
                jobs.append((picks[pi], math.radians(th), mode, 48))
    with Pool(36) as p:
        res = p.map(job, jobs)
    out = []
    for j, r in zip(jobs, res):
        r["pick"] = [0, 5, 3][len(out) // 6]; out.append(r)
    json.dump(out, open(r"C:\lightweight_mom\results\phase30\refine\rotation_check.json", "w"), indent=1)
    for r in out: print(r)
