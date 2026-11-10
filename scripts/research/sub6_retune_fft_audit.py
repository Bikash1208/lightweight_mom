"""Check existing FFT assembler against unchanged production sum at real Bloch k."""
import json
import os
import sys
from pathlib import Path
from types import SimpleNamespace

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0,str(Path(__file__).resolve().parent))
import sub6_iter1_island_retune as R
import numpy as np
from src.dispersion.leaky_supercell_mom import assemble_z_complex

geo = R.Island(R.BASE,96)
cell = R.UnitCell(R.RectangularLattice(.0116,.0116),geo,R.S.make_sub(R.SUB),R.Metal(),backing="pec")
basis = R.RooftopBasis(cell,96,96)
indices = np.concatenate([np.linspace(0,basis.n_x-1,16,dtype=int),
                          np.linspace(basis.n_x,basis.n_unknowns-1,16,dtype=int)])
view = SimpleNamespace(n_unknowns=len(indices),n_x=16,lattice=basis.lattice,
                       dx=basis.dx,dy=basis.dy,center_x=basis.center_x[indices],
                       center_y=basis.center_y[indices])
view.fourier_transform = lambda kx,ky:R.RooftopBasis.fourier_transform(view,kx,ky)
rows = []
for kt,d,order in [(1.017,"x",128),(1.097,"y",128),(1.017,"y",192),(1.097,"x",192)]:
    k0 = R.free_space_wavenumber(5.8e9)
    kx,ky = (kt*k0,0) if d=="x" else (0,kt*k0)
    direct = R.assemble_impedance_matrix(view,R.LayeredGreenFunction(cell.substrate,5.8e9,backing="pec"),kx,ky,order)
    fast = assemble_z_complex(basis,cell.substrate,5.8e9,kx,ky,order,order)
    sample = fast[np.ix_(indices,indices)]
    rel = float(np.max(abs(sample-direct))/np.max(abs(direct)))
    rows.append(dict(kt=kt,direction=d,floquet_order=order,sampled_entries=int(direct.size),max_abs_Z_difference=float(np.max(abs(sample-direct))),max_norm_relative_Z_difference=rel))
    del fast
    print(rows[-1],flush=True)
    if rel > 1e-10:
        raise RuntimeError("HARD_STOP_FFT_ASSEMBLY_AUDIT_FAILED")
(R.OUT/"fft_assembly_audit.json").write_text(json.dumps(dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED",method="existing src.dispersion.leaky_supercell_mom.assemble_z_complex; unchanged src",rows=rows),indent=2),encoding="utf-8")
