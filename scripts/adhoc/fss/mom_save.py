import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))
from fssgeom import H, unit_cell

from src.basis_functions import MeshSettings, RooftopBasis
from src.excitation import PlaneWaveExcitation, ReferencePlanes, stack_tangential_field
from src.green_function import LayeredGreenFunction
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.solver import solve
from src.sparameters import compute_s_parameters

n = int(sys.argv[1])
out = sys.argv[3] if len(sys.argv) > 3 else "coef"
uc = unit_cell()
planes = ReferencePlanes(top_air_thickness=0.0, bottom_air_thickness=0.0, substrate_thickness=H)
for f in [float(x) * 1e9 for x in sys.argv[2].split(",")]:
    t = time.time()
    green = LayeredGreenFunction(uc.substrate, f, backing="air")
    basis = RooftopBasis(uc, n, n)
    Z = assemble_impedance_matrix(basis, green, 0.0, 0.0, n)
    exc = PlaneWaveExcitation(f, polarization="x")
    e_stack = stack_tangential_field(green, exc, "top")
    v = assemble_excitation_vector(basis, e_stack, 0.0, 0.0)
    I = solve(Z, v[:, None])[:, 0]
    J00 = basis.current_coefficients(I, 0.0, 0.0)
    from src.sparameters import scattering_from_fields

    mean = np.zeros((2, 2, 2), complex)
    mean[0, 0] = J00
    # full S for reference (same as compute_s_parameters, both pols/ports): reuse function on its own
    v2 = np.stack(
        [
            assemble_excitation_vector(
                basis, stack_tangential_field(green, PlaneWaveExcitation(f, polarization=p), "top"), 0.0, 0.0
            )
            for p in "xy"
        ],
        1,
    )
    I2 = solve(Z, v2)
    mean = np.zeros((2, 2, 2), complex)
    for k in range(2):
        mean[0, k] = basis.current_coefficients(I2[:, k], 0.0, 0.0)
    s = scattering_from_fields(green, planes, mean)
    if len(sys.argv) > 4:
        s_ref = compute_s_parameters(uc, f, MeshSettings(n, n, n), planes, warn_grating=False)
        print("check vs compute_s_parameters", abs(s.s11[0, 0] - s_ref.s11[0, 0]), abs(s.s21[0, 0] - s_ref.s21[0, 0]))
    np.savez(
        f"{out}/f{f / 1e9:.3f}.npz",
        f=f,
        coeffs=I,
        cx=basis.center_x,
        cy=basis.center_y,
        dirn=basis.direction,
        dx=basis.dx,
        dy=basis.dy,
        J00=J00,
        e_stack=np.asarray(e_stack),
        s11=s.s11,
        s21=s.s21,
        n=n,
    )
    print(
        f"saved {f / 1e9:.3f} GHz  S11x {20 * np.log10(abs(s.s11[0, 0])):.3f} dB  S21x {20 * np.log10(abs(s.s21[0, 0])):.3f} dB  {time.time() - t:.0f}s",
        flush=True,
    )
