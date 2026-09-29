"""Load configs/config.yaml and print the parsed unit cell (SI internally, mm/GHz for display).

Run from the project root:  python examples/print_config.py [path/to/config.yaml]
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import units
from src.config import load_config
from src.periodic_bc import normal_wavenumber
from src.solver import build_system, solve
from src.sparameters import s_parameters_from_config


def main(path):
    cfg = load_config(path)
    cell, exc = cfg.unit_cell, cfg.excitation
    lat, geo, sub = cell.lattice, cell.geometry, cell.substrate

    print(f"Config: {path}   (time convention exp(+jwt), SI internally)")
    print(f"  period        : {units.to_mm(lat.period_x):.4g} x {units.to_mm(lat.period_y):.4g} mm")
    print(
        f"  substrate     : {sub.name}, eps_r={sub.eps_r}, tan_d={sub.tan_delta}, h={units.to_mm(sub.thickness):.4g} mm"
    )
    print(
        f"  geometry      : {geo.kind} {units.to_mm(geo.length):.4g} x {units.to_mm(geo.width):.4g} mm, "
        f"rot={units.to_deg(geo.rotation):.4g} deg, metal fill={cell.metal_fill_fraction:.3f}"
    )
    print(
        f"  frequency     : {units.to_ghz(exc.frequency):.6g} GHz  "
        f"(lambda0={units.to_mm(exc.wavelength):.4g} mm, k0={exc.k0:.5g} rad/m)"
    )
    print(f"  period/lambda0: {lat.period_x / exc.wavelength:.3f} x {lat.period_y / exc.wavelength:.3f}")
    print(f"  substrate k   : {sub.wavenumber(exc.frequency):.5g} rad/m")

    m, n = lat.mode_grid(3)
    kx, ky = lat.transverse_wavenumbers(*exc.transverse_wavenumber, m, n)
    kz = normal_wavenumber(exc.k0, kx, ky)
    print(f"  propagating free-space Floquet orders (|m|,|n| <= 3): {int((kz.imag == 0).sum())}")

    rp = cfg.reference_planes
    print(f"  reference planes: z_top = {units.to_mm(rp.z_top):.4g} mm, z_bottom = {units.to_mm(rp.z_bottom):.4g} mm")
    msh = cfg.mesh
    system = build_system(cfg)
    current = solve(system.z, system.v)
    j00 = system.basis.current_coefficients(current, system.kx0, system.ky0)
    print(
        f"  MoM: {msh.cells_x} x {msh.cells_y} mesh, Floquet order {msh.floquet_order}, "
        f"{system.basis.n_unknowns} unknowns ({system.basis.n_x} x-, {system.basis.n_y} y-directed)"
    )
    print(f"  mean current (0,0 Floquet coeff.): Jx = {j00[0]:.4e} A/m, Jy = {j00[1]:.4e} A/m")

    s = s_parameters_from_config(cfg)
    print("  S-parameters at the reference planes (exp(+jwt); S_ij[received, incident]):")
    for name in ("S11_xx", "S11_yx", "S21_xx", "S21_yx", "S21_yy"):
        print(
            f"    {name}: |S| = {s.magnitude(name):.5f} ({s.magnitude_db(name):8.3f} dB), "
            f"phase = {s.phase_deg(name):8.3f} deg"
        )
    print(f"    lossy-cell power balance |S11|^2+|S21|^2 (x in): {s.power_out_per_incident()[0, 0]:.6f}   (1 - loss)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent / "configs" / "config.yaml")
