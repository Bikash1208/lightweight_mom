"""Phase 12D Task B: the Floquet-order power decomposition, validated on the two cases the brief
requires -- a uniform array (must reproduce |S11|**2 exactly, all other orders exactly zero except
multiples of n) and a trivial case with a known answer (a full-cell PEC sheet: uniform current, so
every order but (0,0) is EXACTLY zero, and |S11| = 1, full reflection)."""

from src import solver as solver_mod
from src.basis_functions import MeshSettings
from src.excitation import PlaneWaveExcitation
from src.geometry import Geometry, UnitCell
from src.green_function import LayeredGreenFunction
from src.phase12c import design_point as dp
from src.phase12c.grounded_supercell import build_grounded_supercell
from src.phase12d.floquet_power import floquet_order_powers
from src.sparameters import compute_s_parameters

N = 6


def test_full_pec_sheet_puts_all_power_in_specular_order():
    geom = Geometry(kind="rectangle", length=dp.PERIOD_M, width=dp.PERIOD_M)
    uc = UnitCell(dp.LATTICE, geom, dp.SUBSTRATE, dp.METAL, backing="pec")
    exc = PlaneWaveExcitation(dp.DESIGN_FREQUENCY_HZ, polarization="x")
    mesh = MeshSettings(N, N, 2 * N)
    system = solver_mod.assemble_system(uc, exc, mesh)
    currents = solver_mod.solve(system.z, system.v)
    green = LayeredGreenFunction(dp.SUBSTRATE, dp.DESIGN_FREQUENCY_HZ, backing="pec")
    result = floquet_order_powers(system.basis, currents, green, dp.LATTICE, max_order=3)

    assert abs(result.fraction_at(0, 0) - 1.0) < 1e-9
    assert abs(result.total_reflected_fraction - 1.0) < 1e-9
    for m in range(-3, 4):
        for n in range(-3, 4):
            if (m, n) == (0, 0):
                continue
            assert abs(result.fraction_at(m, n)) < 1e-20


def test_uniform_supercell_reproduces_isolated_s11_exactly():
    n_cells = 3
    l_m = 6.0e-3
    grid = [[(l_m, l_m)] * n_cells for _ in range(n_cells)]
    supercell = build_grounded_supercell(n_cells, dp.PERIOD_M, dp.PERIOD_M, dp.SUBSTRATE, dp.METAL, grid)
    exc = PlaneWaveExcitation(dp.DESIGN_FREQUENCY_HZ, polarization="x")
    mesh_cells = n_cells * N
    mesh = MeshSettings(mesh_cells, mesh_cells, 2 * mesh_cells)
    system = solver_mod.assemble_system(supercell, exc, mesh)
    currents = solver_mod.solve(system.z, system.v)
    green = LayeredGreenFunction(dp.SUBSTRATE, dp.DESIGN_FREQUENCY_HZ, backing="pec")
    result = floquet_order_powers(system.basis, currents, green, supercell.lattice, max_order=2 * n_cells)

    sp = compute_s_parameters(
        dp.unit_cell_at(l_m),
        dp.DESIGN_FREQUENCY_HZ,
        MeshSettings(N, N, 2 * N),
        dp.REFERENCE_PLANES,
        warn_grating=False,
        warn_surface_wave=False,
    )
    true_fraction = sp.magnitude("S11_xx") ** 2
    assert abs(result.fraction_at(0, 0) - true_fraction) < 1e-9

    for m in range(-2 * n_cells, 2 * n_cells + 1):
        for n in range(-2 * n_cells, 2 * n_cells + 1):
            if m % n_cells != 0 or n % n_cells != 0:
                assert abs(result.fraction_at(m, n)) < 1e-15


def test_phase_deg_at_is_defined_and_finite():
    n_cells = 3
    l_m = 6.0e-3
    grid = [[(l_m, l_m)] * n_cells for _ in range(n_cells)]
    supercell = build_grounded_supercell(n_cells, dp.PERIOD_M, dp.PERIOD_M, dp.SUBSTRATE, dp.METAL, grid)
    exc = PlaneWaveExcitation(dp.DESIGN_FREQUENCY_HZ, polarization="x")
    mesh_cells = n_cells * N
    mesh = MeshSettings(mesh_cells, mesh_cells, 2 * mesh_cells)
    system = solver_mod.assemble_system(supercell, exc, mesh)
    currents = solver_mod.solve(system.z, system.v)
    green = LayeredGreenFunction(dp.SUBSTRATE, dp.DESIGN_FREQUENCY_HZ, backing="pec")
    result = floquet_order_powers(system.basis, currents, green, supercell.lattice, max_order=2 * n_cells)
    phase = result.phase_deg_at(0, 0)
    assert -180.0 < phase <= 180.0
