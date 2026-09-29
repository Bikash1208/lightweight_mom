"""Phase 12D Task A: the dipole/current -> equivalent-S11 conversion must reproduce
``compute_s_parameters`` exactly when fed the exact Floquet coefficient, and closely (small,
documented Riemann-sum residual) when fed Phase 9/12C's own dipole-moment convention."""

from src import solver as solver_mod
from src.basis_functions import MeshSettings
from src.excitation import PlaneWaveExcitation
from src.phase12c import design_point as dp
from src.phase12c.grounded_supercell import solve_isolated_grounded_cell
from src.phase12d.reflection_metric import (
    equivalent_reflection_from_dipole_moment,
    equivalent_reflection_xx,
    reflection_phase_amplitude_error,
)
from src.sparameters import compute_s_parameters

L = 2.4e-3
N = 8


def _exact_j00():
    uc = dp.unit_cell_at(L)
    exc = PlaneWaveExcitation(dp.DESIGN_FREQUENCY_HZ, polarization="x")
    mesh = MeshSettings(N, N, 2 * N)
    system = solver_mod.assemble_system(uc, exc, mesh)
    currents = solver_mod.solve(system.z, system.v)
    return system.basis.current_coefficients(currents, 0.0, 0.0)[0], uc, mesh


def test_exact_floquet_coefficient_reproduces_s11_to_machine_precision():
    j00, uc, mesh = _exact_j00()
    eq = equivalent_reflection_xx(complex(j00), dp.SUBSTRATE, dp.DESIGN_FREQUENCY_HZ, dp.REFERENCE_PLANES)
    sp = compute_s_parameters(
        uc, dp.DESIGN_FREQUENCY_HZ, mesh, dp.REFERENCE_PLANES, warn_grating=False, warn_surface_wave=False
    )
    assert abs(eq.s11_xx - sp.component("S11_xx")) < 1e-12 * abs(sp.component("S11_xx"))


def test_dipole_moment_conversion_matches_true_s11_within_riemann_sum_tolerance():
    iso = solve_isolated_grounded_cell(
        dp.PERIOD_M, dp.PERIOD_M, L, L, dp.SUBSTRATE, dp.METAL, dp.DESIGN_FREQUENCY_HZ, N, 2 * N
    )
    eq = equivalent_reflection_from_dipole_moment(
        iso.dipole.px, dp.PERIOD_M, dp.SUBSTRATE, dp.DESIGN_FREQUENCY_HZ, dp.REFERENCE_PLANES
    )
    sp = compute_s_parameters(
        dp.unit_cell_at(L),
        dp.DESIGN_FREQUENCY_HZ,
        MeshSettings(N, N, 2 * N),
        dp.REFERENCE_PLANES,
        warn_grating=False,
        warn_surface_wave=False,
    )
    assert abs(eq.magnitude - sp.magnitude("S11_xx")) < 0.01  # documented small Riemann-sum residual
    from src.phase import wrapped_difference_deg

    assert abs(wrapped_difference_deg(eq.phase_deg, sp.phase_deg("S11_xx"))) < 1.0


def test_zero_deviation_gives_zero_error():
    phase_err, mag_db_err, iso, emb = reflection_phase_amplitude_error(
        1e-4 + 2e-4j,
        1e-4 + 2e-4j,
        dp.PERIOD_M,
        dp.SUBSTRATE,
        dp.DESIGN_FREQUENCY_HZ,
        dp.REFERENCE_PLANES,
    )
    assert phase_err == 0.0
    assert mag_db_err == 0.0
