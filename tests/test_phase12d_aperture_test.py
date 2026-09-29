"""Phase 12D Task B, at reduced scale (small N, small mesh) so this stays fast while still
exercising the full aperture-level pipeline: lookup-table curve, periodic supercell construction
(reused from phase12c.locality_study, unmodified), full-supercell solve, Floquet-order power
decomposition, and the per-cell phase-error profile."""

from src.phase12d.aperture_test import run_aperture_test


def _small_l_family(mesh_cells_per_subcell=8):
    unit_mm = 2.0 * 8.0 / mesh_cells_per_subcell
    return tuple(round(unit_mm * k, 6) for k in range(1, 4))   # stays well inside an 8mm period


def test_run_aperture_test_reduced_scale_produces_every_field():
    r = run_aperture_test(n_cells=3, mesh_cells_per_subcell=8, l_values_mm=_small_l_family())
    assert r.n_cells == 3
    assert r.n_unknowns > 0
    assert r.runtime_s > 0.0
    assert r.intended_order == (1, 0)
    assert 0.0 <= r.fraction_specular <= 1.0 + 1e-9
    assert 0.0 <= r.fraction_intended <= 1.0 + 1e-9
    assert 0.0 <= r.total_propagating_fraction <= 1.0 + 1e-6
    assert len(r.per_cell_phase_errors) == 3


def test_run_aperture_test_total_propagating_fraction_is_energy_bounded():
    """Never more reflected power than a lossy grounded cell can physically return -- a basic
    sanity bound on the whole Floquet-power machinery at aperture scale."""
    r = run_aperture_test(n_cells=3, mesh_cells_per_subcell=8, l_values_mm=_small_l_family())
    assert r.total_propagating_fraction <= 1.0 + 1e-6


def test_run_aperture_test_spurious_orders_exclude_specular_and_intended():
    r = run_aperture_test(n_cells=3, mesh_cells_per_subcell=8, l_values_mm=_small_l_family())
    flagged = {mn for mn, _ in r.spurious_orders}
    assert (0, 0) not in flagged
    assert r.intended_order not in flagged
