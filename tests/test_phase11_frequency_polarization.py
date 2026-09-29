"""Sec. 14-15, 28: frequency and polarization dependence of the coupling extent."""

from src.materials import Metal, Substrate
from src.phase11 import frequency_study as freq_mod
from src.phase11 import polarization_study as pol_mod

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
CENTER = (2.8, 2.8)


def test_frequency_study_covers_every_requested_frequency():
    freqs = (10.0e9, 12.5e9, 15.0e9)
    results = freq_mod.run_frequency_study("C3", CENTER, P, SUB, MET, mesh_cells_per_subcell=8,
                                            frequencies_hz=freqs, sizes=(1, 3))
    assert [r.frequency_hz for r in results] == list(freqs)
    for r in results:
        assert set(r.responses_by_n.keys()) == {1, 3}
        assert r.max_shell_increment_relative >= 0


def test_polarization_study_covers_x_and_y():
    results = pol_mod.run_polarization_study(CENTER, P, SUB, MET, 12.5e9, mesh_cells_per_subcell=8,
                                              sizes=(1, 3))
    assert {r.polarization for r in results} == {"x", "y"}


def test_circular_status_is_explicit_not_available():
    assert pol_mod.CIRCULAR_STATUS == "NOT_SUPPORTED_AS_INDEPENDENT_EXCITATION"
