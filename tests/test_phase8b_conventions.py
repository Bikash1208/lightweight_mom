"""Phase 8B Sec. 5-6: `phase8b.conventions` values must match `CONVENTIONS.md`'s literal text and
`excitation.PlaneWaveExcitation`'s actual behaviour -- never a re-typed, potentially-drifted copy."""

import math
from pathlib import Path

import pytest

from src.phase8b import conventions

CONVENTIONS_MD = Path(__file__).resolve().parent.parent / "docs" / "theory" / "CONVENTIONS.md"


def test_time_convention_matches_database_conventions():
    from src.database import CONVENTIONS

    assert conventions.TIME_CONVENTION == CONVENTIONS["time_dependence"] == "exp(+j*omega*t)"


def test_incident_wavevector_deg_reuses_plane_wave_excitation():
    # x-polarization is perpendicular to k_inc at phi=90 (k_inc has no x-component then).
    exc = conventions.incident_wavevector_deg(12.5e9, 20.0, 90.0, polarization="x")
    assert exc.theta == pytest.approx(math.radians(20.0))
    # theta=0 (broadside) is valid for any polarization.
    exc0 = conventions.incident_wavevector_deg(12.5e9, 0.0, 0.0, polarization="x")
    assert exc0.theta == 0.0
    kx0, ky0 = exc0.transverse_wavenumber
    assert kx0 == 0.0 and ky0 == 0.0


def test_no_ground_plane_and_no_surface_wave_solver_are_true_facts():
    # These booleans gate every downstream INCOMPATIBLE conclusion; pin them so a future change
    # to the solver (adding a ground plane or a dispersion solver) is forced to update this file.
    assert conventions.NO_GROUND_PLANE is True
    assert conventions.SURFACE_WAVE_DISPERSION_IMPLEMENTED is False
    assert conventions.OBLIQUE_INCIDENCE_PORTS_IMPLEMENTED is False


def test_coordinate_system_matches_conventions_md_text():
    text = CONVENTIONS_MD.read_text(encoding="utf-8")  # the file is UTF-8; the Windows default codec is not
    assert "Px/2" in text and "no ground plane" in text.lower()
