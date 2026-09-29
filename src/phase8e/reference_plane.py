"""Sec. 6: reference-plane compatibility between the reference field and the unit-cell's own
reference planes, and (if genuinely required) the angular-spectrum propagation that bridges them --
never an arbitrary phase offset (Sec. 6, 26's explicit prohibitions).

The pattern plane (z=0) is where `phase8b.conventions.COORDINATE_SYSTEM` already fixes the metal
geometry and the S-parameter phase reference (`pattern_plane` field, quoting `CONVENTIONS.md`
"Coordinates and stack"); this is the one fixed fact this module reuses, not re-derives.
"""

from dataclasses import dataclass

import numpy as np

from ..frequency import free_space_wavenumber
from ..periodic_bc import normal_wavenumber

PATTERN_PLANE_Z_M = 0.0   # phase8b.conventions.COORDINATE_SYSTEM.pattern_plane; not re-derived here


@dataclass(frozen=True)
class ReferencePlaneCompatibility:
    z_reference_m: float | None
    z_unit_cell_pattern_m: float
    delta_z_m: float | None
    status: str                # "COMPATIBLE" | "BRIDGE_REFERENCE_PLANE_MISMATCH" | "UNKNOWN_REFERENCE_PLANE"
    required_transform: str


def check_reference_plane_compatibility(z_reference_m, rtol=1e-9):
    """Sec. 6: compares the reference field's own recorded height ``z_reference_m`` (``None`` if
    unknown -- e.g. a synthetic field with no stated height) against the pattern plane z=0."""
    if z_reference_m is None:
        return ReferencePlaneCompatibility(
            None, PATTERN_PLANE_Z_M, None, "UNKNOWN_REFERENCE_PLANE",
            "z_reference_m was not supplied with this reference field; compatibility cannot be assessed.",
        )
    delta_z = z_reference_m - PATTERN_PLANE_Z_M
    if abs(delta_z) <= rtol * max(abs(z_reference_m), 1.0):
        return ReferencePlaneCompatibility(z_reference_m, PATTERN_PLANE_Z_M, delta_z, "COMPATIBLE", "none required")
    return ReferencePlaneCompatibility(
        z_reference_m, PATTERN_PLANE_Z_M, delta_z, "BRIDGE_REFERENCE_PLANE_MISMATCH",
        f"angular-spectrum propagation by delta_z={delta_z*1e3:.6g} mm is required "
        "(propagate_angular_spectrum) -- a scalar phase offset is not sufficient whenever "
        "evanescent content is present (spectral.reference_field_spectrum reports the fraction).",
    )


def propagate_angular_spectrum(field_xy, dx_m, dy_m, frequency_hz, z_reference_m):
    """Sec. 6: the standard angular-spectrum method -- propagates a field known at z=z_reference_m
    down to z=0 (the pattern plane) for a wave travelling toward -z (this project's one incident-
    wave direction convention, `excitation.PlaneWaveExcitation`/`CONVENTIONS.md`: k_inc=(kx0,ky0,
    -kz0) under `exp(-j*k.r)`, so U(0) = U(z_ref)*exp(-j*kz*z_ref) with `kz` on the physical branch
    (`periodic_bc.normal_wavenumber`, unmodified) -- evanescent components (kt > k0) then decay
    correctly (Im(kz)<0) rather than needing a separate special case. This is a real derivation
    from an existing project convention, not an invented/arbitrary phase shift (Sec. 6, 26)."""
    field_xy = np.asarray(field_xy, dtype=complex)
    ny, nx = field_xy.shape
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx_m)
    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=dy_m)
    kx_grid, ky_grid = np.meshgrid(kx, ky, indexing="xy")
    k0 = free_space_wavenumber(frequency_hz)
    kz = normal_wavenumber(k0, kx_grid, ky_grid)

    spectrum = np.fft.fft2(field_xy)
    propagated_spectrum = spectrum * np.exp(-1j * kz * z_reference_m)
    return np.fft.ifft2(propagated_spectrum)
