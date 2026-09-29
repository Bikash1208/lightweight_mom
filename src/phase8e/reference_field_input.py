"""Sec. 3-4: reference-field and object-wave inputs.

This project has no real HFSS/measured reference-field export anywhere in the repository
(grep-verified, same finding as `src/reference_field.py`'s module docstring). Sec. 3's own
"known historical reference-grid format" (``Grid Output Min/Max``, ``Grid Size``, ``CS =
RelativeCS2``) is quoted context for what a *real* file would look like -- it is deliberately
**not** hard-coded here as a default grid (Sec. 3's explicit instruction); every grid/height value
below is a caller-supplied argument.

Two reference-field sources are supported, both explicit about what they are:

- a real file, loaded via `src.reference_field.validate_reference_field` (Phase 8A, unmodified) and
  gridded here -- tagged ``source="FILE:<path>"``, amplitude always available (the file schema is
  complex E only, `REQUIRED_REFERENCE_FIELD_COLUMNS`);
- a synthetic demonstration field (no real file exists to load), tagged ``source="SYNTHETIC_DEMO_*"``
  so it can never be mistaken for real hardware/measurement data anywhere downstream.

Candidate A (Fong et al. 2010) needs a **guided surface wave** reference (`monopole_reference_field`,
reused unmodified from `phase8b.reference_wave`, still requiring an externally-supplied ``k_swg``).
Candidates B/C (reflectarray/transmitarray-style phase holograms) use the standard **spherical
point-source feed** model instead (`spherical_feed_reference_field`) -- a different, also standard,
architecture-appropriate reference field, not a substitute for Candidate A's.
"""

from dataclasses import dataclass

import numpy as np

from ..frequency import free_space_wavenumber
from ..phase8b.reference_wave import monopole_reference_field
from ..reference_field import validate_reference_field


@dataclass(frozen=True)
class ReferenceFieldGrid:
    x_m: np.ndarray  # 2-D
    y_m: np.ndarray  # 2-D
    Ex: np.ndarray | None  # 2-D complex, or None if unavailable
    Ey: np.ndarray | None
    amplitude_available: bool
    source: str  # "FILE:<path>" | "SYNTHETIC_DEMO_SPHERICAL_FEED" | "SYNTHETIC_DEMO_MONOPOLE" | "PHASE_ONLY_INPUT"
    dx_m: float | None
    dy_m: float | None
    z_reference_m: float | None  # height of this field above the pattern plane z=0, if known
    phase_deg: np.ndarray | None = None  # arg(Ex) if amplitude is available, else the supplied phase-only input

    def phase(self):
        """arg(E_ref) [deg], the one quantity every candidate formulation's phase-only path
        (Formulations B/C) actually needs -- available whether or not amplitude is."""
        if self.phase_deg is not None:
            return self.phase_deg
        if self.Ex is not None:
            return np.degrees(np.angle(self.Ex))
        raise ValueError("ReferenceFieldGrid: neither Ex nor phase_deg is available")


def load_reference_field_from_csv(path, expected_frequency_ghz=None):
    """Sec. 3: load a real reference-field CSV via the existing, unmodified
    `reference_field.validate_reference_field`; raises ValueError (not a silent fallback to
    synthetic data) if the file fails validation."""
    import csv
    from pathlib import Path

    validation = validate_reference_field(path, expected_frequency_ghz=expected_frequency_ghz)
    if not validation.ok:
        raise ValueError(f"reference field file {path} failed validation: {validation.errors}")

    with open(Path(path), newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if expected_frequency_ghz is not None:
        rows = [
            r
            for r in rows
            if abs(float(r["frequency_GHz"]) - expected_frequency_ghz) <= 1e-6 * max(abs(expected_frequency_ghz), 1.0)
        ]

    xs = sorted({round(float(r["x_mm"]), 9) for r in rows})
    ys = sorted({round(float(r["y_mm"]), 9) for r in rows})
    x_grid, y_grid = np.meshgrid(np.array(xs) * 1e-3, np.array(ys) * 1e-3, indexing="xy")
    ex_grid = np.full(x_grid.shape, np.nan + 0j, dtype=complex)
    ey_grid = np.full(x_grid.shape, np.nan + 0j, dtype=complex)
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}
    for r in rows:
        i, j = yi[round(float(r["y_mm"]), 9)], xi[round(float(r["x_mm"]), 9)]
        ex_grid[i, j] = complex(float(r["Ex_real"]), float(r["Ex_imag"]))
        ey_grid[i, j] = complex(float(r["Ey_real"]), float(r["Ey_imag"]))
    if np.any(np.isnan(ex_grid)) or np.any(np.isnan(ey_grid)):
        raise ValueError(f"reference field file {path}: not a complete rectangular grid at the requested frequency")

    return ReferenceFieldGrid(
        x_m=x_grid,
        y_m=y_grid,
        Ex=ex_grid,
        Ey=ey_grid,
        amplitude_available=True,
        source=f"FILE:{path}",
        dx_m=(validation.dx_mm * 1e-3) if validation.dx_mm else None,
        dy_m=(validation.dy_mm * 1e-3) if validation.dy_mm else None,
        z_reference_m=None,
    )


def spherical_feed_reference_field(x_m, y_m, feed_height_m, frequency_hz, amplitude_at_1m=1.0):
    """Sec. 3-4 (synthetic demonstration, no real data exists): the standard spherical point-
    source feed model used for reflectarray/transmitarray phase-hologram reference fields
    (Candidates B/C; e.g. Huang & Encinar, 'Reflectarray Antennas,' Ch. 3): a feed on-axis at
    height ``feed_height_m`` above the surface radiates `E_ref(x,y) = A/R * exp(-j*k0*R)`,
    `R = sqrt(x^2+y^2+h^2)` (`exp(-j*k.r)` convention, `phase8b.conventions.
    PLANE_WAVE_PHASE_CONVENTION` extended to a spherical wave -- the standard far-field-of-a-
    point-source form), **evaluated at the pattern plane z=0** (``R`` already measures the
    distance from the feed to a point ON that plane) -- so this field's own ``z_reference_m`` is
    0.0, i.e. already reference-plane-*compatible* by construction, unlike a field recorded at a
    genuinely different height (see ``propagate_angular_spectrum`` for that case). Scalar field,
    x-polarized by convention (returned as ``Ex``; ``Ey`` is exactly zero -- an idealized feed,
    not a claim about any real antenna's cross-polarization)."""
    if feed_height_m <= 0.0:
        raise ValueError(f"feed_height_m must be > 0, got {feed_height_m!r}")
    x_m, y_m = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
    k0 = free_space_wavenumber(frequency_hz)
    r = np.sqrt(x_m**2 + y_m**2 + feed_height_m**2)
    ex = amplitude_at_1m / r * np.exp(-1j * k0 * r)
    return ReferenceFieldGrid(
        x_m=x_m,
        y_m=y_m,
        Ex=ex,
        Ey=np.zeros_like(ex),
        amplitude_available=True,
        source="SYNTHETIC_DEMO_SPHERICAL_FEED",
        dx_m=float(x_m[0, 1] - x_m[0, 0]) if x_m.ndim == 2 and x_m.shape[1] > 1 else None,
        dy_m=float(y_m[1, 0] - y_m[0, 0]) if y_m.ndim == 2 and y_m.shape[0] > 1 else None,
        z_reference_m=0.0,
    )


def monopole_demo_reference_field(x_m, y_m, k_swg, amplitude_at_1m=1.0):
    """Candidate A's reference field (`phase8b.reference_wave.monopole_reference_field`, reused
    unmodified). ``k_swg`` still has no default -- this project cannot compute it."""
    x_m, y_m = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
    ex = monopole_reference_field(x_m, y_m, k_swg, amplitude_at_1m)
    return ReferenceFieldGrid(
        x_m=x_m,
        y_m=y_m,
        Ex=ex,
        Ey=np.zeros_like(ex),
        amplitude_available=True,
        source="SYNTHETIC_DEMO_MONOPOLE",
        dx_m=float(x_m[0, 1] - x_m[0, 0]) if x_m.ndim == 2 and x_m.shape[1] > 1 else None,
        dy_m=float(y_m[1, 0] - y_m[0, 0]) if y_m.ndim == 2 and y_m.shape[0] > 1 else None,
        z_reference_m=0.0,  # a surface wave propagates *on* the surface, z=0 by construction
    )


def phase_only_reference_field(x_m, y_m, phase_deg):
    """Sec. 3's phase-only case: amplitude is explicitly marked unavailable (``amplitude_available
    =False``) and never silently reconstructed (``Ex``/``Ey`` are left `None` -- callers that need
    a complex field for a candidate formulation that needs amplitude must fail loudly, not receive
    a fabricated unit-amplitude field)."""
    x_m, y_m, phase_deg = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float), np.asarray(phase_deg, dtype=float)
    return ReferenceFieldGrid(
        x_m=x_m,
        y_m=y_m,
        Ex=None,
        Ey=None,
        amplitude_available=False,
        source="PHASE_ONLY_INPUT",
        dx_m=None,
        dy_m=None,
        z_reference_m=None,
        phase_deg=phase_deg,
    )
