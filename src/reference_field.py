"""Phase 8A Sec. 7: generic reference-field and beam-definition input validation for Mode B
(design-specific) feasibility analysis.

This project has no established reference-field file format or beam-steering convention of its
own (grep-verified: no `aperture`/`beam`/`reference_field` module anywhere in `src/`;
`README.md`: "Holographic synthesis is not implemented"). The minimal schema below is therefore
defined *here*, for this gate's own input-validation purposes only -- it is not a claim that this
is "the project's" reference-field convention, because no such convention exists yet.
`feasibility_spatial.check_mapping_method_defined` still fails regardless of whether a reference
field parses cleanly against this schema, because no *validated* method exists to map its samples
onto unit-cell centres (Sec. 31).

Reference-field file (CSV), one row per spatial sample:
    x_mm, y_mm, frequency_GHz, Ex_real, Ex_imag, Ey_real, Ey_imag
"""

import csv
import math
from dataclasses import dataclass
from pathlib import Path

REQUIRED_REFERENCE_FIELD_COLUMNS = (
    "x_mm", "y_mm", "frequency_GHz", "Ex_real", "Ex_imag", "Ey_real", "Ey_imag",
)


@dataclass(frozen=True)
class ReferenceFieldValidation:
    ok: bool
    errors: tuple
    warnings: tuple
    n_samples: int
    x_min_mm: float | None
    x_max_mm: float | None
    y_min_mm: float | None
    y_max_mm: float | None
    dx_mm: float | None
    dy_mm: float | None
    frequencies_ghz: tuple


def _invalid(errors, n_samples=0, frequencies_ghz=()):
    return ReferenceFieldValidation(
        False, tuple(errors), (), n_samples, None, None, None, None, None, None, tuple(frequencies_ghz),
    )


def validate_reference_field(path, expected_frequency_ghz=None, rtol=1e-6):
    """Sec. 7: file format, coordinate/grid sanity, missing/duplicate samples, NaN/Inf, and
    (if requested) that the field actually has data at the requested frequency. Never repairs,
    resamples, or fills in missing data -- only reports."""
    path = Path(path)
    if not path.exists():
        return _invalid([f"reference-field file not found: {path}"])

    try:
        with open(path, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
    except UnicodeDecodeError as e:
        return _invalid([f"{path}: not a UTF-8 text file: {e}"])

    if not rows:
        return _invalid(["reference-field file has no data rows"])

    missing_cols = [c for c in REQUIRED_REFERENCE_FIELD_COLUMNS if c not in rows[0]]
    if missing_cols:
        return _invalid([f"missing required column(s): {missing_cols}"], n_samples=len(rows))

    errors = []
    xs, ys, freqs = [], [], set()
    seen_xy = {}
    for lineno, row in enumerate(rows, start=2):
        try:
            x = float(row["x_mm"]); y = float(row["y_mm"]); f = float(row["frequency_GHz"])
            ex_r, ex_i = float(row["Ex_real"]), float(row["Ex_imag"])
            ey_r, ey_i = float(row["Ey_real"]), float(row["Ey_imag"])
        except (TypeError, ValueError):
            errors.append(f"row {lineno}: non-numeric or missing field")
            continue
        for name, v in (
            ("x_mm", x), ("y_mm", y), ("frequency_GHz", f),
            ("Ex_real", ex_r), ("Ex_imag", ex_i), ("Ey_real", ey_r), ("Ey_imag", ey_i),
        ):
            if math.isnan(v) or math.isinf(v):
                errors.append(f"row {lineno}: {name} is NaN/Inf")
        key = (round(x, 9), round(y, 9))
        if key in seen_xy:
            errors.append(f"row {lineno}: duplicate sample at (x={x}, y={y}) mm (first seen row {seen_xy[key]})")
        else:
            seen_xy[key] = lineno
        xs.append(x); ys.append(y); freqs.add(round(f, 9))

    if expected_frequency_ghz is not None and not any(
        abs(f - expected_frequency_ghz) <= rtol * max(abs(expected_frequency_ghz), 1.0) for f in freqs
    ):
        errors.append(
            f"reference field has no samples at the requested frequency {expected_frequency_ghz} "
            f"GHz (found: {sorted(freqs)})"
        )

    if not xs:
        return _invalid(errors, n_samples=len(rows), frequencies_ghz=sorted(freqs))

    xs_u, ys_u = sorted(set(round(x, 9) for x in xs)), sorted(set(round(y, 9) for y in ys))
    dx = min(b - a for a, b in zip(xs_u, xs_u[1:])) if len(xs_u) > 1 else None
    dy = min(b - a for a, b in zip(ys_u, ys_u[1:])) if len(ys_u) > 1 else None

    return ReferenceFieldValidation(
        ok=len(errors) == 0, errors=tuple(errors), warnings=(), n_samples=len(rows),
        x_min_mm=min(xs), x_max_mm=max(xs), y_min_mm=min(ys), y_max_mm=max(ys),
        dx_mm=dx, dy_mm=dy, frequencies_ghz=tuple(sorted(freqs)),
    )


# ------------------------------------------------------------------------------ beam definition

SUPPORTED_BEAM_POLARIZATIONS = ("x", "y")
_REQUIRED_BEAM_FIELDS = ("theta_deg", "phi_deg", "weight", "initial_phase_deg", "polarization")


@dataclass(frozen=True)
class BeamValidation:
    ok: bool
    errors: tuple


def validate_beam_definition(beam):
    """Generic sanity checks on one beam spec. theta/phi use the standard spherical convention
    (theta from +z boresight, phi from +x in the x-y plane) -- this project defines no
    beam-steering convention of its own (see module docstring), so this is stated as a minimal,
    generic, physically standard convention, not a claim that the project has settled on it."""
    missing = [k for k in _REQUIRED_BEAM_FIELDS if k not in beam]
    if missing:
        return BeamValidation(False, (f"missing beam field(s): {missing}",))

    errors = []
    theta, phi = beam["theta_deg"], beam["phi_deg"]
    weight, phase0, pol = beam["weight"], beam["initial_phase_deg"], beam["polarization"]
    if not (isinstance(theta, (int, float)) and math.isfinite(theta) and 0.0 <= theta <= 90.0):
        errors.append(f"theta_deg={theta!r} must be finite and in [0, 90]")
    if not (isinstance(phi, (int, float)) and math.isfinite(phi) and 0.0 <= phi < 360.0):
        errors.append(f"phi_deg={phi!r} must be finite and in [0, 360)")
    if not (isinstance(weight, (int, float)) and math.isfinite(weight) and weight > 0.0):
        errors.append(f"weight={weight!r} must be finite and > 0")
    if not (isinstance(phase0, (int, float)) and math.isfinite(phase0)):
        errors.append(f"initial_phase_deg={phase0!r} must be finite")
    if pol not in SUPPORTED_BEAM_POLARIZATIONS:
        errors.append(f"polarization={pol!r} not one of {SUPPORTED_BEAM_POLARIZATIONS}")
    return BeamValidation(len(errors) == 0, tuple(errors))
