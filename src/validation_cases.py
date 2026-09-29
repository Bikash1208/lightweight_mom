"""Definitions of the representative unit cells exported and plotted in Phase 6.

Every case is an ordinary `config.yaml`-shaped dict, built and validated by the **existing**
`config.config_from_dict` (Phase 1) — nothing here defines geometry, materials or mesh rules of
its own. A case adds only what Phase 1-5 configs do not carry: an id, a category (for the
selection checklist below) and a one-line description, plus its own frequency grid (`sweep.
FrequencyGrid`, Phase 5, unmodified) since different cases need different bands to show their
behaviour (Sec. "Case selection" below).

All ten cases share one lattice, substrate and reference-plane configuration so the exported
substrate solid is identical across cases and the numbers are easy to compare; only the metal
geometry, mesh (case 3 needs a finer mesh — see below) and frequency grid vary.

Case selection (the six categories requested for Phase 6), chosen by *measuring* the solver's
own output (`scratch scans, not committed`; the concrete numbers below and the qualitative
behaviour are what those scans found, not assumed):

* **non-resonant** (case01): a patch much smaller than the 8 mm cell; |S11_xx| stays below 0.23
  and |S21_xx| above 0.97 over the whole 8-12 GHz band — no resonance in this band.
* **resonant** (case02): a thin dipole whose length was tuned by scanning length at 10 GHz and
  frequency at fixed length until reflection peaked; at length 7 mm, width 1 mm, |S11_xx| rises
  to 0.998 (|S21_xx| falls to 0.013) at 16.4 GHz, well clear of the 6-20 GHz band edges chosen
  for this case specifically to show the full resonance curve.
* **narrow** (case03): a dipole of aspect ratio 4:1 (6 mm x 1.5 mm); needs a finer mesh (32x32)
  than the other cases so the rooftop mesh resolves the narrow dimension at all.
* **wide** (case04): width larger than length (6 mm x 3 mm), the opposite aspect ratio.
* different lengths/widths: cases 1, 4, 5, 6 span four unrelated (length, width) pairs.
* **rotated** (case07, case08): non-zero `rotation_deg`; case08 additionally offsets the feature
  from the cell centre. Both show measurable cross-polarization (|S21_yx| ~ 0.06-0.08), unlike
  every unrotated case (|S21_yx| = 0 to machine precision) — evidence the rotation is doing
  something, not just present in the config.
* **slot / complement** (case09, case10): `geometry.type: slot`, the metal-complement aperture
  formulation (`geometry.py`, Phase 1); case10 is additionally rotated.

Not selected: any case whose mesh a-priori checks (`discretization.py`) would fail is not
included; every case here has 0-3% metal-area staircase error against the analytic feature area
(`check_cell`, no case reaches the default 5% policy threshold) and stays specular-only across
its whole frequency grid (grating onset `c0/8mm` = 37.5 GHz, above every case's top frequency).
"""

from dataclasses import dataclass

from .config import SimulationConfig, config_from_dict
from .sweep import FrequencyGrid

# Shared across every case (see module docstring).
_PERIOD_M = 8.0e-3
_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3}
_REFERENCE_PLANES = {"top_air_thickness": 5.0e-3, "bottom_air_thickness": 5.0e-3}
_DEFAULT_MESH = {"cells_x": 16, "cells_y": 16, "floquet_order": 32}
_FREQ_STEP_HZ = 0.5e9


@dataclass(frozen=True)
class ValidationCase:
    case_id: str
    category: str
    description: str
    config: SimulationConfig
    frequency: FrequencyGrid


def _case(case_id, category, description, geometry, f_start_hz, f_stop_hz, mesh=None):
    data = {
        "period_x": _PERIOD_M,
        "period_y": _PERIOD_M,
        "substrate": _SUBSTRATE,
        "geometry": geometry,
        "excitation": {"frequency": (f_start_hz + f_stop_hz) / 2.0},
        "mesh": mesh or _DEFAULT_MESH,
        "reference_planes": _REFERENCE_PLANES,
    }
    config = config_from_dict(data)
    frequency = FrequencyGrid(f_start_hz, f_stop_hz, _FREQ_STEP_HZ)
    return ValidationCase(case_id, category, description, config, frequency)


CASES = (
    _case(
        "case01_non_resonant_patch",
        "non-resonant",
        "Small rectangular patch (2 x 1 mm); reflection stays low across the whole band.",
        {"type": "rectangle", "length": 2.0e-3, "width": 1.0e-3},
        8.0e9,
        12.0e9,
    ),
    _case(
        "case02_resonant_dipole",
        "resonant",
        "Thin dipole (7 x 1 mm); near-total reflection / transmission null around 16.4 GHz.",
        {"type": "dipole", "length": 7.0e-3, "width": 1.0e-3},
        6.0e9,
        20.0e9,
    ),
    _case(
        "case03_narrow_dipole",
        "narrow",
        "High-aspect-ratio dipole (6 x 1.5 mm, 4:1); uses a finer 32x32 mesh to resolve the width.",
        {"type": "dipole", "length": 6.0e-3, "width": 1.5e-3},
        8.0e9,
        12.0e9,
        mesh={"cells_x": 32, "cells_y": 32, "floquet_order": 48},
    ),
    _case(
        "case04_wide_patch",
        "wide",
        "Wide patch (3 x 6 mm, width > length): the opposite aspect ratio from case03.",
        {"type": "rectangle", "length": 3.0e-3, "width": 6.0e-3},
        8.0e9,
        12.0e9,
    ),
    _case(
        "case05_small_square_patch",
        "different lengths/widths",
        "Square patch (2 x 2 mm) built with the 'patch' geometry alias for 'rectangle'.",
        {"type": "patch", "length": 2.0e-3, "width": 2.0e-3},
        8.0e9,
        12.0e9,
    ),
    _case(
        "case06_large_rect_patch",
        "different lengths/widths",
        "Larger rectangular patch (6 x 3 mm), a length/width pair distinct from every other case.",
        {"type": "rectangle", "length": 6.0e-3, "width": 3.0e-3},
        8.0e9,
        12.0e9,
    ),
    _case(
        "case07_rotated_dipole",
        "rotated",
        "Dipole (6 x 1.2 mm) rotated 30 deg about the cell centre; produces measurable cross-pol.",
        {"type": "dipole", "length": 6.0e-3, "width": 1.2e-3, "rotation_deg": 30.0},
        8.0e9,
        12.0e9,
    ),
    _case(
        "case08_rotated_offset_patch",
        "rotated",
        "Patch (5 x 2 mm) rotated 60 deg and offset from the cell centre (0.5, -0.3 mm).",
        {
            "type": "rectangle",
            "length": 5.0e-3,
            "width": 2.0e-3,
            "rotation_deg": 60.0,
            "center_x": 0.5e-3,
            "center_y": -0.3e-3,
        },
        8.0e9,
        12.0e9,
    ),
    _case(
        "case09_slot_rectangular",
        "slot",
        "Rectangular slot (5 x 2 mm), metal-complement formulation: a solid plate with an aperture.",
        {"type": "slot", "length": 5.0e-3, "width": 2.0e-3},
        8.0e9,
        12.0e9,
    ),
    _case(
        "case10_slot_rotated",
        "slot",
        "Rotated rectangular slot (4 x 1.5 mm, 20 deg); slot geometry combined with rotation.",
        {"type": "slot", "length": 4.0e-3, "width": 1.5e-3, "rotation_deg": 20.0},
        8.0e9,
        12.0e9,
    ),
)

CASES_BY_ID = {c.case_id: c for c in CASES}
CATEGORIES = tuple(sorted({c.category for c in CASES}))
