"""Part 3, HMSA item 2: a MoM surface-wave-eigenmode capability -- **NOT IMPLEMENTED**.

PHASE13_WEBAPP_PROMPT.md's explicit fallback: "If it can't be finished cleanly in this phase, ship
(1) [`hfss_table.py`] and leave (2) as a documented, tested scaffold -- don't block the app on it."
This module is that scaffold: the physics, the plan for a real implementation, and a function that
raises clearly (never silently returns a fabricated number) so the app can offer the menu item,
explain why it's unavailable, and point at `hfss_table.py` as the working alternative.

**What a real implementation needs** (for the record, not attempted here):
A surface-wave eigenmode of the periodic grounded unit cell is a solution of this project's own
Z(kx, ky) I = 0 (source-free) system at a REAL (kx, ky) with kt = sqrt(kx^2+ky^2) > k0 (bound,
non-radiating -- `CONVENTIONS.md`'s `Im(kz) <= 0` branch choice already fixes the right analytic
continuation for kt > k0, since `Im(kz)` there is already negative/decaying by the existing branch
rule; no new sign convention needed). This is generally NOT the same (kx0, ky0) = (0,0) normal-
incidence excitation point `sparameters.compute_s_parameters` solves at -- it requires:

1. Assembling `mom_matrix.assemble_impedance_matrix(basis, green, kx, ky, floquet_order)` (already
   supports non-zero (kx, ky) -- reusable unmodified) at a family of trial (kx, ky) along the
   direction of interest.
2. Locating where `det(Z(kx, ky))` has a zero (or, in practice, an extremum of the smallest
   singular value of Z, more numerically stable than a raw determinant for a matrix this size) --
   a genuine root-finding/optimization problem over one or two real parameters, needing a stable
   bracketing strategy this project has not yet developed or validated.
3. Validating the found root against `dispersion.analytic_tm0` (the bare-slab case, patch pattern
   set to fully absent) BEFORE trusting it on a real patterned cell -- required by the brief before
   this label may read anything but `EXPERIMENTAL`.

None of this exists yet. `find_surface_wave_mode` below raises immediately so nothing downstream
can silently receive a fabricated dispersion point.
"""

from ..errors import LwmomError, ErrorCode
from ..labels import DISPERSION_FROM_MOM, EXPERIMENTAL

NOT_IMPLEMENTED = ErrorCode(
    "DISPERSION_MOM_NOT_IMPLEMENTED",
    "MoM surface-wave eigenmode search is not implemented in this phase.",
    "Use an imported HFSS eigenmode table instead (dispersion.hfss_table.read_dispersion_csv, "
    "label DISPERSION_FROM_HFSS).",
)

SCAFFOLD_LABELS = (DISPERSION_FROM_MOM, EXPERIMENTAL)


def find_surface_wave_mode(unit_cell, frequency_hz, direction_deg=0.0):
    """Raises `LwmomError` (code `DISPERSION_MOM_NOT_IMPLEMENTED`) unconditionally. Kept as a
    real, importable, tested function (not just a comment) so the app can present this menu item,
    catch the specific error code, and show the remedy -- rather than the feature silently not
    existing at all."""
    raise LwmomError(NOT_IMPLEMENTED)
