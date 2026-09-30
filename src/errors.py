"""Phase 13 Part 2: a typed error-code taxonomy for the app layer (and any script that wants it).
The solver core keeps raising its own exceptions (`ConfigError`, `ValueError`, ...) exactly as
before -- unchanged, per CLAUDE.md's "physics/numerics untouched" rule; `LwmomError` and its
subclasses below are a thin, additive WRAPPING layer the app uses at its own boundary (job
runner / API handlers) to attach a stable machine-readable `code`, a short user-facing `message`,
a `remedy` (what to do about it) and optional `debug_detail` (the original exception's own text) --
"no silent fallbacks": every one of these is raised, never swallowed into a default value.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class ErrorCode:
    code: str
    default_message: str
    default_remedy: str


CONFIG_INVALID = ErrorCode(
    "CONFIG_INVALID",
    "The configuration is invalid or incomplete.",
    "Check the reported field against docs/theory/CONVENTIONS.md; no electrical value has a code default.",
)
GEOMETRY_OUT_OF_CELL = ErrorCode(
    "GEOMETRY_OUT_OF_CELL",
    "The requested geometry does not fit inside the unit cell.",
    "Reduce the patch/feature size or increase the period.",
)
MESH_UNALIGNED = ErrorCode(
    "MESH_UNALIGNED",
    "A feature edge does not sit on a mesh line at this resolution.",
    "Increase the mesh resolution or adjust the dimension to a mesh-aligned value "
    "(see `src/discretization.py::feature_edges_aligned`); the solve still runs, but with staircasing error.",
)
FLOQUET_RANK_DEFICIENT = ErrorCode(
    "FLOQUET_RANK_DEFICIENT",
    "The Floquet truncation order is too small for this mesh.",
    "Increase `floquet_order` so that 2*floquet_order + 1 >= max(cells_x, cells_y).",
)
SOLVE_FAILED = ErrorCode(
    "SOLVE_FAILED",
    "The linear solve did not complete.",
    "Check the impedance matrix condition (see the run's SolveDiagnostics); an extremely fine mesh "
    "or degenerate geometry can make Z singular or numerically unstable.",
)
NOT_CONVERGED = ErrorCode(
    "NOT_CONVERGED",
    "The mesh-convergence (Richardson) fit did not converge.",
    "Report the finest fixed-mesh result with a MESH_NOT_CONVERGED label instead of an "
    "extrapolated one; see `src/phase13/convergence.py` for the fit's own failure modes.",
)
PHASE_UNREACHABLE = ErrorCode(
    "PHASE_UNREACHABLE",
    "The requested reflection phase is outside this library's achievable range.",
    "Widen the geometry-parameter sweep, accept a quantization error, or use a different library.",
)
CAD_UNAVAILABLE = ErrorCode(
    "CAD_UNAVAILABLE",
    "STEP export / CAD-dependent features need the optional 'cad' extra.",
    "Install with `pip install -e '.[cad]'` (cadquery, matplotlib).",
)
FIELD_EXPRESSION_UNSUPPORTED = ErrorCode(
    "FIELD_EXPRESSION_UNSUPPORTED",
    "This reference-field file's field-calculator expression is not recognised.",
    "Supported: CmplxPhase/CmplxMag/CmplxReal/CmplxImag of ScalarX/Y/Z(<Ex,Ey,Ez>). "
    "Re-export with one of these, or extend `src/fields/fld_parser.py`'s expression table.",
)
IMPEDANCE_UNREACHABLE = ErrorCode(
    "IMPEDANCE_UNREACHABLE",
    "A required surface reactance value is outside the dispersion library's range.",
    "Widen the dispersion library's geometry sweep, or accept the nearest achievable value "
    "(the app reports which cells this affects).",
)
LIBRARY_FREQUENCY_UNAVAILABLE = ErrorCode(
    "LIBRARY_FREQUENCY_UNAVAILABLE",
    "The requested frequency has no stored record in this library.",
    "Choose one of the library's own stored frequencies, or build/select a library that includes "
    "the requested one -- `UnitCellLibrary.response` never interpolates across frequency.",
)
DISPERSION_MOM_NO_CLEAR_MODE = ErrorCode(
    "DISPERSION_MOM_NO_CLEAR_MODE",
    "The MoM surface-wave search found no clear bound mode for this cell/frequency.",
    "Try a different geometry_param, direction, or mesh; or import an HFSS eigenmode table instead "
    "(dispersion.hfss_table.read_dispersion_csv, label DISPERSION_FROM_HFSS).",
)
DISPERSION_SHEET_MODEL_NO_ROOT = ErrorCode(
    "DISPERSION_SHEET_MODEL_NO_ROOT",
    "The homogenized sheet-admittance transverse-resonance equation has no bound-mode root for this cell/frequency.",
    "This patch size may be too close to the cell boundary (mesh-limited) or past this model's usable "
    "range; try a different geometry_param, a finer mesh, or import an HFSS eigenmode table instead.",
)

ALL_ERROR_CODES = (
    CONFIG_INVALID,
    GEOMETRY_OUT_OF_CELL,
    MESH_UNALIGNED,
    FLOQUET_RANK_DEFICIENT,
    SOLVE_FAILED,
    NOT_CONVERGED,
    PHASE_UNREACHABLE,
    CAD_UNAVAILABLE,
    FIELD_EXPRESSION_UNSUPPORTED,
    IMPEDANCE_UNREACHABLE,
    LIBRARY_FREQUENCY_UNAVAILABLE,
    DISPERSION_MOM_NO_CLEAR_MODE,
    DISPERSION_SHEET_MODEL_NO_ROOT,
)


class LwmomError(Exception):
    """Base class for every app-boundary error. `code` is one of `ALL_ERROR_CODES`'s `.code`
    strings (stable, safe to match on in the UI/tests); `message`/`remedy` default to the code's
    own text but can be overridden with case-specific detail; `debug_detail` carries the original
    exception's `str()` when this wraps one (never silently discarded)."""

    def __init__(self, error_code: ErrorCode, message=None, remedy=None, debug_detail=None):
        self.error_code = error_code
        self.code = error_code.code
        self.message = message or error_code.default_message
        self.remedy = remedy or error_code.default_remedy
        self.debug_detail = debug_detail
        super().__init__(f"[{self.code}] {self.message}")

    def to_dict(self):
        return {"code": self.code, "message": self.message, "remedy": self.remedy, "debug_detail": self.debug_detail}

    @classmethod
    def wrap(cls, error_code: ErrorCode, original_exception, message=None, remedy=None):
        """Wraps an already-raised exception (from the solver core or elsewhere), preserving its
        text as `debug_detail` -- the standard way app code turns a low-level exception into one
        of this taxonomy's stable codes without losing the original diagnostic information."""
        return cls(error_code, message=message, remedy=remedy, debug_detail=str(original_exception))


def config_invalid(original_exception=None, message=None):
    if original_exception is not None:
        return LwmomError.wrap(CONFIG_INVALID, original_exception, message=message)
    return LwmomError(CONFIG_INVALID, message=message)


def geometry_out_of_cell(original_exception=None, message=None):
    if original_exception is not None:
        return LwmomError.wrap(GEOMETRY_OUT_OF_CELL, original_exception, message=message)
    return LwmomError(GEOMETRY_OUT_OF_CELL, message=message)


def mesh_unaligned(message=None):
    return LwmomError(MESH_UNALIGNED, message=message)


def floquet_rank_deficient(message=None):
    return LwmomError(FLOQUET_RANK_DEFICIENT, message=message)


def solve_failed(original_exception=None, message=None):
    if original_exception is not None:
        return LwmomError.wrap(SOLVE_FAILED, original_exception, message=message)
    return LwmomError(SOLVE_FAILED, message=message)


def not_converged(message=None):
    return LwmomError(NOT_CONVERGED, message=message)


def phase_unreachable(message=None):
    return LwmomError(PHASE_UNREACHABLE, message=message)


def cad_unavailable(original_exception=None, message=None):
    if original_exception is not None:
        return LwmomError.wrap(CAD_UNAVAILABLE, original_exception, message=message)
    return LwmomError(CAD_UNAVAILABLE, message=message)


def field_expression_unsupported(expression, message=None):
    return LwmomError(FIELD_EXPRESSION_UNSUPPORTED, message=message or f"Unsupported field expression: {expression!r}")


def impedance_unreachable(message=None):
    return LwmomError(IMPEDANCE_UNREACHABLE, message=message)


def library_frequency_unavailable(message=None):
    return LwmomError(LIBRARY_FREQUENCY_UNAVAILABLE, message=message)


def dispersion_mom_no_clear_mode(message=None):
    return LwmomError(DISPERSION_MOM_NO_CLEAR_MODE, message=message)


def dispersion_sheet_model_no_root(message=None):
    return LwmomError(DISPERSION_SHEET_MODEL_NO_ROOT, message=message)
