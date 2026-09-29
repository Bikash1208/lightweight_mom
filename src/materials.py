"""Material data structures.

Nothing here has default electrical values: eps_r, tan_delta and thickness must
come from configuration.

Phase 1 assumptions (also listed in README):
* non-magnetic (mu_r = 1),
* eps_r and tan_delta are constant over frequency (no dispersion model),
* the metal is PEC.
"""

from dataclasses import dataclass

from ._validate import finite_real
from .frequency import free_space_wavenumber


@dataclass(frozen=True)
class Substrate:
    """Single homogeneous, isotropic dielectric layer.

    Complex relative permittivity under exp(+j*omega*t):

        eps_c = eps_r * (1 - j*tan_delta)

    so Im(eps_c) <= 0 for a passive (lossy) medium.
    """

    name: str
    eps_r: float  # real relative permittivity, > 0
    tan_delta: float  # loss tangent, >= 0
    thickness: float  # [m], > 0

    def __post_init__(self):
        for name in ("eps_r", "tan_delta", "thickness"):
            object.__setattr__(self, name, finite_real(name, getattr(self, name)))
        if self.eps_r <= 0.0:
            raise ValueError(f"eps_r must be > 0, got {self.eps_r}")
        if self.tan_delta < 0.0:
            raise ValueError(f"tan_delta must be >= 0 (passive medium), got {self.tan_delta}")
        if self.thickness <= 0.0:
            raise ValueError(f"thickness must be > 0 m, got {self.thickness}")

    @property
    def eps_complex(self):
        """Complex relative permittivity eps_r*(1 - j*tan_delta)."""
        return self.eps_r * (1.0 - 1j * self.tan_delta)

    @property
    def refractive_index(self):
        """sqrt(eps_c), principal branch => Re > 0 and Im <= 0."""
        return complex(self.eps_complex) ** 0.5

    def wavenumber(self, f_hz):
        """k = k0 * sqrt(eps_c)  [rad/m], complex with Im(k) <= 0."""
        return free_space_wavenumber(f_hz) * self.refractive_index


@dataclass(frozen=True)
class Metal:
    """Metallic pattern material. Phase 1 supports only ``pec``.

    Finite conductivity (e.g. copper surface impedance) is deliberately not
    modelled yet; requesting it is an error rather than being silently treated
    as PEC.
    """

    model: str = "pec"

    def __post_init__(self):
        if self.model != "pec":
            raise ValueError(f"metal model {self.model!r} is not supported in Phase 1; only 'pec'")
