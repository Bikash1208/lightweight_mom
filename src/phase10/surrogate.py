"""Sec. 22-24, 27, 35-37: Baseline 3, "the 5-cell/9-cell nonlinear surrogate" -- the simplest
model that extends Baseline 2 (Sec. 23) to capture the nonlinearity Sec. 11 measures: a per-center
baseline plus a RIDGE-REGULARIZED polynomial (quadratic self- and opposite-pair cross-terms)
correction in the 8 neighbors' (dL, dW) offsets from the center, fit separately for Re(S)/Im(S)
(Rule 5 -- cartesian, never phase-alone). Plain least squares/ridge regression via numpy only: no
new dependency, and "the simplest model that meets the validation requirements should be
preferred" (Sec. 22) -- a small neural network or GP was not justified by anything measured so far
in this project.
"""

import json
from dataclasses import dataclass

import numpy as np

from . import environment_geometry as eg
from .dataset import ALL_POSITIONS
from .isolated_reference import IsolatedReferenceCache

RIDGE_CANDIDATES = (0.0, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 1.0)


def _linear_deltas(sample):
    cl, cw = sample.center_lw_mm
    out = []
    for pos in ALL_POSITIONS:
        nl, nw = sample.neighbors_lw_mm[pos]
        out.extend((nl - cl, nw - cw))
    return out


def nonlinear_features(sample):
    """16 linear + 24 self-quadratic (dL^2, dW^2, dL*dW per position) + 8 opposite-pair cross
    terms (N-S, E-W, NE-SW, NW-SE) = 48 features total. Opposite pairs only (not every pairwise
    combination of all 8 positions, which would be 28 unphysically undirected cross terms for a
    dataset of a few hundred samples) -- Sec. 11's own nonlinearity test specifically probes N+S
    and E+W interactions, so those (and their diagonal analogues) are the physically motivated
    cross terms to include."""
    linear = _linear_deltas(sample)
    idx = {pos: k for k, pos in enumerate(ALL_POSITIONS)}

    def dl(pos):
        return linear[2 * idx[pos]]

    def dw(pos):
        return linear[2 * idx[pos] + 1]

    quad = []
    for pos in ALL_POSITIONS:
        d_l, d_w = dl(pos), dw(pos)
        quad.extend((d_l * d_l, d_w * d_w, d_l * d_w))

    cross = []
    for a, b in (("N", "S"), ("E", "W"), ("NE", "SW"), ("NW", "SE")):
        cross.extend((dl(a) * dl(b), dw(a) * dw(b)))

    return np.asarray(linear + quad + cross, dtype=float)


N_FEATURES = 16 + 24 + 8


@dataclass(frozen=True)
class DomainSpec:
    l_min_mm: float = eg.DOMAIN_L_MIN_MM
    l_max_mm: float = eg.DOMAIN_L_MAX_MM
    w_min_mm: float = eg.DOMAIN_W_MIN_MM
    w_max_mm: float = eg.DOMAIN_W_MAX_MM
    frequency_hz: float = 12.5e9
    neighborhood: str = "9_cell"
    polarization: str = "x"

    def contains(self, sample):
        if abs(sample.frequency_hz - self.frequency_hz) > 1.0:
            return False
        l, w = sample.center_lw_mm
        if not (self.l_min_mm - 1e-9 <= l <= self.l_max_mm + 1e-9 and
                self.w_min_mm - 1e-9 <= w <= self.w_max_mm + 1e-9):
            return False
        for pos in ALL_POSITIONS:
            nl, nw = sample.neighbors_lw_mm[pos]
            if not (self.l_min_mm - 1e-9 <= nl <= self.l_max_mm + 1e-9 and
                    self.w_min_mm - 1e-9 <= nw <= self.w_max_mm + 1e-9):
                return False
        return True


OUT_OF_DOMAIN = "OUT_OF_VALIDATED_DOMAIN"


@dataclass(frozen=True)
class EnvironmentSurrogate:
    """Sec. 27 point 7: ``domain`` is machine-readable (``DomainSpec``, also exported to JSON by
    ``to_dict``). Sec. 27 point 6: ``predict`` never silently extrapolates -- an out-of-domain
    input returns the literal string ``OUT_OF_VALIDATED_DOMAIN`` instead of a number. ``S_C``
    (the zero-neighbor-perturbation reference) is looked up from ``isolated_cache`` -- a real,
    on-demand MoM solve (see ``isolated_reference.py``), not a training-data average, so it is
    exact and available for ANY center geometry in the validated continuous range, seen in
    training or not (Sec. 21)."""
    isolated_cache: IsolatedReferenceCache
    beta_re: np.ndarray
    beta_im: np.ndarray
    ridge_lambda: float
    domain: DomainSpec

    def predict(self, sample):
        if not self.domain.contains(sample):
            return OUT_OF_DOMAIN
        s0 = self.isolated_cache.get(sample.center_lw_mm, sample.frequency_hz)
        x = nonlinear_features(sample)
        return complex(s0.real + float(x @ self.beta_re), s0.imag + float(x @ self.beta_im))

    def to_dict(self):
        cache = self.isolated_cache
        return {
            "isolated_reference_config": {
                "period_m": cache.period_m, "mesh_cells_per_subcell": cache.mesh_cells_per_subcell,
                "substrate": {"name": cache.substrate.name, "eps_r": cache.substrate.eps_r,
                              "tan_delta": cache.substrate.tan_delta, "thickness_m": cache.substrate.thickness},
                "metal_model": cache.metal.model,
            },
            "beta_re": list(self.beta_re), "beta_im": list(self.beta_im),
            "ridge_lambda": self.ridge_lambda,
            "domain": {
                "l_min_mm": self.domain.l_min_mm, "l_max_mm": self.domain.l_max_mm,
                "w_min_mm": self.domain.w_min_mm, "w_max_mm": self.domain.w_max_mm,
                "frequency_hz": self.domain.frequency_hz, "neighborhood": self.domain.neighborhood,
                "polarization": self.domain.polarization,
            },
            "n_features": N_FEATURES,
            "feature_order": "16 linear (per ALL_POSITIONS dL,dW) + 24 self-quadratic + 8 opposite-pair cross",
        }

    @classmethod
    def from_dict(cls, doc):
        from ..materials import Metal, Substrate

        c = doc["isolated_reference_config"]
        sub = c["substrate"]
        substrate = Substrate(sub["name"], sub["eps_r"], sub["tan_delta"], sub["thickness_m"])  # positional: (name, eps_r, tan_delta, thickness)
        metal = Metal(c["metal_model"])
        cache = IsolatedReferenceCache(c["period_m"], substrate, metal, c["mesh_cells_per_subcell"])
        d = doc["domain"]
        domain = DomainSpec(d["l_min_mm"], d["l_max_mm"], d["w_min_mm"], d["w_max_mm"],
                             d["frequency_hz"], d["neighborhood"], d["polarization"])
        return cls(
            cache, np.asarray(doc["beta_re"], dtype=float), np.asarray(doc["beta_im"], dtype=float),
            doc["ridge_lambda"], domain,
        )

    def save(self, path):
        from pathlib import Path
        Path(path).write_text(json.dumps(self.to_dict(), indent=1) + "\n")

    @classmethod
    def load(cls, path):
        from pathlib import Path
        return cls.from_dict(json.loads(Path(path).read_text()))


def _ridge_fit(x, y, lam):
    n_features = x.shape[1]
    a = x.T @ x + lam * np.eye(n_features)
    b = x.T @ y
    return np.linalg.solve(a, b)


def _fit_for_lambda(train_samples, lam, frequency_hz, isolated_cache):
    x_rows, y_re, y_im = [], [], []
    for s in train_samples:
        x_rows.append(nonlinear_features(s))
        s0 = isolated_cache.get(s.center_lw_mm, s.frequency_hz)
        residual = s.s_complex - s0
        y_re.append(residual.real)
        y_im.append(residual.imag)
    x = np.asarray(x_rows, dtype=float)
    beta_re = _ridge_fit(x, np.asarray(y_re), lam)
    beta_im = _ridge_fit(x, np.asarray(y_im), lam)
    domain = DomainSpec(frequency_hz=frequency_hz)
    return EnvironmentSurrogate(isolated_cache, beta_re, beta_im, lam, domain)


def fit_surrogate(train_samples, validation_samples, frequency_hz, isolated_cache,
                   ridge_candidates=RIDGE_CANDIDATES):
    """Sec. 20: the ridge strength is selected on the VALIDATION split (minimizing mean complex
    error there), never on the blind test set -- hyperparameter selection is exactly what the
    validation split exists for."""
    from .error_metrics import compute_error_distribution

    best_model, best_lambda, best_val_error = None, None, float("inf")
    for lam in ridge_candidates:
        model = _fit_for_lambda(train_samples, lam, frequency_hz, isolated_cache)
        preds = [model.predict(s) for s in validation_samples]
        preds = [p if isinstance(p, complex) else complex("nan") for p in preds]
        refs = [s.s_complex for s in validation_samples]
        report = compute_error_distribution(preds, refs)
        if report.mean_complex_error < best_val_error:
            best_val_error, best_model, best_lambda = report.mean_complex_error, model, lam
    return best_model, best_lambda, best_val_error
