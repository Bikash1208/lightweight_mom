"""Sec. 23: the three required baselines, evaluated before any more complex surrogate is
considered (Sec. 22's "the simplest model that meets the validation requirements should be
preferred").
"""

from dataclasses import dataclass

from .isolated_reference import IsolatedReferenceCache

ALL_POSITIONS = ("N", "E", "S", "W", "NE", "NW", "SE", "SW")


def build_baseline_1(isolated_cache):
    """S = S(G_C): the OLD, now-invalidated (Phase 9) periodic-cell lookup, ignoring every
    neighbor. Quantifies exactly how much error remains from the LPA approach Phase 9 invalidated.
    Generalizes to any center geometry automatically (a real MoM solve, not a training lookup)."""
    return lambda sample: isolated_cache.get(sample.center_lw_mm, sample.frequency_hz)


def _delta_features(sample):
    """[dL_N, dW_N, dL_E, dW_E, ..., dL_SW, dW_SW] -- 16 values, each neighbor's offset from the
    center geometry. No intercept per position: DeltaS_i(0, 0) = 0 by construction (an unperturbed
    neighbor contributes nothing), which is the physical content of "additive correction"."""
    cl, cw = sample.center_lw_mm
    feats = []
    for pos in ALL_POSITIONS:
        nl, nw = sample.neighbors_lw_mm[pos]
        feats.extend((nl - cl, nw - cw))
    return feats


@dataclass(frozen=True)
class Baseline2AdditiveCorrection:
    """S = S_C + sum_i DeltaS_i(G_i), each DeltaS_i linear in the neighbor's (dL, dW) offset from
    the center -- fit as ONE joint least-squares regression (Re/Im separately, Rule 5) over all 16
    neighbor-offset features at once (not 8 independent marginal regressions, which would let each
    position's fit silently absorb the other positions' simultaneous effects). This is exactly the
    additivity assumption Sec. 11's nonlinearity test measures the (in)validity of; fitting it
    explicitly and scoring it on the same blind test as every other model quantifies that directly.
    ``S_C`` itself is looked up from ``isolated_cache`` (a real MoM solve, memoized), not averaged
    from training samples -- it is a physical quantity available for ANY center geometry, including
    ones never seen during training (Sec. 21's whole point)."""

    isolated_cache: IsolatedReferenceCache
    beta_re: object  # length-16 ndarray
    beta_im: object

    def predict(self, sample):
        s = self.isolated_cache.get(sample.center_lw_mm, sample.frequency_hz)
        feats = _delta_features(sample)
        import numpy as np

        x = np.asarray(feats, dtype=float)
        return complex(s.real + float(x @ self.beta_re), s.imag + float(x @ self.beta_im))


def fit_baseline_2(train_samples, isolated_cache):
    """Sec. 23. The correction coefficients are fit on the training split only (Sec. 20: never
    touches validation/blind data); ``isolated_cache`` supplies S_C for every sample (training and,
    later, at prediction time, any other sample) from the real solver, on demand."""
    import numpy as np

    x_rows, y_re, y_im = [], [], []
    for s in train_samples:
        x_rows.append(_delta_features(s))
        s0 = isolated_cache.get(s.center_lw_mm, s.frequency_hz)
        residual = s.s_complex - s0
        y_re.append(residual.real)
        y_im.append(residual.imag)
    x = np.asarray(x_rows, dtype=float)
    beta_re, *_ = np.linalg.lstsq(x, np.asarray(y_re), rcond=None)
    beta_im, *_ = np.linalg.lstsq(x, np.asarray(y_im), rcond=None)
    return Baseline2AdditiveCorrection(isolated_cache, beta_re, beta_im)
