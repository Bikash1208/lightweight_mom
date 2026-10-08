"""Monte Carlo error bar on the lens prediction: propagate the local-periodicity error (Phase 19 Part D-bis).

The prediction in ``pipeline.predict_pattern`` uses each cell's *isolated* response.  Exact supercell RCWA shows the cell's
effective field differs from it by a complex perturbation eps_c that grows with the phase step to its neighbours (see
``scripts/research/phase19_coupling_samples.py``).  Here eps_c is drawn per cell from the measured relation, applied as
E_c' = E_c (1 + eps_c), and the pattern is recomputed.  The far field is linear in E, so all draws go through one
precomputed operator (chunked) -- thousands of draws cost seconds.

Metrics per draw: both beam peaks (gain, direction), inter-beam level difference, SLL (max co-pol outside 6 deg cones around
the NOMINAL beam directions, relative to the lower beam peak -- the definition of ``pattern.beam_metrics``).
"""

from __future__ import annotations

import math

import numpy as np


def neighbour_lists(x_mm, y_mm, pitch_mm):
    """For every cell the indices of its 4-neighbours: list of (index, axis) with axis 'x' or 'y' (E_x polarisation)."""
    idx = {
        (math.floor(x / pitch_mm), math.floor(y / pitch_mm)): i for i, (x, y) in enumerate(zip(x_mm, y_mm, strict=True))
    }
    out = []
    for x, y in zip(x_mm, y_mm, strict=True):
        ix, iy = math.floor(x / pitch_mm), math.floor(y / pitch_mm)
        nb = []
        for dx, dy, axis in ((1, 0, "x"), (-1, 0, "x"), (0, 1, "y"), (0, -1, "y")):
            j = idx.get((ix + dx, iy + dy))
            if j is not None:
                nb.append((j, axis))
        out.append(nb)
    return out


def phase_steps_deg(t_co, neighbours):
    """|wrapped phase difference| (deg) of isolated co-pol responses for each (cell, neighbour) edge; list of arrays."""
    ph = np.degrees(np.angle(t_co))
    return [
        np.array([abs((ph[j] - ph[i] + 180.0) % 360.0 - 180.0) for j, _ax in nb]) for i, nb in enumerate(neighbours)
    ]


class PatternOperator:
    """Linear map from cell fields to co-pol U on (a) two fine beam grids and (b) a global grid, chunked."""

    def __init__(self, x_m, y_m, k0, pitch_m, u_inc, beam_deg, theta_step=0.25, phi_step=2.0, chunk=1500):
        self.x, self.y, self.k0, self.p, self.uinc, self.chunk = x_m, y_m, k0, pitch_m, u_inc, chunk
        self.beam_deg = beam_deg
        u0 = math.sin(math.radians(beam_deg))
        ux = u0 + np.arange(-0.06, 0.0601, 0.001)
        uy = np.arange(-0.06, 0.0601, 0.004)
        gx, gy = np.meshgrid(ux, uy, indexing="ij")
        self.beam_u = [(gx.ravel() * s, gy.ravel()) for s in (+1, -1)]
        th = np.arange(0.0, 90.0001, theta_step)
        ph = np.arange(0.0, 360.0, phi_step)
        T, Ph = np.meshgrid(th, ph, indexing="ij")
        self.gth, self.gph = T.ravel(), Ph.ravel()
        self.gux = np.sin(np.radians(self.gth)) * np.cos(np.radians(self.gph))
        self.guy = np.sin(np.radians(self.gth)) * np.sin(np.radians(self.gph))
        # outside-cone mask relative to the nominal beams (6 deg cone, angle between direction vectors)
        cosb = math.cos(math.radians(beam_deg))
        ct = np.cos(np.radians(self.gth))
        mask = np.ones(self.gth.size, bool)
        for s in (+1, -1):
            ang = np.degrees(np.arccos(np.clip(self.gux * s * u0 + ct * cosb, -1, 1)))
            mask &= ang > 6.0
        self.sll_mask = mask

    def _matrix(self, ux, uy):
        k0, p = self.k0, self.p
        ph = np.exp(1j * k0 * (np.outer(ux, self.x) + np.outer(uy, self.y)))
        patch = (
            np.sinc(k0 * (ux[:, None] - self.uinc[0][None, :]) * p / 2 / math.pi)
            * np.sinc(k0 * (uy[:, None] - self.uinc[1][None, :]) * p / 2 / math.pi)
            * p**2
        )
        return ph * patch

    @staticmethod
    def _kfac(ux, uy, k0):
        ct = np.sqrt(np.clip(1 - ux**2 - uy**2, 0, 1))
        return (k0 / (2 * math.pi)) ** 2 * ((1 + ct) / 2) ** 2

    def run(self, ex):
        """ex: (n_cells, n_draws) complex co-pol fields.  Returns dict of per-draw arrays."""
        nd = ex.shape[1]
        peak, pos = [], []
        for ux, uy in self.beam_u:
            f = self._matrix(ux, uy) @ ex
            u = self._kfac(ux, uy, self.k0)[:, None] * np.abs(f) ** 2
            i = np.argmax(u, axis=0)
            peak.append(u[i, np.arange(nd)])
            pos.append(np.stack([ux[i], uy[i]], axis=1))
        sll_peak = np.zeros(nd)
        for s in range(0, self.gth.size, self.chunk):
            sl = slice(s, s + self.chunk)
            sel = self.sll_mask[sl]
            if not sel.any():
                continue
            ux, uy = self.gux[sl][sel], self.guy[sl][sel]
            f = self._matrix(ux, uy) @ ex
            u = self._kfac(ux, uy, self.k0)[:, None] * np.abs(f) ** 2
            sll_peak = np.maximum(sll_peak, u.max(axis=0))
        g = [10 * np.log10(4 * math.pi * pk) for pk in peak]
        theta = [np.degrees(np.arcsin(np.hypot(p[:, 0], p[:, 1]))) for p in pos]
        phi = [np.degrees(np.arctan2(p[:, 1], p[:, 0])) for p in pos]
        nominal = [(self.beam_deg, 0.0), (self.beam_deg, 180.0)]
        point_err = []
        for p, (t0, p0) in zip(pos, nominal, strict=True):
            tt = np.arcsin(np.clip(np.hypot(p[:, 0], p[:, 1]), 0, 1))
            pp = np.arctan2(p[:, 1], p[:, 0])
            t0r, p0r = math.radians(t0), math.radians(p0)
            cosang = np.cos(tt) * math.cos(t0r) + np.sin(tt) * math.sin(t0r) * np.cos(pp - p0r)
            point_err.append(np.degrees(np.arccos(np.clip(cosang, -1, 1))))
        lower = np.minimum(g[0], g[1])
        sll_db = 10 * np.log10(4 * math.pi * sll_peak) - lower
        return dict(
            gain_dbi=np.stack(g, axis=1),
            theta_deg=np.stack(theta, axis=1),
            phi_deg=np.stack(phi, axis=1),
            pointing_error_deg=np.stack(point_err, axis=1),
            level_diff_db=g[0] - g[1],
            sll_db=sll_db,
            mean_gain_dbi=0.5 * (g[0] + g[1]),
        )


def percentiles(a, qs=(10, 50, 90)):
    return {f"p{q}": float(np.percentile(a, q)) for q in qs}
