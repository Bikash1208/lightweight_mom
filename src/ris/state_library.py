"""CP-2 state library (Phase 20 R4): tabulated MoM reflection of the four one-hot states over (f, theta, phi), a mesh-correction map, and the file format.

Storage frame. The raw MoM output is the 2x2 power-wave reflection in the incidence-plane (TE, TM) frame (`ObliqueSParameters.s11`). The library stores the
circular reflection `C[a, b]` (a = reflected, b = incident handedness, propagation-relative, `docs/theory/CONVENTIONS.md`) referred to the FIXED structure axes
(`ris.polarization.to_structure_frame`): then the Pancharatnam-Berry phase of a turned element is explicit and the dependence on the incidence azimuth is weak
(smooth to interpolate; the mesh correction is nearly phi independent, R4 pilot). `S` in the (TE, TM) frame is recovered exactly (`te_tm_from_circular`).

Mesh correction (no surrogate, no fit to data of another kind): at a sparse correction set the same points are solved at meshes 60 (the grid mesh), 90 and 120, ALL with the
Floquet cutoff F = 2n (a smaller cutoff leaves isolated numerical lines, see `results/phase20_r4/FLOQUET_CONVERGENCE.md`). With a converged cutoff the mesh convergence is regular
(close to second order), so the converged value is the Richardson extrapolation `C_inf = C120 + R2 (C120 - C90)` with `R2 = 120^-2 / (90^-2 - 120^-2)`; the correction is
`dC = C_inf - C60`, averaged over the two azimuths of the set and interpolated in (f, theta), then added to every grid value. The uncertainty of every library entry is the distance between
the p = 1 and p = 2 extrapolations, `|C_inf(p=1) - C_inf(p=2)|` (R2/R3 practice), plus the azimuth spread of `dC`; both are stored. The model uncertainty (staircased octagon, L_s = 0 assumed,
datasheet extrapolation) is NOT in these numbers.

Arrays (`cp2_library_v1.npz`; axes `f_ghz` (41), `theta_deg` (10), `phi_deg` (24), `state` (4)):
  `c_raw`   (f, theta, phi, state, 2, 2)  structure-frame circular reflection at mesh 60 (the grid mesh)
  `c_corr`  same shape                    `c_raw + correction`  (the library values)
  `correction`, `uncertainty`             same shape, complex / real
  `s_corr`  (f, theta, phi, state, 2, 2)  the corrected reflection in the incidence-plane (TE, TM) frame
"""

import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .polarization import (
    cartesian_reflection,
    circular_reflection,
    from_structure_frame,
    jones_to_te_tm,
    te_tm_from_circular,
    to_structure_frame,
)

N_STATES = 4
PERIOD_MM = 15.0
C0_MM_GHZ = 299.792458  # speed of light in mm * GHz
WOOD_GUARD = 1.03  # flag |k_t + G| / k0 below this: a Wood-anomaly dip sits just below onset (R4 grid QA: 1.0158 and 1.0227 seen)
R2_RICHARDSON = (120.0**-2) / (90.0**-2 - 120.0**-2)  # p = 2 extrapolation weight on (C120 - C90)
R1_RICHARDSON = (120.0**-1) / (90.0**-1 - 120.0**-1)  # p = 1 (= 3)
THETA_LABEL_LIMIT_DEG = 30.0  # beyond the R3 targets: EXTRAPOLATED_BEYOND_R3_TARGETS


def higher_order_margin(f_ghz, theta_deg, phi_deg, period_mm=PERIOD_MM):
    """min over the first Floquet orders (m, n) != (0, 0) of |k_t + G_mn| / k0 for the lattice of `period_mm`. A value < 1 means a higher diffraction order
    PROPAGATES (the stored (0,0) reflection is then not the whole response and absorbed power cannot be read from it); < `WOOD_GUARD` is the guard band below onset."""
    f, th, ph = np.broadcast_arrays(*(np.asarray(x, dtype=float) for x in (f_ghz, theta_deg, phi_deg)))
    a = C0_MM_GHZ / (f * period_mm)  # lambda / P
    ux, uy = np.sin(np.radians(th)) * np.cos(np.radians(ph)), np.sin(np.radians(th)) * np.sin(np.radians(ph))
    out = np.full(f.shape, np.inf)
    for m in (-1, 0, 1):
        for n in (-1, 0, 1):
            if (m, n) != (0, 0):
                out = np.minimum(out, np.hypot(ux + a * m, uy + a * n))
    return out


def fill_normal_incidence(s_grid, theta_deg, phi_deg):
    """`s_grid[f, theta, phi, state, 2, 2]` with only the phi = 0 entry of theta = 0 filled: the other phi follow by the exact frame rotation (at theta = 0 the
    azimuth only relabels the TE/TM axes). Returns a new array."""
    out = np.array(s_grid, dtype=complex)
    i0 = int(np.argmin(np.abs(np.asarray(theta_deg))))
    if abs(theta_deg[i0]) > 1e-9:
        raise ValueError("theta grid must contain 0")
    for j, phi in enumerate(phi_deg):
        if abs(phi) < 1e-12:
            continue
        for f in range(out.shape[0]):
            for k in range(out.shape[3]):
                out[f, i0, j, k] = jones_to_te_tm(cartesian_reflection(out[f, i0, 0, k], 0.0), np.radians(phi))
    return out


def structure_frame_circular(s_grid, phi_deg):
    """TE/TM -> structure-frame circular, same leading shape."""
    out = np.zeros_like(s_grid, dtype=complex)
    for idx in np.ndindex(s_grid.shape[:-2]):
        out[idx] = to_structure_frame(circular_reflection(s_grid[idx]), np.radians(phi_deg[idx[2]]))
    return out


def te_tm_from_structure_circular(c_grid, phi_deg):
    out = np.zeros_like(c_grid, dtype=complex)
    for idx in np.ndindex(c_grid.shape[:-2]):
        out[idx] = te_tm_from_circular(from_structure_frame(c_grid[idx], np.radians(phi_deg[idx[2]])))
    return out


def correction_map(f_corr, points, c60, c90, c120):
    """From the correction set. `points` is [(theta, phi)], `cXX[f, point, state, 2, 2]` the structure-frame circular reflection at that mesh. Returns
    `(theta_nodes, dC, spread, mesh_unc)` with shape `[f, theta, state, 2, 2]`: `dC = C_inf(p=2) - C60` (phi-mean), the phi spread (max |dC(phi_i) - dC(phi_j)|), and the mesh
    uncertainty `|C_inf(p=1) - C_inf(p=2)|` (phi-mean), per theta node."""
    pts = np.asarray(points, dtype=float)
    thetas = sorted({float(t) for t in pts[:, 0]})
    c2 = c120 + R2_RICHARDSON * (c120 - c90)
    unc = np.abs((R1_RICHARDSON - R2_RICHARDSON) * (c120 - c90))
    d, sp, st = [], [], []
    for th in thetas:
        sel = np.flatnonzero(np.isclose(pts[:, 0], th))
        dc = (c2 - c60)[:, sel]
        d.append(dc.mean(axis=1))
        st.append(unc[:, sel].mean(axis=1))
        sp.append(
            np.max(np.abs(dc[:, :, None] - dc[:, None, :]), axis=(1, 2))
            if len(sel) > 1
            else np.zeros_like(dc[:, 0], dtype=float)
        )
    mv = lambda a: np.stack(a, axis=1)  # noqa: E731
    return np.array(thetas), mv(d), mv(sp), mv(st)


def interpolate_correction(f_nodes, theta_nodes, values, f_grid, theta_grid, f_method="cubic"):
    """Cubic interpolation of `values[f, theta, state, 2, 2]` (complex or real) onto the library (f, theta) grid, per state and element (re/im separately).
    `f_method="pchip"` (opt-in, Phase 20 R4 v2.1): shape-preserving PCHIP in frequency for each theta node (no overshoot between unevenly spaced nodes, the 0.5 GHz / 0.05 GHz mix
    of the band-refined correction map), then the same theta interpolation as before; the default "cubic" is unchanged."""
    kind = "cubic" if len(theta_nodes) >= 4 else "linear"
    if f_method == "pchip":
        from scipy.interpolate import PchipInterpolator, interp1d

        along_f = np.empty((len(f_grid), len(theta_nodes), *values.shape[2:]), dtype=values.dtype)
        for idx in np.ndindex(values.shape[2:]):
            v = values[(slice(None), slice(None), *idx)]
            parts = (v.real, v.imag) if np.iscomplexobj(v) else (v,)
            res = [PchipInterpolator(np.asarray(f_nodes), p, axis=0, extrapolate=True)(np.asarray(f_grid)) for p in parts]
            along_f[(slice(None), slice(None), *idx)] = res[0] + 1j * res[1] if np.iscomplexobj(v) else res[0]
        out = np.empty((len(f_grid), len(theta_grid), *values.shape[2:]), dtype=values.dtype)
        for idx in np.ndindex(values.shape[2:]):
            v = along_f[(slice(None), slice(None), *idx)]
            parts = (v.real, v.imag) if np.iscomplexobj(v) else (v,)
            res = [interp1d(np.asarray(theta_nodes), p, axis=1, kind=kind, fill_value="extrapolate")(np.asarray(theta_grid)) for p in parts]
            out[(slice(None), slice(None), *idx)] = res[0] + 1j * res[1] if np.iscomplexobj(v) else res[0]
        return out
    if f_method != "cubic":
        raise ValueError("f_method must be 'cubic' or 'pchip'")
    out = np.zeros((len(f_grid), len(theta_grid), *values.shape[2:]), dtype=values.dtype)
    ff, tt = np.meshgrid(f_grid, theta_grid, indexing="ij")
    pts = np.stack([ff.ravel(), tt.ravel()], axis=-1)
    for idx in np.ndindex(values.shape[2:]):
        v = values[(slice(None), slice(None), *idx)]
        parts = (v.real, v.imag) if np.iscomplexobj(v) else (v,)
        res = [
            RegularGridInterpolator((f_nodes, theta_nodes), p, method=kind, bounds_error=False, fill_value=None)(pts)
            for p in parts
        ]
        out[(slice(None), slice(None), *idx)] = (res[0] + 1j * res[1] if np.iscomplexobj(v) else res[0]).reshape(
            ff.shape
        )
    return out


def assemble(f_ghz, theta_deg, phi_deg, s_raw, corr_f, corr_points, c60, c90, c120, f_method="cubic"):
    """All the arrays of the library file from the raw grid (theta = 0 already filled) and the correction set."""
    f_ghz, theta_deg, phi_deg = (np.asarray(a, dtype=float) for a in (f_ghz, theta_deg, phi_deg))
    c_raw = structure_frame_circular(s_raw, phi_deg)
    th_nodes, d, spread, step = correction_map(corr_f, corr_points, c60, c90, c120)
    corr = interpolate_correction(np.asarray(corr_f), th_nodes, d, f_ghz, theta_deg, f_method)  # (f, theta, state, 2, 2)
    unc = interpolate_correction(np.asarray(corr_f), th_nodes, step + spread, f_ghz, theta_deg, f_method)
    corr = np.broadcast_to(
        corr[:, :, None], c_raw.shape
    ).copy()  # same correction at every phi (supported by the pilot)
    unc = np.broadcast_to(np.abs(unc)[:, :, None], c_raw.shape).copy()
    c_corr = c_raw + corr
    return {
        "f_ghz": f_ghz,
        "theta_deg": theta_deg,
        "phi_deg": phi_deg,
        "state": np.arange(1, N_STATES + 1),
        "c_raw": c_raw,
        "c_corr": c_corr,
        "correction": corr,
        "uncertainty": unc,
        "s_corr": te_tm_from_structure_circular(c_corr, phi_deg),
        "s_raw": np.asarray(s_raw, dtype=complex),
        "higher_order_margin": higher_order_margin(
            f_ghz[:, None, None], theta_deg[None, :, None], phi_deg[None, None, :]
        ),
    }


def save_library(path, arrays, manifest):
    path = Path(path)
    np.savez_compressed(path, **arrays)
    path.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8", newline="\n")


# ---------------------------------------------------------------------------------------------------------------------- query
class OutOfDomain(ValueError):
    """The query lies outside the tabulated grid: `NOT_AVAILABLE_FROM_CURRENT_DATA` (no extrapolation)."""

    state = "NOT_AVAILABLE_FROM_CURRENT_DATA"


class StateLibrary:
    """The CP-2 library file with a smooth interpolant in (f, theta, phi).

    Method (justified): the library stores the STRUCTURE-frame circular reflection, a slowly varying complex field of (f, theta, phi) (the Pancharatnam-Berry phase is
    already explicit in that frame, so no unwrapping is needed and entries that pass near zero, like the plate-like terms, are fine). Each of the 4 states x 4 elements
    is interpolated as its real and imaginary parts, tensor-product cubic in (f, theta) and cubic with PERIODIC padding in phi. No value outside the grid is returned.
    """

    def __init__(self, arrays, manifest=None):
        self.a = {k: np.asarray(v) for k, v in arrays.items()}
        self.manifest = manifest or {}
        self.f, self.theta, self.phi = self.a["f_ghz"], self.a["theta_deg"], self.a["phi_deg"]
        self._pad = 3
        self._cache = {}
        step = float(self.phi[1] - self.phi[0])
        self._phi_pad = np.concatenate([self.phi[-self._pad :] - 360.0, self.phi, self.phi[: self._pad] + 360.0])
        assert abs(self.phi[-1] + step - 360.0) < 1e-9, (
            "phi grid must be periodic with the last point one step below 360"
        )

    @classmethod
    def load(cls, path):
        path = Path(path)
        z = np.load(path)
        man = path.with_suffix(".manifest.json")
        return cls({k: z[k] for k in z.files}, json.loads(man.read_text(encoding="utf-8")) if man.exists() else None)

    def _interp(self, key, state, a, b):
        ck = (key, state, a, b)
        if ck not in self._cache:
            v = self.a[key][:, :, :, state - 1, a, b]
            vp = np.concatenate([v[:, :, -self._pad :], v, v[:, :, : self._pad]], axis=2)
            axes = (self.f, self.theta, self._phi_pad)
            self._cache[ck] = (
                RegularGridInterpolator(axes, vp.real, method="cubic"),
                RegularGridInterpolator(axes, vp.imag, method="cubic") if np.iscomplexobj(vp) else None,
            )
        return self._cache[ck]

    def _check_orders(self, f, th, ph):
        margin = higher_order_margin(f, th, ph)
        if np.any(margin < 1.0):
            i = int(np.argmin(margin))
            raise OutOfDomain(
                "NOT_AVAILABLE_FROM_CURRENT_DATA: a higher diffraction order propagates at "
                f"(f, theta, phi) = ({f[i]:g} GHz, {th[i]:g}, {ph[i]:g}); the stored (0,0) reflection is not the whole response"
            )

    def guard(self, f_ghz, theta_deg, phi_deg):
        """True where the query lies in the Wood-anomaly guard band just below diffraction onset (the response varies faster than the grid resolves)."""
        return higher_order_margin(f_ghz, theta_deg, phi_deg) < WOOD_GUARD

    def _check(self, f_ghz, theta_deg):
        f_ghz, theta_deg = np.atleast_1d(f_ghz).astype(float), np.atleast_1d(theta_deg).astype(float)
        bad = []
        if f_ghz.min() < self.f[0] - 1e-9 or f_ghz.max() > self.f[-1] + 1e-9:
            bad.append(f"frequency outside {self.f[0]:g}-{self.f[-1]:g} GHz")
        if theta_deg.min() < self.theta[0] - 1e-9 or theta_deg.max() > self.theta[-1] + 1e-9:
            bad.append(f"theta outside {self.theta[0]:g}-{self.theta[-1]:g} deg")
        if bad:
            raise OutOfDomain("NOT_AVAILABLE_FROM_CURRENT_DATA: " + "; ".join(bad))

    def _query(self, f_ghz, theta_deg, phi_deg):
        """Broadcast the query, check the domain, clip the 1e-9 boundary slack onto the grid, wrap phi: returns raveled (f, theta, phi) arrays."""
        f, th, ph = np.broadcast_arrays(
            *(np.atleast_1d(np.asarray(x, dtype=float)) for x in (f_ghz, theta_deg, phi_deg))
        )
        self._check(f, th)
        self._check_orders(f.ravel(), th.ravel(), ph.ravel())
        f = np.clip(f.ravel(), self.f[0], self.f[-1])
        th = np.clip(th.ravel(), self.theta[0], self.theta[-1])
        return f, th, ph.ravel()

    def structure_circular(self, f_ghz, theta_deg, phi_deg, state, key="c_corr"):
        """Interpolated structure-frame circular reflection, shape (n, 2, 2) for the n points of the broadcast query (raveled in C order)."""
        f, th, ph = self._query(f_ghz, theta_deg, phi_deg)
        q = np.stack([f, th, np.mod(ph, 360.0)], axis=-1)
        out = np.zeros((q.shape[0], 2, 2), dtype=complex)
        for a in (0, 1):
            for b in (0, 1):
                re, im = self._interp(key, state, a, b)
                out[:, a, b] = re(q) + (1j * im(q) if im is not None else 0.0)
        return out

    def te_tm(self, f_ghz, theta_deg, phi_deg, state, key="c_corr"):
        """The same response as the (TE, TM) power-wave reflection in the incidence-plane frame, shape (n, 2, 2) in the same point order."""
        f, th, ph = self._query(f_ghz, theta_deg, phi_deg)
        c = self.structure_circular(f, th, ph, state, key)
        return np.stack([te_tm_from_circular(from_structure_frame(c[i], np.radians(ph[i]))) for i in range(c.shape[0])])

    def uncertainty(self, f_ghz, theta_deg, phi_deg, state):
        """Interpolated absolute uncertainty per circular element (distance between the p = 1 and p = 2 mesh extrapolations + azimuth spread of the correction; NOT the model uncertainty)."""
        z = self.structure_circular(f_ghz, theta_deg, phi_deg, state, key="uncertainty")
        return np.abs(z.real)

    def label(self, theta_deg):
        return "EXTRAPOLATED_BEYOND_R3_TARGETS" if float(theta_deg) > THETA_LABEL_LIMIT_DEG else "within_R3_targets"
