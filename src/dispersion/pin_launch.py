"""Phase 26 (launch fraction): how a monopole / coax pin on a grounded slab splits its radiated power between
the TM0 surface wave and direct space-wave radiation.

Model. The pin is a vertical electric current filament `I(z)` from the ground plane (z = -h) up to its top at
`z_top` (z = 0 is the substrate surface), optionally ending in a top-hat disk of radius `a` that carries a
radial surface current `I_rho(rho) = I(z_top) (1 - rho^2/a^2)` (uniform charge). The layered medium is the
grounded slab of this project (air | eps_c, h | PEC), unmodified: the same TM transmission line along z that
`green_function.LayeredGreenFunction` is built on, written here for a VERTICAL current.

Plane-wave (kt) form, exp(+j omega t), fields exp(-j kt x), TM_z fields E_x = V, H_y = I:

    dV/dz = -j kz^2/(omega eps) I + kt/(omega eps) J_z      dI/dz = -j omega eps V - J_x
    E_z   = -kt/(omega eps) I + j J_z/(omega eps)

A vertical segment `I_n dz` is a voltage jump `kt I_n dz / (omega eps_n)`; the disk is a current jump
`-j F(kt)`, `F = int_0^a I_rho(rho) J1(kt rho) drho`. With the left (ground-side) and right (outgoing)
homogeneous solutions `(V_L, I_L)`, `(V_R, I_R)` and `D = V_R I_L - V_L I_R = Z_0 cos(k_1 h) + j Z_1 sin(k_1 h)`
(`D = 0` is the TM0 dispersion equation `k_1 tan(k_1 h) = eps_r alpha_0`), the complex power integrand is

    S(kt) = sum_n I_n dz E_z,n + conj(J_x) V_d,        P = -1/2 (1/2 pi) int Re[S] kt dkt.

* `P_space`: the integral over `0 <= kt <= k0` (`kt = k0 sin theta`).
* `P_sw`: the TM0 pole at `kt = beta > k0`, where `S` is purely imaginary except for the pole itself:
  `P_sw = beta |Im Res S| / 4` (residue from a contour around `beta`).
Their sum is the radiated power; `launch_fraction = P_sw / (P_sw + P_space)`. An independent lossy-slab
numerical integral of `Re S` over the whole real axis (`power_balance`) checks the split.

Approximations (documented, not hidden): thin filament (radius enters only the characteristic impedance of
the current-distribution model and the top-hat capacitance), no coax aperture / gap, no feed-line detail,
disk radial current profile assumed, the sinusoidal current distribution carries a thin-wire transmission-line
phase `k(z)` (k0 in air, k0 sqrt(eps_r) inside the slab). The input impedance is therefore an ESTIMATE
(`R` is the radiated-power resistance, `X` the thin-wire transmission-line reactance).
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import j1

from ..constants import EPS0, ETA0
from ..frequency import angular_frequency, free_space_wavenumber


def _root_kz(k2):
    """`sqrt(k2)` on the physical branch Im <= 0 (Re >= 0 when real)."""
    r = np.sqrt(k2 + 0j)
    flip = (r.imag > 1e-300) | ((abs(r.imag) <= 1e-300) & (r.real < 0.0))
    return -r if flip else r


def _sin_over_k(kz, z):
    return z if abs(kz * z) < 1e-9 else np.sin(kz * z) / kz


@dataclass(frozen=True)
class Launcher:
    """Pin geometry. `top_above_surface_m`: height of the pin top above the substrate surface (negative:
    fully inside the substrate, but above the ground). `disk_radius_m`: top-hat disk, 0 = none."""

    top_above_surface_m: float
    pin_radius_m: float = 0.305e-3
    disk_radius_m: float = 0.0


class _Layer:
    """The grounded slab at one `kt`: left / right solutions and `D`."""

    def __init__(self, kt, k0, eps_c, h, y_sheet=0.0):
        self.kt, self.h = kt, h
        self.y_sheet = y_sheet
        self.omega = k0 * 299792458.0
        self.k0 = k0
        self.eps_c = eps_c
        self.kz0 = _root_kz(k0 * k0 - kt * kt)
        self.kz1 = _root_kz(eps_c * k0 * k0 - kt * kt)
        w_e0 = self.omega * EPS0
        self.z0 = self.kz0 / w_e0
        self.z1 = self.kz1 / (w_e0 * eps_c)
        self.w_e0 = w_e0
        self.il0 = np.cos(self.kz1 * h)
        self.vl0 = -1j * self.z1 * np.sin(self.kz1 * h)
        if y_sheet == 0.0:
            self.d = self.z0 * self.il0 + 1j * self.z1 * np.sin(self.kz1 * h)
            self.il0_air = self.il0
        else:  # a shunt sheet admittance at z = 0: the H_y jump is -Y V
            self.il0_air = self.il0 - y_sheet * self.vl0
            self.d = self.z0 * self.il0_air - self.vl0

    @staticmethod
    def _propagate(v0, i0, kz, w_eps, z):
        """(V, I) a distance `z` from a plane carrying (v0, i0), in a medium of `kz`, `omega eps`."""
        c = np.cos(kz * z)
        return v0 * c - 1j * i0 * kz * np.sin(kz * z) / w_eps, i0 * c - 1j * v0 * w_eps * _sin_over_k(kz, z)

    def left(self, z):
        """(I, V) of the solution that vanishes (E_x = 0) on the ground plane."""
        if z <= 0.0:
            s = z + self.h
            return np.cos(self.kz1 * s), -1j * self.z1 * np.sin(self.kz1 * s)
        v, i = self._propagate(self.vl0, self.il0_air, self.kz0, self.w_e0, z)
        return i, v

    def right(self, z):
        """(I, V) of the outgoing solution (decaying / radiating upward)."""
        if z >= 0.0:
            e = np.exp(-1j * self.kz0 * z)
            return e, self.z0 * e
        v, i = self._propagate(self.z0, 1.0 + self.y_sheet * self.z0, self.kz1, self.w_e0 * self.eps_c, z)
        return i, v


def _current_profile(z, launcher, f, eps_r, h, ground_z):
    """Sinusoidal thin-wire current of the pin, normalised to 1 at the feed (ground), with the top-hat
    capacitive end condition `tan(Phi) = omega C Z_w`, `C = 8 eps0 a`."""
    k0 = free_space_wavenumber(f)
    z_top = launcher.top_above_surface_m
    length = z_top - ground_z

    def phase(zz):  # integral of k(z') from zz up to the top
        sl = np.maximum(min(0.0, z_top) - np.asarray(zz, dtype=float), 0.0)  # length inside the slab
        ab = np.maximum(z_top - np.maximum(np.asarray(zz, dtype=float), 0.0), 0.0)  # length in air
        return k0 * np.sqrt(eps_r) * sl + k0 * ab

    z_w = 60.0 * max(np.log(2.0 * max(length, 4 * launcher.pin_radius_m) / launcher.pin_radius_m) - 1.0, 1.0)
    cap = 8.0 * EPS0 * launcher.disk_radius_m
    phi = np.arctan(angular_frequency(f) * cap * z_w) if cap > 0 else 0.0
    prof = np.sin(phase(z) + phi)
    norm = np.sin(phase(ground_z) + phi)
    return prof / norm, phi, z_w


def _disk_f(kt, a, i_top):
    if a <= 0.0 or i_top == 0.0:
        return 0.0
    val, _ = quad(lambda r: (1.0 - r * r / (a * a)) * j1(kt * r), 0.0, a, limit=50)
    return i_top * val


class PinLaunch:
    """Reusable evaluator for one `(frequency, substrate, launcher)`."""

    def __init__(self, frequency_hz, eps_r, thickness_m, launcher: Launcher, tan_delta=0.0, n_seg=24, y_sheet_s=0.0):
        self.f, self.eps_r, self.h = frequency_hz, eps_r, thickness_m
        self.y_sheet = y_sheet_s  # patch layer around the feed, shunt sheet admittance [S] (0 = bare slab)
        self.k0 = free_space_wavenumber(frequency_hz)
        self.omega = angular_frequency(frequency_hz)
        self.eps_c = eps_r * (1.0 - 1j * tan_delta)
        self.launcher = launcher
        if launcher.top_above_surface_m <= -thickness_m:
            raise ValueError("the pin top must be above the ground plane")
        top = launcher.top_above_surface_m
        edges = np.linspace(-thickness_m, top, n_seg + 1)
        self.dz = np.diff(edges)
        self.zs = 0.5 * (edges[1:] + edges[:-1])
        # split a segment that straddles the surface? edges are uniform; the slab interface lies in one
        # segment at worst -- refine by inserting z = 0 as an edge when it falls inside
        if -1e-12 < top and not np.any(np.isclose(edges, 0.0)):
            e = np.sort(np.append(edges, 0.0))
            self.dz = np.diff(e)
            self.zs = 0.5 * (e[1:] + e[:-1])
        self.i_seg, self.phi, self.z_w = _current_profile(self.zs, launcher, frequency_hz, eps_r, thickness_m, -thickness_m)
        self.eps_seg = np.where(self.zs < 0.0, self.eps_c, 1.0) * EPS0
        self.i_top = float(_current_profile(top, launcher, frequency_hz, eps_r, thickness_m, -thickness_m)[0])

    # ------------------------------------------------------------------------------------------
    def s_of_kt(self, kt):
        lay = _Layer(kt, self.k0, self.eps_c, self.h, self.y_sheet)
        zs = list(self.zs)
        a_disk = self.launcher.disk_radius_m
        dv = list(kt * self.i_seg * self.dz / (self.omega * self.eps_seg))
        di = [0.0] * len(zs)
        if a_disk > 0.0:
            zs.append(self.launcher.top_above_surface_m)
            dv.append(0.0)
            di.append(-1j * _disk_f(kt.real if np.iscomplexobj(kt) else kt, a_disk, self.i_top))
        zs = np.array(zs)
        dv, di = np.array(dv, dtype=complex), np.array(di, dtype=complex)
        ls = [lay.left(z) for z in zs]
        rs = [lay.right(z) for z in zs]
        n = len(zs)
        a_r = np.array([(dv[s] * ls[s][0] - di[s] * ls[s][1]) / lay.d for s in range(n)])
        a_l = np.array([(dv[s] * rs[s][0] - di[s] * rs[s][1]) / lay.d for s in range(n)])
        i_tl = np.zeros(n, dtype=complex)
        v_tl = np.zeros(n, dtype=complex)
        for o in range(n):
            for s in range(n):
                if zs[o] > zs[s]:
                    i_tl[o] += a_r[s] * rs[o][0]
                    v_tl[o] += a_r[s] * rs[o][1]
                elif zs[o] < zs[s]:
                    i_tl[o] += a_l[s] * ls[o][0]
                    v_tl[o] += a_l[s] * ls[o][1]
                else:
                    i_tl[o] += 0.5 * (a_r[s] * rs[o][0] + a_l[s] * ls[o][0])
                    v_tl[o] += 0.5 * (a_r[s] * rs[o][1] + a_l[s] * ls[o][1])
        nseg = len(self.zs)
        s_val = np.sum(self.i_seg * self.dz * (-kt / (self.omega * self.eps_seg)) * i_tl[:nseg])
        if a_disk > 0.0:
            s_val += np.conj(1j * _disk_f(kt.real if np.iscomplexobj(kt) else kt, a_disk, self.i_top)) * v_tl[nseg]
        return complex(s_val)

    # ------------------------------------------------------------------------------------------
    def beta(self):
        """TM0 constant of the bare slab (lossless)."""
        eps = self.eps_r

        def f(kt):
            lay = _Layer(kt, self.k0, eps, self.h, self.y_sheet)
            return (lay.d / 1j).real

        return float(brentq(f, self.k0 * (1 + 1e-9), self.k0 * np.sqrt(eps) * (1 - 1e-9)))

    def space_power(self, n_theta=400):
        th = (np.arange(n_theta) + 0.5) * (np.pi / 2) / n_theta
        kt = self.k0 * np.sin(th)
        vals = np.array([self.s_of_kt(k).real * k * self.k0 * np.cos(t) for k, t in zip(kt, th, strict=True)])
        return float(-0.5 / (2 * np.pi) * np.sum(vals) * (np.pi / 2) / n_theta)

    def sw_power(self, radius_rel=2e-3, n_circle=24):
        beta = self.beta()
        r = min(radius_rel * self.k0, 0.3 * (beta - self.k0))  # never enclose the branch point at k0
        acc = 0.0 + 0.0j
        for m in range(n_circle):
            ph = 2 * np.pi * (m + 0.5) / n_circle
            kt = beta + r * np.exp(1j * ph)
            acc += self.s_of_kt(kt) * r * np.exp(1j * ph)  # residue = (1/2 pi j) oint S dkt -> mean of S r e^{j ph}
        res = acc / n_circle
        return float(beta * abs(res.imag) / 4.0), beta

    def split(self):
        p_space = self.space_power()
        p_sw, beta = self.sw_power()
        tot = p_space + p_sw
        out = {
            "p_space_w_per_a2": p_space,
            "p_sw_w_per_a2": p_sw,
            "launch_fraction": p_sw / tot if tot > 0 else float("nan"),
            "beta_over_k0": beta / self.k0,
            "r_rad_ohm": 2.0 * tot,
            "z_wire_ohm": self.z_w,
        }
        out["x_est_ohm"] = self._reactance_estimate()
        return out

    def _reactance_estimate(self):
        """Thin-wire transmission-line reactance seen at the feed: `-Z_w cot(Phi_total)` where `Phi_total`
        is the end-condition phase plus the electrical length of the pin."""
        k0 = self.k0
        top = self.launcher.top_above_surface_m
        inside = min(top, 0.0) + self.h
        above = max(top, 0.0)
        theta = k0 * np.sqrt(self.eps_r) * inside + k0 * above + self.phi
        return float(-self.z_w / np.tan(theta))

    def power_balance(self, tan_delta=1e-3, kt_max_over_k0=None, n_per_k0=4000):
        """Independent numerical total: `-1/2 (1/2 pi) int Re[S] kt dkt` over the whole real axis with a
        lossy slab (the pole is then a finite-width resonance). Compare with `p_space + p_sw`."""
        lossy = PinLaunch(
            self.f, self.eps_r, self.h, self.launcher, tan_delta=tan_delta, n_seg=len(self.zs), y_sheet_s=self.y_sheet
        )
        k_hi = (kt_max_over_k0 or 1.6 * np.sqrt(self.eps_r)) * self.k0
        beta = self.beta()
        pts = np.unique(np.concatenate([np.linspace(0, k_hi, n_per_k0), beta + np.linspace(-0.02, 0.02, 4001) * self.k0]))
        vals = np.array([lossy.s_of_kt(k).real * k for k in pts])
        return float(-0.5 / (2 * np.pi) * np.trapezoid(vals, pts))
