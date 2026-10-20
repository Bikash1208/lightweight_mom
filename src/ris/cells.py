"""Parametrised 2-bit PIN-diode RIS cells (Phase 20 R3 item 2): the CP-2 octagonal PB cell, the LP-1 three-section patch and the LP-2
stub-loaded patch. Each cell is a frozen dataclass of dimensions (STARTING geometries -- item 3 optimises them, none is a design) with
`unit_cell(state, kit)` returning the `UnitCell` (composite/polygon geometry + diode loads) the MoM solves.

Common modelling choices (all labelled approximations, `VIA_OMITTED_APPROXIMATION` among them)
* RT/duroid 5880 1.575 mm, PEC ground, zero-thickness metal; the DC bias via at the patch centre / bias layer behind the ground are omitted.
* A diode is the ideal series model of `src.ris.diodes` on a small conductor region (`diode_span` along the conduction line, the line
  width across it), `transverse_conducts=False`; the region is metal in the pattern (R1 finding). The HFSS chain representation of R2 is
  for the check packages (item 5), not for the MoM here.
* States are tuples of "ON"/"OFF" per diode site: LP cells use (d1, d2) in {OFF,ON}^2 (4 states); the CP cell uses a one-hot code (one of
  four diodes ON, three OFF) = element rotation 0/45/90/135 deg.

`DiodeKit` pairs an ON and an OFF `DiodeModel` (SMP1345-079LF, or MADP-000907 with an assumed L_s) so a cell can be evaluated with either.
"""

import math
from dataclasses import dataclass

from ..composite_geometry import AxisRect, CompositeGeometry
from ..geometry import UnitCell
from ..materials import Metal, Substrate
from ..periodic_bc import RectangularLattice
from ..polygon_geometry import PolyComposite, oriented_rect, regular_polygon
from .diodes import diode_load, madp_000907_off, madp_000907_on, smp1345_off, smp1345_on

RT5880 = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
ON, OFF = "ON", "OFF"


@dataclass(frozen=True)
class DiodeKit:
    name: str
    on: object
    off: object

    def model(self, state):
        return self.on if state == ON else self.off


def smp1345_kit(bias="10mA", r_p=None):
    return DiodeKit("SMP1345-079LF", smp1345_on(bias), smp1345_off(r_p))


def madp_kit(l_s=None):
    """MADP-000907-14020; `l_s=None` is the NOT_AVAILABLE nominal (L_s = 0), a number an ASSUMED sensitivity value."""
    return DiodeKit("MADP-000907-14020", madp_000907_on(l_s), madp_000907_off(l_s))


def ideal_kit(r_on=1e-3, c_off=None):
    """Ideal switch for symmetry tests: ON = (almost) short, OFF = open (`open_circuit`-like, a 1e-6 pF capacitor)."""
    from .diodes import DiodeModel

    return DiodeKit("ideal", DiodeModel("ideal ON", R_s=r_on), DiodeModel("ideal OFF", C_T=c_off or 1e-18))


def _check_states(states, n):
    if len(states) != n or any(s not in (ON, OFF) for s in states):
        raise ValueError(f"state must be {n} entries of 'ON'/'OFF', got {states!r}")


def _cell(period, geometry, loads, backing="pec"):
    return UnitCell(RectangularLattice(period, period), geometry, RT5880, Metal(), backing=backing, loads=tuple(loads))


# ======================================================================================================================
@dataclass(frozen=True)
class Lp1Cell:
    """LP-1: three-section patch (left | middle | right along x, the polarisation across the gaps), a diode across each of the two gaps.
    Asymmetric (section lengths, strip widths and positions all differ) so the four (d1, d2) combinations are four distinct states. The
    DEFAULTS are the best of a coarse 70-point random scan (min pairwise |dGamma_TM| 0.15 at 10 GHz, SMP1345, phases only ~35 deg apart):
    a valid starting point, nowhere near a 2-bit cell. Item 3 optimises."""

    period: float = 12.0e-3
    width: float = 6.0e-3  # patch extent in y
    left: float = 3.5e-3  # section lengths along x
    middle: float = 2.5e-3
    right: float = 1.5e-3
    diode_span: float = 1.0e-3  # conductor region carrying the diode, along x
    pad: float = 0.25e-3  # strip conductor each side of the diode span (inside the gap)
    strip_width1: float = 0.5e-3  # width of the diode bridge strip across gap 1 (left | middle)
    strip_width2: float = 2.0e-3  # ... across gap 2 (middle | right)
    strip_y1: float = 0.75e-3  # y of the bridge strip across gap 1
    strip_y2: float = 0.0  # y of the bridge strip across gap 2

    kind = "lp1"
    n_diodes = 2

    @property
    def gap(self):
        return self.diode_span + 2 * self.pad

    def states(self):
        return [(a, b) for a in (OFF, ON) for b in (OFF, ON)]

    def _x_edges(self):
        total = self.left + self.gap + self.middle + self.gap + self.right
        x0 = -total / 2
        l_end = x0 + self.left
        m_start = l_end + self.gap
        m_end = m_start + self.middle
        r_start = m_end + self.gap
        return x0, l_end, m_start, m_end, r_start, r_start + self.right

    def geometry(self):
        x0, l_end, m_start, m_end, r_start, x1 = self._x_edges()
        rect = lambda a, b, cy, w: AxisRect(0.5 * (a + b), cy, b - a, w)  # noqa: E731
        metal = (
            rect(x0, l_end, 0.0, self.width),
            rect(m_start, m_end, 0.0, self.width),
            rect(r_start, x1, 0.0, self.width),
            rect(l_end, m_start, self.strip_y1, self.strip_width1),
            rect(m_end, r_start, self.strip_y2, self.strip_width2),
        )
        return CompositeGeometry(metal)

    def diode_sites(self):
        """[(x, y, length, width, direction)] of the diode regions."""
        _, l_end, m_start, m_end, r_start, _ = self._x_edges()
        return [
            (0.5 * (l_end + m_start), self.strip_y1, self.diode_span, self.strip_width1, "x"),
            (0.5 * (m_end + r_start), self.strip_y2, self.diode_span, self.strip_width2, "x"),
        ]

    def unit_cell(self, state, kit):
        _check_states(state, 2)
        loads = [diode_load(kit.model(s), *site) for s, site in zip(state, self.diode_sites(), strict=True)]
        return _cell(self.period, self.geometry(), loads)


# ======================================================================================================================
@dataclass(frozen=True)
class Lp2Cell:
    """LP-2 (proposal of this project, no published precedent): a rectangular patch with two diode-switched stubs of different length.
    Stub 1 continues the patch along +x, stub 2 along -x (polarisation along x); a diode in series with each stub (OFF: the stub is cut
    by the diode capacitance, ON: it extends the resonant length). Single layer, in-plane diodes, no via other than the omitted bias via.
    Four distinct states because stub 1 != stub 2. DEFAULTS: best of a coarse 70-point random scan (min pairwise |dGamma_TM| 0.20 at 10 GHz,
    SMP1345, phases 101/90/72/45 deg: a ~55 deg spread, a valid starting point, nowhere near a 2-bit cell; item 3 optimises)."""

    period: float = 12.0e-3
    patch_length: float = 6.0e-3
    patch_width: float = 3.5e-3
    stub_width: float = 1.0e-3
    diode_span: float = 1.0e-3
    pad: float = 0.25e-3
    stub1_length: float = 0.75e-3  # conductor beyond the diode
    stub2_length: float = 0.25e-3
    stub1_y: float = 0.5e-3
    stub2_y: float = 1.0e-3

    kind = "lp2"
    n_diodes = 2

    def states(self):
        return [(a, b) for a in (OFF, ON) for b in (OFF, ON)]

    def _stub(self, sign, length, y):
        x_edge = sign * self.patch_length / 2
        x_pad0 = x_edge  # pad starts at the patch edge
        x_d0 = x_edge + sign * self.pad
        x_d1 = x_d0 + sign * self.diode_span
        x_end = x_d1 + sign * (self.pad + length)
        rect = AxisRect(0.5 * (x_pad0 + x_end), y, abs(x_end - x_pad0), self.stub_width)
        site = (0.5 * (x_d0 + x_d1), y, self.diode_span, self.stub_width, "x")
        return rect, site

    def geometry(self):
        s1, _ = self._stub(+1, self.stub1_length, self.stub1_y)
        s2, _ = self._stub(-1, self.stub2_length, self.stub2_y)
        return CompositeGeometry((AxisRect(0.0, 0.0, self.patch_length, self.patch_width), s1, s2))

    def diode_sites(self):
        return [self._stub(+1, self.stub1_length, self.stub1_y)[1], self._stub(-1, self.stub2_length, self.stub2_y)[1]]

    def unit_cell(self, state, kit):
        _check_states(state, 2)
        loads = [diode_load(kit.model(s), *site) for s, site in zip(state, self.diode_sites(), strict=True)]
        return _cell(self.period, self.geometry(), loads)


# ======================================================================================================================
@dataclass(frozen=True)
class Lp3Cell:
    """LP-3 (proposal of this project; the brief's 3-diode fallback if no 2-diode LP cell reaches 90 +- 15 deg steps): a rectangular patch with
    one diode-switched stub on the -x side (diode d1) and a CHAIN of two diode-switched segments on the +x side (d2 connects segment 1, d3
    connects segment 2 behind it). 8 diode codes (d1, d2, d3); the physically distinct ones are 6 ((d2, d3) = (0,1) is (0,0) with a floating
    segment), from which the optimiser picks the 4 closest to the 0/90/180/270 lattice. Single layer, in-plane diodes."""

    period: float = 14.0e-3
    patch_length: float = 5.0e-3
    patch_width: float = 3.5e-3
    stub_width: float = 1.0e-3
    diode_span: float = 1.0e-3
    pad: float = 0.25e-3
    left_length: float = 0.5e-3  # conductor beyond d1
    right1_length: float = 0.5e-3  # conductor beyond d2
    right2_length: float = 0.5e-3  # conductor beyond d3
    left_y: float = 0.0
    right_y: float = 0.0

    kind = "lp3"
    n_diodes = 3

    def states(self):
        return [(a, b, c) for a in (OFF, ON) for b in (OFF, ON) for c in (OFF, ON)]

    def _chain(self, sign, y, lengths):
        """Rectangles and diode sites of a chain of (pad | diode | pad | metal) segments starting at the patch edge, direction `sign`."""
        rects, sites, cur = [], [], sign * self.patch_length / 2
        for length in lengths:
            d0 = cur + sign * self.pad
            d1 = d0 + sign * self.diode_span
            end = d1 + sign * (self.pad + length)
            rects.append(AxisRect(0.5 * (cur + end), y, abs(end - cur), self.stub_width))
            sites.append((0.5 * (d0 + d1), y, self.diode_span, self.stub_width, "x"))
            cur = end
        return rects, sites

    def geometry(self):
        left, _ = self._chain(-1, self.left_y, [self.left_length])
        right, _ = self._chain(+1, self.right_y, [self.right1_length, self.right2_length])
        return CompositeGeometry((AxisRect(0.0, 0.0, self.patch_length, self.patch_width), *left, *right))

    def diode_sites(self):
        _, ls = self._chain(-1, self.left_y, [self.left_length])
        _, rs = self._chain(+1, self.right_y, [self.right1_length, self.right2_length])
        return [*ls, *rs]

    def unit_cell(self, state, kit):
        _check_states(state, 3)
        loads = [diode_load(kit.model(s), *site) for s, site in zip(state, self.diode_sites(), strict=True)]
        return _cell(self.period, self.geometry(), loads)


# ======================================================================================================================
@dataclass(frozen=True)
class Cp2Cell:
    """CP-2: regular octagonal patch (flat edges facing 0/45/90/135/... deg) with four radial delay-line stubs at 45 deg spacing,
    one diode each; one-hot state k = the diode of stub k ON, the other three OFF. Pancharatnam-Berry: the ON stub sets the anisotropy axis, so
    the four states are the same cell rotated by 45 deg steps (CP phase = 2 x rotation = 0/90/180/270 deg). The centre via is omitted
    (`VIA_OMITTED_APPROXIMATION`). The radial layout of a stub, from the octagon's flat edge (apothem a) outward:
    `overlap` into the patch | pad | diode span | line (`line`); everything but the diode span is conductor."""

    period: float = 13.5e-3
    circumradius: float = 4.5e-3
    stub_angles_deg: tuple = (0.0, 45.0, 90.0, 135.0)
    stub_width: float = 0.5e-3
    overlap: float = 0.1e-3
    pad: float = 0.25e-3
    diode_span: float = 0.75e-3
    line: float = 0.5e-3
    transverse_conducts: bool = (
        True  # False blocks transverse current in the diode region: on a 45-degree strip the Cartesian rooftops
    )
    # then need Jx = Jy exactly and the block kills the longitudinal current too (see tests); True = a full-width sheet in the diode region

    kind = "cp2"
    n_diodes = 4

    @property
    def apothem(self):
        return self.circumradius * math.cos(math.pi / 8)

    def states(self):
        return [
            tuple(ON if i == k else OFF for i in range(len(self.stub_angles_deg)))
            for k in range(len(self.stub_angles_deg))
        ]

    def stub_extent(self):
        """Radial distance from the patch centre to the tip of a stub."""
        return self.apothem + self.pad + self.diode_span + self.pad + self.line

    def _stub(self, ang_deg):
        a = math.radians(ang_deg)
        r0 = self.apothem - self.overlap
        r_d0 = self.apothem + self.pad
        r_d1 = r_d0 + self.diode_span
        r_end = r_d1 + self.pad + self.line
        c = lambda r: (r * math.cos(a), r * math.sin(a))  # noqa: E731
        mid = c(0.5 * (r0 + r_end))
        poly = oriented_rect(mid[0], mid[1], r_end - r0, self.stub_width, a)
        dm = c(0.5 * (r_d0 + r_d1))
        return poly, (dm[0], dm[1], self.diode_span, self.stub_width, ang_deg)

    def geometry(self):
        octagon = regular_polygon(0.0, 0.0, self.circumradius, 8, rotation=math.pi / 8)
        return PolyComposite((octagon, *[self._stub(a)[0] for a in self.stub_angles_deg]))

    def diode_sites(self):
        return [self._stub(a)[1] for a in self.stub_angles_deg]

    def unit_cell(self, state, kit):
        _check_states(state, len(self.stub_angles_deg))
        if self.stub_extent() > self.period / 2 - 0.25e-3 and any(
            abs(math.cos(math.radians(a))) > 0.99 or abs(math.sin(math.radians(a))) > 0.99 for a in self.stub_angles_deg
        ):
            raise ValueError("an axis-aligned stub reaches the cell boundary: shorten `line` or enlarge `period`")
        loads = [
            diode_load(kit.model(s), *site, transverse_conducts=self.transverse_conducts)
            for s, site in zip(state, self.diode_sites(), strict=True)
        ]
        return _cell(self.period, self.geometry(), loads)
