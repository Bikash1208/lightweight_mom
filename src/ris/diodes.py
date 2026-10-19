"""PIN-diode equivalent circuits as frequency-dependent series impedances Z_L(f) (Phase 20 R2 item 1).

ON  : R_s + j w L_s.
OFF : j w L_s + Z_par,  Z_par = C_T in parallel with R_p (R_p = None: lossless C_T; the datasheet gives no R_p).
Generic: any of R_s, L_s, C_T (series C), R_p (parallel to C_T).

Every value carries its provenance. Datasheet values are `PROVISIONAL_DATASHEET` (no measured de-embedding): they are
specified at 1 MHz / 100 MHz and the datasheet covers 10 MHz - 6 GHz, so use at 8 - 12 GHz is an EXTRAPOLATION. Values the
datasheet does not give are `NOT_AVAILABLE_FROM_CURRENT_DATA` and are never invented; sweeps over them are labelled
`ASSUMED_SENSITIVITY_NOT_DEVICE_DATA`.
"""

from dataclasses import dataclass, field, replace

import numpy as np

PROVISIONAL = "PROVISIONAL_DATASHEET"
NOT_AVAILABLE = "NOT_AVAILABLE_FROM_CURRENT_DATA"
ASSUMED = "ASSUMED_SENSITIVITY_NOT_DEVICE_DATA"

DOC = "Skyworks SMP1345 Series datasheet 200046U, 2018-06-19"
RANGE_NOTE = "datasheet covers 10 MHz-6 GHz (Rs at 100 MHz, CT at 1 MHz); use at 8-12 GHz is an extrapolation, no measured de-embedding"


@dataclass(frozen=True)
class Provenance:
    value: "float | None"
    unit: str
    status: str
    source: str
    condition: str = ""
    maximum: "float | None" = None  # datasheet maximum, where one is given (tolerance bounds)


@dataclass(frozen=True)
class DiodeModel:
    name: str
    R_s: float = 0.0  # series resistance [ohm]
    L_s: float = 0.0  # series (package) inductance [H]
    C_T: "float | None" = None  # junction capacitance [F]; None: no capacitor branch
    R_p: "float | None" = None  # resistance in parallel with C_T [ohm]; None: infinite (lossless C_T)
    provenance: "dict[str, Provenance]" = field(default_factory=dict)

    def __post_init__(self):
        for k in ("R_s", "L_s"):
            v = getattr(self, k)
            if not np.isfinite(v) or v < 0:
                raise ValueError(f"{k} must be finite and >= 0, got {v}")
        for k in ("C_T", "R_p"):
            v = getattr(self, k)
            if v is not None and (not np.isfinite(v) or v <= 0):
                raise ValueError(f"{k} must be finite and > 0 (or None), got {v}")
        if self.R_p is not None and self.C_T is None:
            raise ValueError("R_p is in parallel with C_T: C_T is required")

    def impedance(self, frequency):
        """Z_L(f) [ohm] for a frequency (Hz) scalar or array. Passive: Re Z >= 0."""
        f = np.asarray(frequency, float)
        if np.any(f <= 0):
            raise ValueError("frequency must be > 0")
        w = 2 * np.pi * f
        z = self.R_s + 1j * w * self.L_s
        if self.C_T is not None:
            y = 1j * w * self.C_T + (0.0 if self.R_p is None else 1.0 / self.R_p)
            z = z + 1.0 / y
        return complex(z) if np.ndim(z) == 0 else z

    @property
    def series_resonance_hz(self):
        """f0 = 1 / (2 pi sqrt(L_s C_T)) (where Im Z = 0 for R_p = None); None without both L_s and C_T."""
        if self.C_T is None or self.L_s <= 0:
            return None
        return 1.0 / (2 * np.pi * np.sqrt(self.L_s * self.C_T))

    def scaled(self, **factors):
        """Copy with the named parameters multiplied (sensitivity runs), e.g. `scaled(C_T=1.2)`. Provenance is marked as perturbed."""
        kw, prov = {}, dict(self.provenance)
        for k, fac in factors.items():
            if k not in ("R_s", "L_s", "C_T", "R_p"):
                raise ValueError(f"unknown parameter {k!r}")
            v = getattr(self, k)
            if v is None:
                raise ValueError(f"{k} is not set on this model")
            kw[k] = v * fac
            if k in prov:
                p = prov[k]
                prov[k] = replace(
                    p, value=kw[k], status="PERTURBED_SENSITIVITY", condition=f"{p.condition}; x{fac:g} of nominal"
                )
        return replace(
            self, name=f"{self.name} [{', '.join(f'{k} x{v:g}' for k, v in factors.items())}]", provenance=prov, **kw
        )


_LS = Provenance(0.7e-9, "H", PROVISIONAL, f"{DOC}, Table 1 (p. 1), SC-79 package", "package series inductance L_S")


def smp1345_on(bias="10mA"):
    """SMP1345-079LF forward-biased: R_s = 1.5 ohm (10 mA, typ) or 3.5 ohm (1 mA, typ), both at 100 MHz, + L_s = 0.7 nH."""
    table = {"10mA": 1.5, "1mA": 3.5}
    if bias not in table:
        raise ValueError(f"bias must be one of {sorted(table)}")
    rs = Provenance(
        table[bias],
        "ohm",
        PROVISIONAL,
        f"{DOC}, Table 3 (p. 2)",
        f"I_F = {bias[:-2]} mA, F = 100 MHz, typ; {RANGE_NOTE}",
    )
    return DiodeModel(f"SMP1345-079LF ON ({bias})", R_s=rs.value, L_s=_LS.value, provenance={"R_s": rs, "L_s": _LS})


def smp1345_off(r_p=None):
    """SMP1345-079LF reverse-biased: C_T = 0.19 pF (1 V, 1 MHz, typ) in series with L_s = 0.7 nH, lossless.

    The datasheet gives no R_p (`NOT_AVAILABLE_FROM_CURRENT_DATA`); `r_p` is an assumed sensitivity value, not device data. C_T at
    0 V (a RIS board's OFF bias) is not given and may be slightly higher.
    """
    ct = Provenance(
        0.19e-12,
        "F",
        PROVISIONAL,
        f"{DOC}, Table 3 (p. 2)",
        f"V = 1 V, F = 1 MHz, typ; C_T at 0 V not given; {RANGE_NOTE}",
    )
    rp = (
        Provenance(None, "ohm", NOT_AVAILABLE, f"{DOC}: no parallel resistance specified", "OFF modelled lossless")
        if r_p is None
        else Provenance(r_p, "ohm", ASSUMED, "assumed sensitivity value, not device data", "")
    )
    name = "SMP1345-079LF OFF" + ("" if r_p is None else f" (R_p = {r_p:g} ohm assumed)")
    return DiodeModel(name, L_s=_LS.value, C_T=ct.value, R_p=r_p, provenance={"C_T": ct, "L_s": _LS, "R_p": rp})


MACOM_DOC = "MACOM MADP-000907-14020x datasheet Rev V8 (DC-0004063)"
MACOM_RANGE_NOTE = (
    "datasheet: Rs at 10 mA / 1 GHz, CT at -10 V / 1 MHz, 'useable up to 70 GHz' (so 8-12 GHz is inside the stated usable range, but "
    "Rs is specified at 1 GHz and no measured de-embedding exists)"
)


def _madp_ls(l_s):
    """L_s of the flip-chip MADP-000907: not in the datasheet (no package/bump/air-bridge inductance is specified)."""
    if l_s is None:
        return 0.0, Provenance(
            0.0,
            "H",
            NOT_AVAILABLE,
            f"{MACOM_DOC}: no series/package inductance specified",
            "L_s = 0 is a stated ASSUMPTION for the nominal runs; sweep it",
        )
    return float(l_s), Provenance(float(l_s), "H", ASSUMED, "assumed sensitivity value, not device data", "")


def madp_000907_on(l_s=None):
    """MADP-000907-14020 forward-biased: R_s = 5.2 ohm typ (7.0 max), 10 mA, 1 GHz (page 2). `l_s=None`: L_s = 0 (NOT_AVAILABLE)."""
    rs = Provenance(
        5.2,
        "ohm",
        PROVISIONAL,
        f"{MACOM_DOC}, Electrical Specifications (p. 2)",
        f"I_F = 10 mA, F = 1 GHz, typ; {MACOM_RANGE_NOTE}",
        maximum=7.0,
    )
    ls, ls_p = _madp_ls(l_s)
    return DiodeModel("MADP-000907-14020 ON (10mA)", R_s=rs.value, L_s=ls, provenance={"R_s": rs, "L_s": ls_p})


def madp_000907_off(l_s=None, r_p=None):
    """MADP-000907-14020 reverse-biased: C_T = 0.025 pF typ (0.030 max), -10 V, 1 MHz (page 2). R_p is not specified."""
    ct = Provenance(
        0.025e-12,
        "F",
        PROVISIONAL,
        f"{MACOM_DOC}, Electrical Specifications (p. 2)",
        f"V = -10 V, F = 1 MHz, typ; {MACOM_RANGE_NOTE}",
        maximum=0.030e-12,
    )
    rp = (
        Provenance(
            None, "ohm", NOT_AVAILABLE, f"{MACOM_DOC}: no parallel resistance specified", "OFF modelled lossless"
        )
        if r_p is None
        else Provenance(r_p, "ohm", ASSUMED, "assumed sensitivity value, not device data", "")
    )
    ls, ls_p = _madp_ls(l_s)
    return DiodeModel(
        "MADP-000907-14020 OFF", L_s=ls, C_T=ct.value, R_p=r_p, provenance={"C_T": ct, "L_s": ls_p, "R_p": rp}
    )


def diode_load(model, x, y, length, width, direction="x", transverse_conducts=False):
    """A `LumpedLoad` for `model` on a region of `length` (along `direction`) x `width` (across), centred at (x, y) [m].

    The region is the diode/pad footprint: a parameter, not a constant. `transverse_conducts=False` is the physical lumped
    element (R1 finding: a transversely conducting strip is not the empty gap). A model without R_p is a plain series RLC and
    becomes the R1 R/L/C load itself (bit-identical); with R_p the load calls `model.impedance(f)` at every frequency.
    """
    from ..loads import LumpedLoad

    kw = dict(x=x, y=y, length=length, width=width, direction=direction, transverse_conducts=transverse_conducts)
    if model.R_p is None:
        return LumpedLoad(R=model.R_s, L=model.L_s, C=model.C_T, **kw)
    return LumpedLoad(z_of_f=model.impedance, **kw)


def hfss_chain(model, x, y, element_length, island_length, width, direction="x"):
    """The HFSS representation of a series diode circuit as a chain of single-type lumped boundaries (Phase 20 R2 item 3).

    HFSS lumped RLC boundaries are not guaranteed to combine R, L, C in series (the semantics differ by version/setup), so the
    series circuit is built from SEPARATE boundaries, each carrying ONE element type (R, or L, or C, or C parallel R_p -- an
    unambiguous parallel pair), joined end to end by small perfect-conductor islands. Order along `direction`: R_s, C_T(||R_p), L_s
    (elements with value 0 / None are omitted). Returns `(loads, islands)`: `LumpedLoad`s (MoM equivalents of the boundaries,
    `transverse_conducts=False`) and `(cx, cy, length, width)` island rectangles (perfect conductor, in the MoM and in HFSS).
    """

    elems = []
    if model.R_s > 0:
        elems.append(dict(R=model.R_s))
    if model.C_T is not None:
        elems.append(dict(C=model.C_T) if model.R_p is None else dict(z_of_f=_cap_par_r(model.C_T, model.R_p)))
    if model.L_s > 0:
        elems.append(dict(L=model.L_s))
    if not elems:
        raise ValueError("model has no elements (an ideal short needs no boundary)")
    n = len(elems)
    total = n * element_length + (n - 1) * island_length
    pos = -total / 2.0
    loads, islands = [], []
    for k, e in enumerate(elems):
        c = pos + element_length / 2.0
        loads.append(_chain_load(direction, x, y, c, element_length, width, e))
        pos += element_length
        if k < n - 1:
            islands.append(_chain_rect(direction, x, y, pos + island_length / 2.0, island_length, width))
            pos += island_length
    return loads, islands


def _cap_par_r(c, r):
    return lambda f: 1.0 / (1.0 / r + 2j * np.pi * f * c)


def _chain_load(direction, x, y, offset, length, width, elem):
    from ..loads import LumpedLoad

    cx, cy = (x + offset, y) if direction == "x" else (x, y + offset)
    return LumpedLoad(cx, cy, length, width, direction, transverse_conducts=False, **elem)


def _chain_rect(direction, x, y, offset, length, width):
    return (x + offset, y, length, width) if direction == "x" else (x, y + offset, width, length)


def bias_current_budget(n_cells, fraction_on, i_forward=10e-3, v_forward=0.89):
    """Board DC bias implication: (total forward current [A], total DC power [W]) for `n_cells` cells, one diode each,
    `fraction_on` of them ON at `i_forward` (datasheet: V_F = 0.89 V typ at 10 mA, Table 3). Bias-network losses are not included."""
    n_on = n_cells * fraction_on
    return n_on * i_forward, n_on * i_forward * v_forward
