"""Phase 27 item 2: cell-family / state-library registry, described by DOMAIN OF VALIDITY, not by code path.

A `CellFamily` states what a family was characterised for -- substrate, backing, period, frequencies with data, incidence
range, what the data *is* (scalar reflection, 2x2 Jones/tensor, supercell beta/alpha table, discrete switchable states) --
plus its validation label, known biases and fabrication rules. `check_request` answers whether a design request lies inside
that domain; outside it the answer is the named state `NOT_AVAILABLE_FROM_CURRENT_DATA` with the reasons, never an
extrapolation.

Discrete-state libraries (the RIS case, roadmap R4) are `CellFamily` entries with `data_kind="discrete_states"` and one
`StateSpec` per state; each state carries a `ResponseGrid` (a 2x2 reflection matrix versus frequency, theta, phi). No schema
change is needed to register N states (`tests/test_cell_family_registry.py` registers a synthetic 4-state family).
"""

from dataclasses import asdict, dataclass, field
from pathlib import Path

NOT_AVAILABLE = "NOT_AVAILABLE_FROM_CURRENT_DATA"
AVAILABLE = "AVAILABLE"

DATA_KINDS = ("scalar_gamma", "jones_2x2", "supercell_beta_alpha_table", "discrete_states")
RESPONSE_KINDS = ("scalar_gamma", "reflection_2x2", "beta_alpha_table")
VALIDATION_LABELS = (
    "HFSS_VALIDATED",
    "HFSS_VALIDATED_ONE_RUN",
    "HFSS_CALIBRATED",
    "SOLVER_SUPERCELL",
    "SOLVER_INTERNAL",
    "MODEL_ONLY",
    "NOT_HFSS_VALIDATED",
)
_TOL_EPS, _TOL_TAN, _TOL_H_MM, _TOL_F_GHZ, _TOL_P_MM = 1e-6, 1e-6, 1e-4, 5e-3, 1e-6

ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class SubstrateSpec:
    name: str
    eps_r: float
    tan_d: float
    h_mm: float


@dataclass(frozen=True)
class ResponseGrid:
    """What a state's data cover. `reflection_2x2`: a 2x2 Jones matrix (linear x/y basis) per (f, theta, phi)."""

    kind: str
    f_ghz: tuple
    theta_deg: tuple = (0.0,)
    phi_deg: tuple = (0.0,)
    basis: str = "linear_xy"
    data_path: "str | None" = None


@dataclass(frozen=True)
class StateSpec:
    state_id: str
    label: str
    response: ResponseGrid
    bias: "dict | None" = None  # e.g. {"diode_1": "on", "diode_2": "off"} or a geometry description
    note: str = ""


@dataclass(frozen=True)
class CellFamily:
    name: str
    geometry_type: str
    substrate: SubstrateSpec
    backing: str  # "pec" (grounded) | "air"
    period_mm: float
    data_kind: str
    validation_label: str
    frequencies_ghz: tuple = ()  # discrete frequencies with data ...
    frequency_band_ghz: "tuple | None" = None  # ... or a continuous band the solver can build on demand
    incidence_max_theta_deg: float = 0.0  # 0: normal incidence (or the beam's kt) only
    states: tuple = ()
    known_biases: tuple = ()
    fabrication_rules: dict = field(default_factory=dict)
    data_paths: tuple = ()
    validation_evidence: tuple = ()
    notes: str = ""
    diffraction_check: bool = False  # discrete-state RIS families: refuse (f, theta, phi) where a higher diffraction order propagates (Phase 20 R4)

    def to_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class RequestCheck:
    state: str  # AVAILABLE | NOT_AVAILABLE_FROM_CURRENT_DATA
    reasons: tuple
    family: str
    validation_label: str
    known_biases: tuple

    def to_dict(self):
        return asdict(self)


_REGISTRY = {}


def validate_family(f):
    """Schema checks; raises ValueError on a malformed entry. Returns `f`."""
    if f.data_kind not in DATA_KINDS:
        raise ValueError(f"{f.name}: data_kind {f.data_kind!r} not in {DATA_KINDS}")
    if f.validation_label not in VALIDATION_LABELS:
        raise ValueError(f"{f.name}: validation_label {f.validation_label!r} not in {VALIDATION_LABELS}")
    if not f.frequencies_ghz and f.frequency_band_ghz is None:
        raise ValueError(f"{f.name}: needs frequencies_ghz or frequency_band_ghz")
    if f.backing not in ("pec", "air"):
        raise ValueError(f"{f.name}: backing must be 'pec' or 'air'")
    if f.data_kind == "discrete_states":
        if len(f.states) < 2:
            raise ValueError(f"{f.name}: a discrete-state family needs >= 2 states")
        ids = [s.state_id for s in f.states]
        if len(set(ids)) != len(ids):
            raise ValueError(f"{f.name}: duplicate state ids {ids}")
        for s in f.states:
            if s.response.kind not in RESPONSE_KINDS:
                raise ValueError(f"{f.name}/{s.state_id}: response kind {s.response.kind!r}")
            grid_f = set(round(v, 6) for v in s.response.f_ghz)
            need = set(round(v, 6) for v in f.frequencies_ghz)
            if need and not need <= grid_f:
                raise ValueError(
                    f"{f.name}/{s.state_id}: response grid misses family frequencies {sorted(need - grid_f)}"
                )
    elif f.states:
        raise ValueError(f"{f.name}: states only belong to data_kind 'discrete_states'")
    return f


def register(family):
    """Add (or replace) a family; returns it."""
    _REGISTRY[family.name] = validate_family(family)
    return family


def unregister(name):
    _REGISTRY.pop(name, None)


def get(name):
    if name not in _REGISTRY:
        raise KeyError(f"no registered cell family {name!r}")
    return _REGISTRY[name]


def list_families():
    return dict(_REGISTRY)


def check_request(name, substrate=None, f_ghz=None, period_mm=None, theta_deg=None, phi_deg=None):
    """Is a design request inside the family's domain of validity? Unspecified items are not checked.
    `substrate`: `(eps_r, tan_d, h_mm)` or a `SubstrateSpec`."""
    f = get(name)
    why = []
    if substrate is not None:
        e, t, h = (
            (substrate.eps_r, substrate.tan_d, substrate.h_mm) if isinstance(substrate, SubstrateSpec) else substrate
        )
        s = f.substrate
        if abs(e - s.eps_r) > _TOL_EPS or abs(t - s.tan_d) > _TOL_TAN or abs(h - s.h_mm) > _TOL_H_MM:
            why.append(
                f"substrate (eps_r {e:g}, tan_d {t:g}, h {h:g} mm) is not the characterised one (eps_r {s.eps_r:g}, tan_d {s.tan_d:g}, h {s.h_mm:g} mm)"
            )
    if period_mm is not None and abs(period_mm - f.period_mm) > _TOL_P_MM:
        why.append(f"period {period_mm:g} mm is not the characterised {f.period_mm:g} mm")
    if f_ghz is not None:
        ok = any(abs(f_ghz - v) <= _TOL_F_GHZ for v in f.frequencies_ghz)
        if not ok and f.frequency_band_ghz is not None:
            ok = f.frequency_band_ghz[0] <= f_ghz <= f.frequency_band_ghz[1]
        if not ok:
            have = list(f.frequencies_ghz) or list(f.frequency_band_ghz)
            why.append(f"frequency {f_ghz:g} GHz has no data (have {have})")
    if theta_deg is not None and abs(theta_deg) > f.incidence_max_theta_deg + 1e-9:
        why.append(f"incidence {theta_deg:g} deg exceeds the characterised {f.incidence_max_theta_deg:g} deg")
    if f.diffraction_check and f_ghz is not None and theta_deg is not None:
        import numpy as np

        from ..ris.state_library import higher_order_margin

        phis = [phi_deg] if phi_deg is not None else list(range(0, 360, 30))
        bad = [p for p in phis if float(higher_order_margin(f_ghz, abs(theta_deg), p)) < 1.0]
        if bad:
            where = (
                "at the requested azimuth"
                if phi_deg is not None
                else f"at azimuth(s) {[float(np.round(b)) for b in bad]} deg"
            )
            why.append(
                f"a higher diffraction order propagates at {f_ghz:g} GHz, theta {theta_deg:g} deg {where}: the stored (0,0) reflection is not the whole response"
            )
    return RequestCheck(NOT_AVAILABLE if why else AVAILABLE, tuple(why), f.name, f.validation_label, f.known_biases)


def _builtin():
    rt5880 = SubstrateSpec("Rogers RT/duroid 5880", 2.2, 0.0009, 1.575)
    register(
        CellFamily(
            name="tensor_slotted_circle_rt5880_6p4",
            geometry_type="slotted circular patch (D 6.0 mm, rotating rectangular slot L x W)",
            substrate=rt5880,
            backing="pec",
            period_mm=6.4,
            data_kind="supercell_beta_alpha_table",
            validation_label="SOLVER_SUPERCELL",
            frequencies_ghz=(10.54,),
            incidence_max_theta_deg=0.0,
            known_biases=(
                "SLL optimistic by about 4 dB (model -22.8 dB vs HFSS -18.2/-18.5 dB, one hexagon)",
                "AR and XPD pessimistic (model 2.5 dB / 16.9 dB vs HFSS 1.3 dB / 23 dB, one hexagon)",
                "gain anchored to one HFSS run (-0.68 dB offset); conversion eta 0.54 (target 0.85 not reachable on this board)",
                "launch fraction is HFSS_CALIBRATED on one board (pin 6.0 mm above the surface, R x 1.30)",
                "GAMMA_0_COLLINEAR_ONLY: beta/alpha table for propagation along the supercell axis",
            ),
            fabrication_rules={"min_bridge_mm": 0.35, "table_slots_min_bridge_mm": 0.43, "diameter_mm": 6.0},
            data_paths=(
                "results/phase26/supercell_table_4p4x1p2.json",
                "results/phase26/supercell_table_4p6x1p6.json",
                "results/phase26/supercell_table_4p8x1p6.json",
                "results/phase26/supercell_table_4p8x2p0.json",
            ),
            validation_evidence=(
                "hfss_results/hex150_tensor_cp_single/: 23.03 dBi RHCP, AR 1.3 dB, XPD 23 dB, S11 -15.0 dB at 10.54 GHz (one run)",
            ),
            notes="table cells chosen nearest-grid (no cross-period interpolation); other substrates were studied model-only (results/phase26/substrate/)",
        )
    )
    register(
        CellFamily(
            name="onebit_empty_square585_rt5880_6p4",
            geometry_type="1-bit pair: empty cell | 5.85 mm square patch (gap 0.55 mm)",
            substrate=rt5880,
            backing="pec",
            period_mm=6.4,
            data_kind="discrete_states",
            validation_label="NOT_HFSS_VALIDATED",
            frequencies_ghz=(10.54,),
            states=(
                StateSpec(
                    "0",
                    "empty",
                    ResponseGrid("scalar_gamma", (10.54,)),
                    note="reflection phase +138.1 deg at 10.54 GHz",
                ),
                StateSpec(
                    "1",
                    "square 5.85 mm",
                    ResponseGrid("scalar_gamma", (10.54,)),
                    note="reflection phase -41.4 deg at 10.54 GHz",
                ),
            ),
            known_biases=(
                "LOCAL_PERIODICITY_POOR_AT_DESIGN_GAP_PLANE_WAVE_CHECK (Phase 25): isolated-cell model overstates the contrast of the 0.55 mm-gap pair",
                "4-port hex150 board in HFSS: 19.85 dBi vs model 25.8-26.0 dBi; ~40 % grazing power",
            ),
            fabrication_rules={"gap_mm": 0.55},
            data_paths=("src/hologram/auto_states.py",),
            validation_evidence=("hfss_results/hex150_gap0/ (4-port board, port 3 alone): 19.85 dBi at (30, 224)",),
            notes="states are computed on demand by auto_states.py (mesh-corrected MoM); only the 10.54 GHz phases are recorded here",
        )
    )
    register(
        CellFamily(
            name="square_patch_grounded_rt5880_8mm",
            geometry_type="rectangular PEC patch L x W, grounded (Phase 12C design point)",
            substrate=rt5880,
            backing="pec",
            period_mm=8.0,
            data_kind="scalar_gamma",
            validation_label="MODEL_ONLY",
            frequency_band_ghz=(10.0, 15.0),
            known_biases=("near resonance: mesh-convergence residual 3.9 deg unexplained (G4 case, ungrounded)",),
            data_paths=("configs/config_phase12c_grounded_8mm.yaml",),
            notes="library built on demand by the app's auto-library; period chosen in Phase 12B",
        )
    )

    _ris_builtin(rt5880)


def _ris_builtin(rt5880):
    """Phase 20 R3 close-out + R4: the two RIS cell families (MODEL_ONLY, discrete states, 2x2 reflection per state). Data: the R4 state library, full-grid v2 `results/phase20_r4/cp2_library_v2_fullgrid.npz` (the reduced-grid v1 stays selectable in `src/ris/calibrated_library.py`)
    (structure-frame circular reflection of the four one-hot states, tabulated MoM at the Floquet cutoff 2n with a mesh correction, interpolated in (f, theta, phi) with a stored
    uncertainty; see `results/phase20_r4/README.md`). Grid (v2): f 8-12 GHz (0.1 GHz), theta 0-45 deg (5 deg), phi 0-345 deg (15 deg)."""
    freqs = tuple(round(8.0 + 0.1 * k, 6) for k in range(41))
    thetas = tuple(5.0 * k for k in range(10))
    phis = tuple(15.0 * k for k in range(24))
    data = "results/phase20_r4/cp2_library_v2_fullgrid.npz"
    cp_grid = ResponseGrid(
        "reflection_2x2", freqs, theta_deg=thetas, phi_deg=phis, basis="circular_structure_frame", data_path=data
    )
    biases = (
        "MODEL_ONLY: no HFSS number yet (check package hfss_exports/phase20_r3_cp_cell/ pending; comparison/calibration hook src/ris/hfss_calibration.py, calibration only with the user's approval); model uncertainty of the phases ~10 deg (staircased octagon)",
        "PROVISIONAL_DATASHEET: MADP-000907-14020 (MACOM Rev V8) R_s 5.2 ohm typ, C_T 0.025 pF typ (max 0.030 = +20 %); L_s NOT_AVAILABLE_FROM_CURRENT_DATA, 0 assumed",
        "flip-chip part: assembly not confirmed with fabrication (fallback: re-optimise around a packaged low-capacitance PIN diode, C_T <= ~0.05 pF)",
        "VIA_OMITTED_APPROXIMATION; bias network not modelled",
        "opposite-handed (plate-like) reflection about -14 dB at theta 0 / 30 deg against the -15 dB target; min amplitude -1.60 dB at theta 30 (target -1.5)",
        "theta > 30 deg is EXTRAPOLATED_BEYOND_R3_TARGETS (the cell was optimised and checked to 30 deg; 60 deg degrades badly and is not tabulated); points where a higher diffraction order propagates (e.g. theta 45 deg above 11.7 GHz) are refused",
        "isolated-cell predictions are reliable for regular gradients but NOT for random codings at this 15 mm pitch (order levels off 1-4 dB, up to ~10 dB): verify codings at array level (supercell/coupling-aware)",
        "state library v2 (full grid, the registry default) is VALIDATED against the R4 gate: GATE_PASS on 60 random off-grid points (p95 0.31 deg / 0.012 dB) and on 40 points at phi 75-105 deg (p95 2.46 deg / 0.109 dB); isolated outliers up to ~19 deg remain inside the 10.4-10.9 GHz resonance band (results/phase20_r4/validation_v2_a*.json). v1 (reduced grid) stays selectable but FAILS the R4 gate near phi = 90 deg (p95 6.65 deg on the phi 75-105 draw)",
        "the solver needs a Floquet cutoff of 2n for converged results (smaller cutoffs leave isolated numerical lines); R2/R3 numbers used the minimal cutoff and shift by a few degrees (results/phase20_r4/FLOQUET_CONVERGENCE.md)",
    )
    register(
        CellFamily(
            name="ris_cp2_octagon_madp_rt5880_15",
            geometry_type="CP-2: regular octagon (circumradius 4.60 mm) + four 45-degree radial stubs (width 1.375 mm, line 1.748 mm), one MADP-000907 diode each; one-hot state = Pancharatnam-Berry turn 0/45/90/135 deg",
            substrate=rt5880,
            backing="pec",
            period_mm=15.0,
            data_kind="discrete_states",
            validation_label="MODEL_ONLY",
            frequency_band_ghz=(8.0, 12.0),
            incidence_max_theta_deg=45.0,
            states=tuple(
                StateSpec(
                    f"cp{k + 1}",
                    f"diode {k + 1} ON (stub {45 * k} deg), others OFF; r_RR phase ~ {90 * k} deg + const",
                    cp_grid,
                    bias={f"diode_{i + 1}": ("on" if i == k else "off") for i in range(4)},
                    note="handedness-preserving r_RR is the wanted term (propagation-relative circular basis, docs/theory/CONVENTIONS.md); the library stores it referred to the fixed structure axes",
                )
                for k in range(4)
            ),
            known_biases=biases,
            fabrication_rules={
                "diode": "MACOM MADP-000907-14020 (flip-chip)",
                "diodes_per_cell": 4,
                "one_hot": True,
                "period_mm_max_allowed": 15.0,
            },
            data_paths=(
                data,
                "results/phase20_r4/README.md",
                "results/phase20_r3/OPTIMISATION.md",
                "results/phase20_r3/LOCALITY.md",
            ),
            validation_evidence=(),
            notes="CP family; effective stub line is clamped by the cell edge (the optimum sits on this constraint: pitch vs stub length is an open item)",
            diffraction_check=True,
        )
    )
    lp_grid = ResponseGrid(
        "reflection_2x2", freqs, theta_deg=thetas, phi_deg=phis, basis="te_tm_power_wave", data_path=data
    )
    register(
        CellFamily(
            name="ris_lp_1bit_from_cp2_madp_rt5880_15",
            geometry_type="LP interim: 1-bit linear-polarisation codes of the CP-2 cell (LP -> co-LP: states 1 vs 3; LP -> cross-LP: states 2 vs 4)",
            substrate=rt5880,
            backing="pec",
            period_mm=15.0,
            data_kind="discrete_states",
            validation_label="MODEL_ONLY",
            frequency_band_ghz=(8.0, 12.0),
            incidence_max_theta_deg=45.0,
            states=(
                StateSpec("co_0", "LP->co-LP code, bit 0 (CP-2 state 1)", lp_grid, bias={"cp2_state": 1}),
                StateSpec("co_1", "LP->co-LP code, bit 1 (CP-2 state 3)", lp_grid, bias={"cp2_state": 3}),
                StateSpec("cross_0", "LP->cross-LP code, bit 0 (CP-2 state 2)", lp_grid, bias={"cp2_state": 2}),
                StateSpec("cross_1", "LP->cross-LP code, bit 1 (CP-2 state 4)", lp_grid, bias={"cp2_state": 4}),
            ),
            known_biases=biases
            + (
                "1-bit only: a 2-bit LP cell did not meet the targets (LP-1/LP-2/LP-3 with MADP); 2-bit LP is a future phase needing a multilayer / slot-coupled cell",
                "the linear codes are meaningful for an incident polarisation along a structure axis (incidence-plane azimuth phi = 0, 90, 180, 270 deg); see results/phase20_r4/README.md for the code metrics",
            ),
            fabrication_rules={"diode": "MACOM MADP-000907-14020 (flip-chip)", "diodes_per_cell": 4, "one_hot": True},
            data_paths=(data, "results/phase20_r3/LP_INTERIM.md"),
            notes="the same physical cell and the same library as ris_cp2_octagon_madp_rt5880_15: the user selects the polarisation family, not a different geometry",
            diffraction_check=True,
        )
    )


_builtin()
