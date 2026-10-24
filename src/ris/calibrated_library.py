"""Phase 20 R5-1: the CP-2 state library with the user's HFSS delta, labels, refusals and a diode-tolerance sensitivity band.

`CalibratedLibrary.response(f, theta, phi, state)` returns the structure-frame circular 2x2 response of one state plus honest labels:

* `HFSS_CALIBRATED` - inside the HFSS-measured domain (theta <= 30 deg, phi = 0, f within the HFSS sweep): the user-approved additive delta `C_hfss - C_library` of the baseline cell is added; at the HFSS points the HFSS values are reproduced exactly.
* `HFSS_CALIBRATED_TRANSFERRED_UNCORRECTED_PHI` - phi != 0: the phi = 0 delta is applied as a phi-independent structure-frame offset (an assumption: HFSS was run at phi = 0 only).
* `HFSS_CALIBRATED_EXTRAPOLATED` - theta beyond the HFSS range (delta held at its last value) and/or `EXTRAPOLATED_BEYOND_R3_TARGETS` (library theta > 30 deg).
* `OutOfDomain` (NOT_AVAILABLE_FROM_CURRENT_DATA) outside the library grid or where a higher Floquet order propagates; never silently extrapolated.
* `MODEL_ONLY` always: no HFSS validation of the library itself beyond the calibrated points. Notes flag the F = 2n mesh-correction residual (~12 deg, 10.3-10.8 GHz) and the Wood guard band.

`ToleranceBand` is the sensitivity band for diode tolerance: per-state response differences relative to the nominal C_off from the user's HFSS cases (CSV in the format of `read_hfss_reim` plus a C_off column), added to the calibrated response; `corner_table()` returns the R4 MoM corner table (r_RR only, MODEL_ONLY) as the fallback statement when no HFSS tolerance data exist.
MODEL_ONLY, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.
"""

import csv
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .hfss_calibration import C0, DeltaTable
from .polarization import circular_reflection, from_structure_frame, te_tm_from_circular, to_structure_frame
from .state_library import OutOfDomain, StateLibrary

ROOT = Path(__file__).resolve().parents[2]
LIBRARIES = {  # selectable: v2 = full grid (41 x 10 x 24, server run) with the delta rebuilt against it; v1 = the reduced R4 grid
    "v1": (
        ROOT / "results" / "phase20_r4" / "cp2_library_v1.npz",
        ROOT / "results" / "phase20_r3b" / "hfss_delta_table_run2.npz",
    ),
    "v2": (
        ROOT / "results" / "phase20_r4" / "cp2_library_v2_fullgrid.npz",
        ROOT / "results" / "phase20_r5" / "hfss_delta_table_run2_lib_v2.npz",
    ),
}
DEFAULT_VERSION = "v2"
LIBRARY_STATUS = {  # R4 item 4 gate, server run r4-fullgrid-validate (results/phase20_r4/validation_*.json)
    "v2": "R4_GATE_PASS: validated library (60 random points: p95 0.31 deg / 0.012 dB; 40 points at phi 75-105 deg: p95 2.46 deg / 0.109 dB; limits 3 deg / 0.2 dB, theta <= 30 deg)",
    "v1": "FAILS_R4_GATE_ON_PHI_75_105_DRAW: the reduced-grid library fails the R4 gate near phi = 90 deg (p95 6.65 deg vs 3 deg on 40 points at phi 75-105 deg; it passed the all-phi draw at 1.81 deg); the large errors sit at 10.4-10.9 GHz",
}
DEFAULT_LIBRARY, DEFAULT_DELTA = LIBRARIES[DEFAULT_VERSION]
MESH_RESIDUAL_BAND_GHZ = (10.3, 10.8)  # R4 README: corrected table vs mesh-120 residual up to ~12 deg
NOMINAL_COFF_PF = 0.025


@dataclass(frozen=True)
class Response:
    c: np.ndarray  # (2, 2) structure-frame circular [reflected, incident]
    labels: tuple
    uncertainty: float  # library-only absolute uncertainty (max over the 2x2), NOT the model uncertainty
    notes: tuple = field(default_factory=tuple)


class CalibratedLibrary:
    def __init__(self, library, delta=None, tolerance=None, version=None):
        self.library, self.delta, self.tolerance, self.version = library, delta, tolerance, version

    @classmethod
    def load(cls, library_path=None, delta_path=DEFAULT_DELTA, tolerance=None, version=None):
        """`version` "v1" / "v2" selects a library with its matching delta table; explicit paths override. Default: v2 (full grid)."""
        if version is not None:
            library_path, delta_path = LIBRARIES[version]
        elif library_path is None:
            library_path = DEFAULT_LIBRARY
            if delta_path == DEFAULT_DELTA:
                version = DEFAULT_VERSION
        delta = DeltaTable.load(delta_path) if delta_path is not None else None
        return cls(StateLibrary.load(library_path), delta, tolerance, version)

    def response(self, f_ghz, theta_deg, phi_deg, state, with_tolerance_coff=None):
        """Calibrated response of one state (1..4). Raises `OutOfDomain` outside the library domain. `with_tolerance_coff` (pF) adds the tolerance band's difference to the nominal."""
        c = self.library.structure_circular(f_ghz, theta_deg, phi_deg, state)[0]
        unc = float(self.library.uncertainty(f_ghz, theta_deg, phi_deg, state)[0].max())
        labels, notes = ["MODEL_ONLY", self.library.label(theta_deg)], []
        if self.version == "v1":
            labels.append("LIBRARY_V1_FAILS_R4_GATE_NEAR_PHI_90")
            notes.append(LIBRARY_STATUS["v1"])
        if self.delta is not None:
            val, dlabel = self.delta.at(f_ghz, theta_deg, state)
            c = c + val
            phi0 = abs(((float(phi_deg) + 180.0) % 360.0) - 180.0) < 1e-9
            if dlabel == "HFSS_CALIBRATED" and phi0:
                labels.append("HFSS_CALIBRATED")
            elif dlabel == "HFSS_CALIBRATED":
                labels.append("HFSS_CALIBRATED_TRANSFERRED_UNCORRECTED_PHI")
                notes.append(
                    "the phi = 0 HFSS delta is applied as a phi-independent structure-frame offset; no HFSS data at this phi"
                )
            else:
                labels.append("HFSS_CALIBRATED_EXTRAPOLATED")
                notes.append("outside the HFSS-measured theta / frequency range: delta held at its last value")
        if MESH_RESIDUAL_BAND_GHZ[0] <= float(f_ghz) <= MESH_RESIDUAL_BAND_GHZ[1]:
            notes.append("F = 2n mesh-correction residual up to ~12 deg in this band (R4 validation b)")
        if bool(np.atleast_1d(self.library.guard(f_ghz, theta_deg, phi_deg))[0]):
            notes.append("Wood-anomaly guard band: the response varies faster than the grid resolves")
        if with_tolerance_coff is not None:
            if self.tolerance is None:
                raise OutOfDomain("NOT_AVAILABLE_FROM_CURRENT_DATA: no tolerance band loaded")
            c = c + self.tolerance.difference(f_ghz, theta_deg, state, with_tolerance_coff)
            labels.append("TOLERANCE_BAND_" + self.tolerance.source_label)
        return Response(c, tuple(labels), unc, tuple(notes))

    def all_states(self, f_ghz, theta_deg, phi_deg, with_tolerance_coff=None):
        """(4, 2, 2) responses of the four states and the union of labels."""
        rs = [self.response(f_ghz, theta_deg, phi_deg, s, with_tolerance_coff) for s in (1, 2, 3, 4)]
        labels = tuple(dict.fromkeys(label for r in rs for label in r.labels))
        return np.array([r.c for r in rs]), labels

    def te_tm(self, f_ghz, theta_deg, phi_deg, state):
        r = self.response(f_ghz, theta_deg, phi_deg, state)
        return te_tm_from_circular(from_structure_frame(r.c, math.radians(float(phi_deg)))), r.labels


# --- diode-tolerance band ---------------------------------------------------------------------------------------------------------------------------------------------------
_NUM = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def _reim(row, i, j):
    return float(row[f"re(S(FloquetPort1:{i},FloquetPort1:{j})) []"]) + 1j * float(
        row[f"im(S(FloquetPort1:{i},FloquetPort1:{j})) []"]
    )


def read_hfss_reim_tolerance(path, plane_offset_mm=0.0, coff_column=None):
    """{(coff_pF, state, f_ghz, theta_deg): (2, 2) TE/TM S} from the user's tolerance export: the `read_hfss_reim` columns plus a C_off column (header containing 'coff' or 'c_off', case-insensitive; value in pF, a unit suffix is ignored).
    The reference-plane correction `plane_offset_mm` is applied exactly as in `read_hfss_reim`."""
    out = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for raw in reader:
            r = {k.strip(): v for k, v in raw.items()}
            col = coff_column or next(
                (k for k in r if "coff" in k.lower().replace("_", "") or "c_off" in k.lower()), None
            )
            if col is None:
                raise ValueError("no C_off column found (header containing 'Coff' / 'C_off')")
            coff = float(_NUM.search(r[col]).group())
            if coff > 1.0:  # given in fF? treat values > 1 as fF
                coff /= 1000.0
            th, st, f = float(r["scan_theta [deg]"]), int(float(r["state []"])), float(r["Freq [GHz]"])
            g = np.exp(1j * 2 * (2 * np.pi * f * 1e9 / C0) * (plane_offset_mm * 1e-3) * math.cos(math.radians(th)))
            m = np.array([[_reim(r, 1, 1), _reim(r, 1, 2)], [_reim(r, 2, 1), _reim(r, 2, 2)]]) * g
            out[(round(coff, 6), st, round(f, 6), th)] = m
    return out


def tolerance_plane_offset_mm(tolerance_csv, run1_csv, bare_csv, bare_reference_csv):
    """De-embedding offset for the tolerance export, by the SAME method as the run-2 delta table: run-1 plane (fit of the bare grounded slab against the analytic
    reference) minus the plane shift of the file relative to run 1 measured from the co-polar data (`run_shift_mm`, taken on the nominal C_off rows).
    Returns (offset_mm, shift_mm, (p16, p84))."""
    from .hfss_calibration import fit_reference_plane_mm, read_hfss_reim, run_shift_mm

    d1, _ = fit_reference_plane_mm(bare_csv, bare_reference_csv)
    raw = read_hfss_reim_tolerance(tolerance_csv, 0.0)
    nominal = {(s, f, th, 0.0): m for (c, s, f, th), m in raw.items() if abs(c - NOMINAL_COFF_PF) < 1e-9}
    shift, lo, hi = run_shift_mm(read_hfss_reim(run1_csv, 0.0), nominal)
    return d1 - shift, shift, (lo, hi)


class ToleranceBand:
    """Structure-frame response difference `C(C_off) - C(nominal)` per state on the HFSS (f, theta) grid, interpolated linearly in f, theta and C_off (extrapolation in C_off is refused)."""

    def __init__(self, coffs, f, theta, delta, source_label):
        self.coffs, self.f, self.theta, self.delta, self.source_label = (
            np.asarray(coffs),
            np.asarray(f),
            np.asarray(theta),
            np.asarray(delta),
            source_label,
        )

    @classmethod
    def from_hfss(cls, rows, nominal_coff=NOMINAL_COFF_PF):
        coffs = sorted({k[0] for k in rows})
        if not any(abs(c - nominal_coff) < 1e-9 for c in coffs):
            raise ValueError(f"the tolerance data need the nominal C_off = {nominal_coff} pF case")
        fs = sorted({k[2] for k in rows})
        ths = sorted({k[3] for k in rows})
        delta = np.zeros((len(coffs), len(fs), len(ths), 4, 2, 2), dtype=complex)
        for ic, coff in enumerate(coffs):
            for i, f in enumerate(fs):
                for j, th in enumerate(ths):
                    for s in range(4):
                        c = to_structure_frame(circular_reflection(rows[(coff, s + 1, f, th)]), 0.0)
                        c0 = to_structure_frame(circular_reflection(rows[(nominal_coff, s + 1, f, th)]), 0.0)
                        delta[ic, i, j, s] = c - c0
        return cls(coffs, fs, ths, delta, "HFSS")

    def difference(self, f_ghz, theta_deg, state, coff_pf):
        if coff_pf < self.coffs[0] - 1e-9 or coff_pf > self.coffs[-1] + 1e-9:
            raise OutOfDomain(
                f"NOT_AVAILABLE_FROM_CURRENT_DATA: C_off {coff_pf} pF outside the tolerance data {self.coffs[0]}-{self.coffs[-1]} pF"
            )
        if not (self.f[0] - 1e-9 <= f_ghz <= self.f[-1] + 1e-9):
            raise OutOfDomain("NOT_AVAILABLE_FROM_CURRENT_DATA: frequency outside the tolerance data")
        ci = np.interp(coff_pf, self.coffs, np.arange(len(self.coffs)))
        fi = np.interp(f_ghz, self.f, np.arange(len(self.f)))
        ti = np.interp(min(theta_deg, self.theta[-1]), self.theta, np.arange(len(self.theta)))
        d = self.delta[..., state - 1, :, :]
        out = 0.0
        for a, wa in (
            (int(np.floor(ci)), 1 - (ci - np.floor(ci))),
            (min(int(np.floor(ci)) + 1, len(self.coffs) - 1), ci - np.floor(ci)),
        ):
            for b, wb in (
                (int(np.floor(fi)), 1 - (fi - np.floor(fi))),
                (min(int(np.floor(fi)) + 1, len(self.f) - 1), fi - np.floor(fi)),
            ):
                for c_, wc in (
                    (int(np.floor(ti)), 1 - (ti - np.floor(ti))),
                    (min(int(np.floor(ti)) + 1, len(self.theta) - 1), ti - np.floor(ti)),
                ):
                    out = out + wa * wb * wc * d[a, b, c_]
        return out


def band_edges(cal, f_ghz, theta_deg, phi_deg, state, coffs):
    """Responses of one state at the given C_off values (needs a loaded tolerance band)."""
    return {float(c): cal.response(f_ghz, theta_deg, phi_deg, state, with_tolerance_coff=c).c for c in coffs}


def corner_table(path=ROOT / "results" / "phase20_r4" / "tolerance_corners.json"):
    """The R4 MoM corner rows (MODEL_ONLY): r_RR phase / level shifts per state for C_T +-20 %, L_s 0.1 / 0.2 nH at 9.5 / 10 / 10.5 GHz, theta 0 / 30 - the fallback sensitivity statement when no HFSS tolerance data exist."""
    return json.loads(Path(path).read_text())
