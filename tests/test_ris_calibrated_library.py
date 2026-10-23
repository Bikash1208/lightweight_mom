"""R5-1: calibrated library layer - exact HFSS reproduction at the HFSS points, labels and refusals, tolerance band reader / interpolation."""

import csv
from pathlib import Path

import numpy as np
import pytest

from src.ris.calibrated_library import (
    CalibratedLibrary,
    ToleranceBand,
    band_edges,
    corner_table,
    read_hfss_reim_tolerance,
)
from src.ris.hfss_calibration import DeltaTable, fit_reference_plane_mm, read_hfss_reim
from src.ris.polarization import circular_reflection, to_structure_frame
from src.ris.state_library import OutOfDomain, StateLibrary

ROOT = Path(__file__).resolve().parents[1]
HF = ROOT / "hfss_results" / "phase20_r3_cp_cell"
RUN1 = HF / "hfss_cp2_states1-4_theta0_30.csv"


@pytest.fixture(scope="module")
def cal():
    lib = StateLibrary.load(ROOT / "results/phase20_r4/cp2_library_v1.npz")
    d, _ = fit_reference_plane_mm(
        HF / "hfss_bare_slab_theta0_30.csv", ROOT / "hfss_exports/phase20_r3_cp_cell/bare_grounded_slab_reference.csv"
    )
    h = read_hfss_reim(RUN1, d)
    t = DeltaTable.build(lib, h, "run 1 (test)", approved_by_user=True)
    return CalibratedLibrary(lib, t), h, d


def test_hfss_points_are_reproduced_exactly_and_labelled(cal):
    c, h, _ = cal
    for st, f, th in [(1, 10.5, 0.0), (3, 10.0, 30.0), (4, 9.0, 0.0)]:
        r = c.response(f, th, 0.0, st)
        want = to_structure_frame(circular_reflection(h[(st, f, th, 0.0)]), 0.0)
        np.testing.assert_allclose(r.c, want, atol=1e-9)
        assert "HFSS_CALIBRATED" in r.labels and "MODEL_ONLY" in r.labels


def test_labels_for_phi_theta_and_notes(cal):
    c = cal[0]
    assert "HFSS_CALIBRATED_TRANSFERRED_UNCORRECTED_PHI" in c.response(10.5, 0.0, 45.0, 1).labels
    r = c.response(10.5, 40.0, 0.0, 2)
    assert "HFSS_CALIBRATED_EXTRAPOLATED" in r.labels and "EXTRAPOLATED_BEYOND_R3_TARGETS" in r.labels
    assert any("mesh-correction residual" in n for n in c.response(10.5, 0.0, 0.0, 1).notes)
    assert not any("mesh-correction residual" in n for n in c.response(9.0, 0.0, 0.0, 1).notes)
    cs, labels = c.all_states(10.5, 0.0, 0.0)
    assert cs.shape == (4, 2, 2) and "HFSS_CALIBRATED" in labels
    te_tm, _ = c.te_tm(10.5, 30.0, 0.0, 1)
    assert te_tm.shape == (2, 2)


def test_refusals_outside_the_library_domain(cal):
    c = cal[0]
    with pytest.raises(OutOfDomain, match="NOT_AVAILABLE_FROM_CURRENT_DATA"):
        c.response(7.0, 0.0, 0.0, 1)  # frequency outside
    with pytest.raises(OutOfDomain, match="NOT_AVAILABLE_FROM_CURRENT_DATA"):
        c.response(10.5, 50.0, 0.0, 1)  # theta outside
    with pytest.raises(OutOfDomain, match="NOT_AVAILABLE_FROM_CURRENT_DATA"):
        c.response(12.0, 45.0, 0.0, 1)  # a higher order propagates
    with pytest.raises(OutOfDomain):
        c.response(10.5, 0.0, 0.0, 1, with_tolerance_coff=0.0225)  # no band loaded


def _write_tolerance_csv(path, h_rows, d):
    """Synthetic C_off export in the user's format (run-1 columns + a Coff column): C_off = 0.025 is the HFSS value, the others multiply every S by a known smooth factor."""
    with open(RUN1, newline="") as fh:
        reader = csv.DictReader(fh)
        header = [k for k in reader.fieldnames] + ["Coff [pF]"]
        base = list(reader)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        w.writeheader()
        for coff, tag in ((0.0225, "0.0225pF"), (0.025, "0.025"), (0.030, "0.030")):
            for r in base:
                row = dict(r)
                g = np.exp(1j * np.radians(400.0 * (coff - 0.025)))
                for i in (1, 2):
                    for j in (1, 2):
                        re_k, im_k = (
                            f"re(S(FloquetPort1:{i},FloquetPort1:{j})) []",
                            f"im(S(FloquetPort1:{i},FloquetPort1:{j})) []",
                        )
                        z = (float(r[re_k]) + 1j * float(r[im_k])) * g
                        row[re_k], row[im_k] = repr(float(z.real)), repr(float(z.imag))
                row["Coff [pF]"] = tag
                w.writerow(row)


def test_tolerance_reader_band_and_interpolation(cal, tmp_path):
    c, _, d = cal
    p = tmp_path / "tol.csv"
    _write_tolerance_csv(p, None, d)
    rows = read_hfss_reim_tolerance(p, d)
    assert {k[0] for k in rows} == {0.0225, 0.025, 0.03}
    band = ToleranceBand.from_hfss(rows)
    tc = CalibratedLibrary(c.library, c.delta, band)
    nominal = tc.response(10.5, 0.0, 0.0, 1).c
    np.testing.assert_allclose(tc.response(10.5, 0.0, 0.0, 1, with_tolerance_coff=0.025).c, nominal, atol=1e-12)
    # a known phase factor exp(j*400*(coff-0.025) deg) on every S rotates the TE/TM matrix by a common phase: check on the 0.0225 case against the direct construction
    got = tc.response(10.5, 0.0, 0.0, 1, with_tolerance_coff=0.0225)
    assert "TOLERANCE_BAND_HFSS" in got.labels
    g = np.exp(1j * np.radians(400.0 * (0.0225 - 0.025)))
    ref = to_structure_frame(circular_reflection(rows[(0.025, 1, 10.5, 0.0)] * g), 0.0) - to_structure_frame(
        circular_reflection(rows[(0.025, 1, 10.5, 0.0)]), 0.0
    )
    np.testing.assert_allclose(got.c - nominal, ref, atol=1e-9)
    half = band.difference(10.5, 0.0, 1, 0.02375)  # midway between 0.0225 and 0.025
    np.testing.assert_allclose(half, 0.5 * band.difference(10.5, 0.0, 1, 0.0225), atol=1e-12)
    with pytest.raises(OutOfDomain, match="C_off"):
        band.difference(10.5, 0.0, 1, 0.04)
    edges = band_edges(tc, 10.5, 0.0, 0.0, 2, (0.0225, 0.025, 0.03))
    assert set(edges) == {0.0225, 0.025, 0.03}


def test_tolerance_needs_the_nominal_case(cal, tmp_path):
    _, _, d = cal
    p = tmp_path / "tol.csv"
    _write_tolerance_csv(p, None, d)
    rows = {k: v for k, v in read_hfss_reim_tolerance(p, d).items() if k[0] != 0.025}
    with pytest.raises(ValueError, match="nominal"):
        ToleranceBand.from_hfss(rows)


def test_corner_table_is_the_r4_mom_fallback():
    rows = corner_table()
    assert len(rows) > 0 and "corner" in rows[0]


def test_tolerance_compare_script_reports_zero_difference_for_consistent_inputs(tmp_path):
    """Synthetic consistent inputs: HFSS rows = the library's corrected TE/TM response, 'MoM' array = the library's raw response (the script adds the library correction): all differences ~0."""
    import sys

    from src.ris.polarization import from_structure_frame, te_tm_from_circular

    sys.path.insert(0, str(ROOT / "scripts" / "research"))
    import phase20_r5_tolerance_compare as cmp

    lib = StateLibrary.load(ROOT / "results/phase20_r4/cp2_library_v1.npz")
    f, ths, coffs = 10.5, (0.0, 30.0), (0.0225, 0.025, 0.030)
    c_struct = np.zeros((1, 3, 1, 2, 4, 2, 2), dtype=complex)
    header = [
        "scan_theta [deg]", "state []", "Freq [GHz]", "Coff [pF]",
        "re(S(FloquetPort1:1,FloquetPort1:1)) []", "im(S(FloquetPort1:1,FloquetPort1:1)) []",
        "re(S(FloquetPort1:1,FloquetPort1:2)) []", "im(S(FloquetPort1:1,FloquetPort1:2)) []",
        "re(S(FloquetPort1:2,FloquetPort1:1)) []", "im(S(FloquetPort1:2,FloquetPort1:1)) []",
        "re(S(FloquetPort1:2,FloquetPort1:2)) []", "im(S(FloquetPort1:2,FloquetPort1:2)) []",
    ]  # fmt: skip
    p = tmp_path / "tol.csv"
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for ic, coff in enumerate(coffs):
            for jt, th in enumerate(ths):
                for s in range(4):
                    c_struct[0, ic, 0, jt, s] = lib.structure_circular(f, th, 0.0, s + 1, key="c_raw")[0]
                    m = te_tm_from_circular(
                        from_structure_frame(lib.structure_circular(f, th, 0.0, s + 1, key="c_corr")[0], 0.0)
                    )
                    w.writerow(
                        [
                            th,
                            s + 1,
                            f,
                            coff,
                            m[0, 0].real,
                            m[0, 0].imag,
                            m[0, 1].real,
                            m[0, 1].imag,
                            m[1, 0].real,
                            m[1, 0].imag,
                            m[1, 1].real,
                            m[1, 1].imag,
                        ]
                    )
    npz = tmp_path / "mom.npz"
    np.savez(
        npz,
        c_structure=c_struct,
        coffs_pf=np.array(coffs),
        f_ghz=np.array([f]),
        theta_deg=np.array(ths),
        meshes=np.array([90]),
    )
    rep = cmp.main(p, npz, 0.0, tmp_path / "r.json")
    assert rep["at_10p5_ghz"]["n"] == 24 and rep["at_10p5_ghz"]["max_abs_phase_diff_deg"] < 1e-6


def test_run_shift_recovers_a_known_plane_move():
    import math

    import numpy as np

    from src.ris.hfss_calibration import C0, run_shift_mm

    rng = np.random.default_rng(1)
    a, b = {}, {}
    for f in (9.0, 10.0, 11.0):
        for th in (0.0, 30.0):
            m = np.diag(np.exp(1j * rng.uniform(-3, 3, 2)) * 0.97).astype(complex)
            a[(1, f, th, 0.0)] = m
            g = np.exp(-1j * 2 * (2 * np.pi * f * 1e9 / C0) * 1.7e-3 * math.cos(math.radians(th)))
            b[(1, f, th, 0.0)] = m * g
    assert run_shift_mm(a, b)[0] == pytest.approx(-1.7, abs=1e-9)
