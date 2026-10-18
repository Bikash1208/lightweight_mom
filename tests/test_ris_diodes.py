"""Phase 20 R2 item 1: PIN-diode equivalent circuits -- passivity, limits, resonance, provenance."""

import numpy as np
import pytest

from src.ris.diodes import (
    ASSUMED,
    NOT_AVAILABLE,
    PROVISIONAL,
    DiodeModel,
    bias_current_budget,
    smp1345_off,
    smp1345_on,
)

F = np.linspace(1e9, 20e9, 191)


@pytest.mark.parametrize(
    "m",
    [
        smp1345_on(),
        smp1345_on("1mA"),
        smp1345_off(),
        smp1345_off(3e3),
        smp1345_off(10e3),
        DiodeModel("g", 2.0, 1e-9, 0.2e-12, 5e3),
    ],
)
def test_passive_over_1_to_20_ghz(m):
    assert np.all(m.impedance(F).real >= 0.0)


def test_limits_give_short_and_open():
    short = DiodeModel("s", R_s=0.0, L_s=0.0).impedance(10e9)
    assert short == 0
    assert abs(DiodeModel("c", C_T=1.0).impedance(10e9)) < 1e-9  # C_T -> inf: short
    assert abs(DiodeModel("c", C_T=1e-21).impedance(10e9)) > 1e7  # C_T -> 0: open
    assert abs(smp1345_off().scaled(C_T=1e-6).impedance(10e9)) > 1e7


def test_off_series_resonance_equals_one_over_two_pi_sqrt_LC_and_is_zero_reactance():
    m = smp1345_off()
    f0 = 1 / (2 * np.pi * np.sqrt(0.7e-9 * 0.19e-12))
    assert m.series_resonance_hz == pytest.approx(f0, rel=1e-12)
    assert abs(m.impedance(f0).imag) < 1e-9 * abs(m.impedance(0.5 * f0))
    assert m.impedance(0.8 * f0).imag < 0 < m.impedance(1.2 * f0).imag  # capacitive below, inductive above
    assert smp1345_on().series_resonance_hz is None


def test_on_state_is_Rs_plus_jwLs():
    z = smp1345_on().impedance(10e9)
    assert z.real == 1.5 and z.imag == pytest.approx(2 * np.pi * 10e9 * 0.7e-9, rel=1e-12)
    assert smp1345_on("1mA").impedance(10e9).real == 3.5


def test_rp_in_parallel_with_C_and_lossless_default():
    m = smp1345_off(3e3)
    w = 2 * np.pi * 10e9
    want = 1j * w * 0.7e-9 + 1 / (1 / 3e3 + 1j * w * 0.19e-12)
    assert m.impedance(10e9) == pytest.approx(want, rel=1e-12)
    assert smp1345_off().impedance(10e9).real == 0.0
    with pytest.raises(ValueError):
        DiodeModel("x", R_p=1e3)


def test_provenance_labels_and_scaling():
    on, off = smp1345_on(), smp1345_off()
    assert all(p.status == PROVISIONAL for p in on.provenance.values())
    assert off.provenance["R_p"].status == NOT_AVAILABLE and off.provenance["R_p"].value is None
    assert smp1345_off(3e3).provenance["R_p"].status == ASSUMED
    assert "extrapolation" in on.provenance["R_s"].condition
    s = off.scaled(C_T=1.2)
    assert pytest.approx(0.19e-12 * 1.2) == s.C_T and s.provenance["C_T"].status == "PERTURBED_SENSITIVITY"
    assert off.C_T == 0.19e-12  # original untouched
    with pytest.raises(ValueError):
        off.scaled(R_p=2.0)


def test_input_validation():
    for bad in (dict(R_s=-1.0), dict(L_s=float("nan")), dict(C_T=0.0), dict(C_T=-1e-12)):
        with pytest.raises(ValueError):
            DiodeModel("b", **bad)
    with pytest.raises(ValueError):
        smp1345_on().impedance(0.0)
    with pytest.raises(ValueError):
        smp1345_on("5mA")


def test_bias_budget():
    i, p = bias_current_budget(5685, 0.5)
    assert i == pytest.approx(5685 * 0.5 * 10e-3) and p == pytest.approx(i * 0.89)


# ---- item 2: diode -> LumpedLoad ----
import math  # noqa: E402
import warnings  # noqa: E402

from src.basis_functions import MeshSettings  # noqa: E402
from src.composite_geometry import AxisRect, CompositeGeometry  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.geometry import UnitCell  # noqa: E402
from src.loads import LumpedLoad  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402
from src.ris.diodes import diode_load  # noqa: E402
from src.sparameters import compute_s_parameters_oblique, sweep_s_parameters_oblique  # noqa: E402

SUB = Substrate("rt5880", 2.2, 0.0, 1.575e-3)
PLANES = ReferencePlanes(2e-3, None, SUB.thickness)
MESH = MeshSettings(20, 20, 10)


def _cell(loads):
    # two arms joined only by the diode bridge (1 mm along x, 1 mm across): the load is the only conduction path across the gap
    geo = CompositeGeometry(
        metal=(
            AxisRect(-2.25e-3, 0.0, 3.5e-3, 6e-3),
            AxisRect(2.25e-3, 0.0, 3.5e-3, 6e-3),
            AxisRect(0.0, 0.0, 1e-3, 1e-3),
        )
    )
    return UnitCell(RectangularLattice(10e-3, 10e-3), geo, SUB, Metal(), backing="pec", loads=tuple(loads))


def _s(cell, f=10e9, theta=0.0, phi=0.0):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return compute_s_parameters_oblique(cell, f, MESH, PLANES, theta, phi)


def test_series_rlc_diode_is_exactly_the_r1_load():
    m = DiodeModel("c", 2.0, 0.7e-9, 0.19e-12)
    ld = diode_load(m, 0, 0, 1e-3, 1e-3)
    assert (ld.R, ld.L, ld.C, ld.z_of_f, ld.transverse_conducts) == (2.0, 0.7e-9, 0.19e-12, None, False)
    r1 = LumpedLoad(0, 0, 1e-3, 1e-3, "x", R=2.0, L=0.7e-9, C=0.19e-12, transverse_conducts=False)
    assert ld == r1
    assert np.array_equal(
        _s(_cell([ld]), 10e9, math.radians(30), 0.2).s11, _s(_cell([r1]), 10e9, math.radians(30), 0.2).s11
    )


def test_parallel_branch_uses_model_impedance_at_each_frequency():
    m = smp1345_off(3e3)
    ld = diode_load(m, 0, 0, 1e-3, 1e-3)
    for f in (8e9, 10e9, 12e9):
        assert ld.impedance(2 * math.pi * f) == pytest.approx(m.impedance(f), rel=1e-12)
    # a constant z_of_f reproduces the equivalent R1 load
    const = LumpedLoad(0, 0, 1e-3, 1e-3, "x", z_of_f=lambda f: 5.0 + 0j, transverse_conducts=False)
    plain = LumpedLoad(0, 0, 1e-3, 1e-3, "x", R=5.0, transverse_conducts=False)
    assert np.array_equal(_s(_cell([const])).s11, _s(_cell([plain])).s11)


def test_sweep_8_to_12_ghz_runs_the_load_per_frequency_and_on_off_differ():
    freqs = [8e9, 10e9, 12e9]
    on, _ = sweep_s_parameters_oblique(_cell([diode_load(smp1345_on(), 0, 0, 1e-3, 1e-3)]), freqs, MESH, PLANES, 0.0)
    off, _ = sweep_s_parameters_oblique(_cell([diode_load(smp1345_off(), 0, 0, 1e-3, 1e-3)]), freqs, MESH, PLANES, 0.0)
    for f, a, b in zip(freqs, on, off, strict=True):
        single = _s(_cell([diode_load(smp1345_on(), 0, 0, 1e-3, 1e-3)]), f)
        assert np.array_equal(a.s11, single.s11)  # the sweep evaluates Z_L(f) at each frequency
        assert abs(a.s11[1, 1] - b.s11[1, 1]) > 0.05  # TM = x, across the x-directed gap


# ---- item 3: HFSS chain representation vs the ideal series model ----
import importlib.util  # noqa: E402
from pathlib import Path  # noqa: E402

from src.ris.diodes import hfss_chain  # noqa: E402


def _repr_script():
    p = Path(__file__).resolve().parents[1] / "scripts" / "research" / "phase20_r2_hfss_representation.py"
    spec = importlib.util.spec_from_file_location("phase20_r2_hfss_representation", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_chain_structure_and_total_impedance():
    for m in (smp1345_on(), smp1345_off(), smp1345_off(3e3)):
        loads, islands = hfss_chain(m, 0.0, 0.0, 0.25e-3, 0.5e-3, 0.5e-3)
        assert len(loads) == 2 and len(islands) == 1 and not any(ld.transverse_conducts for ld in loads)
        z = sum(ld.impedance(2 * np.pi * 10e9) for ld in loads)
        assert z == pytest.approx(m.impedance(10e9), rel=1e-12)  # the series sum of the boundaries is the diode's Z_L
        xs = sorted(ld.x for ld in loads)
        assert xs[0] == pytest.approx(-0.375e-3) and xs[1] == pytest.approx(0.375e-3)
        assert islands[0][0] == pytest.approx(0.0) and islands[0][2] == pytest.approx(0.5e-3)
    with pytest.raises(ValueError):
        hfss_chain(DiodeModel("empty"), 0, 0, 1e-3, 1e-3, 1e-3)


def test_hfss_chain_equals_ideal_series_model_with_pads_and_not_without():
    mod = _repr_script()
    with_pad = mod.compare(0.5e-3, 40)
    no_pad = mod.compare(0.0, 40)
    on = [r for r in with_pad if r["model"].startswith("ON")]
    off = [r for r in with_pad if r["model"].startswith("OFF")]
    assert max(abs(r["d_phase_deg"]) for r in on) < 0.5 and max(abs(r["d_abs"]) for r in on) < 1e-3
    # OFF: a systematic ~6 deg from where the C and L sit along the span (not a mesh effect: 40 and 80 agree to < 1 deg, see results/phase20_r2)
    assert max(abs(r["d_phase_deg"]) for r in off) < 8.0 and max(abs(r["d_abs"]) for r in off) < 0.02
    # elements touching the wide arms are a different (2-D current crowding) problem: the pads are required
    assert max(abs(r["d_phase_deg"]) for r in no_pad) > 5 * max(abs(r["d_phase_deg"]) for r in with_pad)
