"""MA4AGP907 measured diode (Phase 28): the extraction reproduces the advisor's values within 2 % from the vendor files (skipped when they are not on this machine); the committed extracted table gives the same values
and the measured-diode load is passive and frequency dependent."""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))
REFS = Path("C:/hmsa_sync/refs/diodes")
ADVISOR_R = {"ON_5mA": 3.3, "ON_10mA": 2.4, "ON_15mA": 2.0}
ADVISOR_C_FF = {"OFF_0V": 28.2, "OFF_m5V": 27.6}


def _check(doc):
    for k, r in ADVISOR_R.items():
        assert doc["states"][k]["at_10_ghz"]["z_ohm"][0] == pytest.approx(r, rel=0.02)
        assert -6.1 <= min(doc["states"][k]["z_im_ohm"]) and max(doc["states"][k]["z_im_ohm"]) <= -1.8  # reactance -1.9 .. -5.9 ohm over 4 - 20 GHz (no series inductance)
    assert doc["states"]["ON_10mA"]["at_10_ghz"]["z_ohm"][1] == pytest.approx(-3.7, rel=0.02)
    for k, c in ADVISOR_C_FF.items():
        assert doc["states"][k]["at_10_ghz"]["c_ff"] == pytest.approx(c, rel=0.02)
        assert doc["states"][k]["c_ff_mean_4_20_ghz"] == pytest.approx(c, rel=0.02)  # flat 4 - 20 GHz


@pytest.mark.skipif(not REFS.exists(), reason="vendor S-parameter files not on this machine")
def test_extraction_from_the_vendor_files_reproduces_the_advisor_values():
    import phase28_diode_extract as ex

    _check(ex.extract(REFS))


def test_committed_table_gives_the_same_values_and_a_passive_load():
    import json

    from src.ris import diodes_measured as dm

    _check(json.loads(dm.DATA.read_text(encoding="utf-8")))
    for s in ("ON_5mA", "ON_10mA", "OFF_0V"):
        z = dm.z_of_f(s)([8e9, 10e9, 12e9])
        assert (z.real >= 0).all()
    on, off = dm.z_of_f("ON_10mA")(10e9), dm.z_of_f("OFF_0V")(10e9)
    assert on.real == pytest.approx(2.43, rel=0.02) and off.imag == pytest.approx(-564.3, rel=0.01)
    assert dm.z_of_f("OFF_0V", {"C": 0.8})(10e9).imag == pytest.approx(off.imag / 0.8)
    assert dm.provenance("ON_10mA")["label"] == "MEASURED_VENDOR_SPARAM"
