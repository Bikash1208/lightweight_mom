"""Phase 29: seven-hexagon cluster layout and pair-excitation model (MODEL_ONLY): geometry, superposition and classifier properties."""

import numpy as np
import pytest

from src.hologram import seven_hex_cluster as c
from src.hologram import tensor_cp_design as t
from src.hologram.tensor_leaky_cp import SupercellTable


def test_layout_shares_flat_edges_and_outer_azimuths_are_outward():
    cen = c.tile_centres_mm(150.0)
    assert cen[0][:2] == (0.0, 0.0) and cen[0][2] is None
    for x, y, az in cen[1:]:
        assert np.hypot(x, y) == pytest.approx(np.sqrt(3) * 150.0)
        assert np.degrees(np.arctan2(y, x)) % 360.0 == pytest.approx(az % 360.0)
    # neighbours (outer k, k+1) are also sqrt(3) R apart: the seven hexagons tile the plane without overlap or gap
    d = [np.hypot(cen[k][0] - cen[k % 6 + 1][0], cen[k][1] - cen[k % 6 + 1][1]) for k in range(1, 7)]
    assert np.allclose(d, np.sqrt(3) * 150.0)


def test_each_lattice_point_belongs_to_at_most_one_tile():
    cen = c.tile_centres_mm(60.0)
    raw = [set(zip(*(np.round(a, 6) for a in c.tile_cells_mm(q[:2], 60.0)), strict=True)) for q in cen]
    assert any(
        raw[i] & raw[j] for i in range(7) for j in range(i + 1, 7)
    )  # the raw lattice does double-count edge points
    pts = [set(zip(*(np.round(a, 6) for a in xy), strict=True)) for xy in c.owned_cells_mm(cen, 60.0)]
    assert sum(len(p) for p in pts) == len(set().union(*raw))  # ... and ownership loses none
    for i in range(7):
        for j in range(i + 1, 7):
            assert not (pts[i] & pts[j])
    assert all(len(p) > 50 for p in pts)


def test_centre_tile_cells_equal_the_phase26_single_hexagon_lattice():
    x, y = c.tile_cells_mm((0.0, 0.0), 150.0)
    xr, yr = t.lattice_cells(t.hexagon())
    assert x.size == xr.size == 1430
    assert np.allclose(np.sort(x), np.sort(xr * 1e3)) and np.allclose(np.sort(y), np.sort(yr * 1e3))


def test_pair_scan_superposition_identities():
    theta, phi = np.arange(0.0, 90.1, 2.0), np.arange(0.0, 360.0, 4.0)
    rng = np.random.default_rng(1)
    shape = (theta.size, phi.size)
    fa = tuple(rng.normal(size=shape) + 1j * rng.normal(size=shape) for _ in range(2))
    # a copy of the same field excited in phase doubles the amplitude: 6.02 dB over one port, 3.01 dB over the 2 W accepted
    one = c.rhcp_gain_dbi(*fa)
    two = c.rhcp_gain_dbi(2 * fa[0], 2 * fa[1], 2.0)
    assert np.allclose(two - one, 10 * np.log10(2.0))
    # opposite phase cancels
    rows = c.pair_scan(fa, fa, theta, phi, (10.0, 0.0), (10.0, 0.0), phases_deg=[180.0])
    assert rows[0]["peak_gain_dbi"] < one.max() - 60.0


def test_lobes_and_classifier_known_answers():
    theta, phi = np.arange(0.0, 90.1, 1.0), np.arange(0.0, 360.0, 2.0)
    tg, pg = np.meshgrid(theta, phi, indexing="ij")

    def beam(th0, ph0, w=6.0):
        ang = np.degrees(np.arccos(np.clip(np.sum(c._unit(tg, pg) * c._unit(th0, ph0), axis=-1), -1, 1)))
        return -12.0 * (ang / w) ** 2  # Gaussian in dB, -3 dB at w/2

    one = beam(20.0, 90.0)
    lb = c.lobes(one, theta, phi)
    assert len(lb) == 1 and lb[0][1] == pytest.approx(20.0, abs=0.5) and lb[0][2] == pytest.approx(90.0, abs=1.0)
    a, b = (20.0, 90.0), (30.0, 90.0)
    assert c.classify_pair(lb, a, b) == "SINGLE_BEAM_AT_ONE_TILE"
    mid = c.lobes(beam(25.0, 90.0), theta, phi)
    assert c.classify_pair(mid, a, b) == "SINGLE_BEAM_IN_BETWEEN"
    two = np.maximum(beam(20.0, 90.0), beam(60.0, 90.0))
    lt = c.lobes(two, theta, phi)
    assert len(lt) == 2 and c.classify_pair(lt, a, (60.0, 90.0)) == "MULTI_LOBE"


def test_pole_peak_is_reported_with_nonnegative_theta():
    theta, phi = np.arange(0.0, 90.1, 1.0), np.arange(0.0, 360.0, 2.0)
    tg, _ = np.meshgrid(theta, phi, indexing="ij")
    lb = c.lobes(-12.0 * (tg / 6.0) ** 2, theta, phi)
    assert lb[0][1] >= 0.0 and lb[0][1] < 1.0


@pytest.mark.slow
def test_tile_field_power_normalisation_and_tilted_beam_pointing():
    table = SupercellTable.from_files(t.TABLE_GLOB)
    green = c.green_function()
    theta, phi = np.arange(0.0, 90.1, 2.0), np.arange(0.0, 360.0, 4.0)
    tile = c.design_tile(table, (200.0, 0.0), 60.0, 20.0, 0.0, n_psi=360)
    et, ep = c.tile_field(tile, green, theta, phi)
    assert t.prad(et, ep, theta, phi) == pytest.approx(t.LAUNCH * tile["eta_conv"], rel=1e-9)
    lb = c.lobes(c.rhcp_gain_dbi(et, ep), theta, phi)
    assert lb[0][1] == pytest.approx(20.0, abs=3.0) and min(lb[0][2], 360 - lb[0][2]) < 6.0
