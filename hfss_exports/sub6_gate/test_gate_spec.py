"""Tests for gate_spec. Pure python: must never import pyaedt or start AEDT."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gate_spec as g


class T(unittest.TestCase):
    def test_no_aedt_import(self):
        self.assertNotIn("ansys", sys.modules)
        self.assertNotIn("pyaedt", sys.modules)

    def test_counts(self):
        rows = g.cases()
        self.assertEqual(len(rows), 22)
        self.assertEqual(sum(r["kind"] == "REF" for r in rows), 2)
        self.assertEqual(sum(r["lid_mm"] == 15.0 for r in rows), 5)
        n3 = [r for r in rows if r["n_gate"] == 3 and r["lid_mm"] == 20.0]
        self.assertEqual(sorted({r["state"] for r in n3}), ["OFF", "OFF_C013", "OFF_C018", "ON"])
        self.assertEqual({r["state"] for r in rows if r["n_gate"] in (1, 2)}, {"OFF", "ON"})

    def test_states(self):
        self.assertEqual(g.STATES["OFF"], (22.23, -335.4))
        self.assertEqual(g.STATES["ON"], (1.60, 29.2))
        self.assertAlmostEqual(g.bar90_sheet(0.13)[1], -390.9, delta=0.15)
        self.assertAlmostEqual(g.bar90_sheet(0.18)[1], -275.0, delta=0.15)

    def test_freqs(self):
        self.assertEqual(len(g.FREQS_GHZ), 13)
        self.assertEqual((g.FREQS_GHZ[0], g.FREQS_GHZ[-1]), (5.5, 6.1))
        self.assertIn(5.8, g.FREQS_GHZ)

    def test_layout(self):
        for n in (1, 2, 3):
            cells, total = g.layout(n)
            self.assertEqual([c["kind"] for c in cells], ["static"] * 3 + ["gate"] * n + ["static"] * 3)
            self.assertAlmostEqual(total, (6 + n) * 11.6 + 2.0)
            self.assertAlmostEqual(cells[0]["xc"] - 5.8 + total / 2, g.SLAB)      # 1.0 mm bare slab to the port
            self.assertAlmostEqual(sum(c["xc"] for c in cells), 0, places=9)       # centred

    def test_dc_gate(self):
        for n in (1, 2, 3):
            r = g.dc_gate(n)
            self.assertTrue(r["ok"], r)
            for c in r["cells"]:
                if c["kind"] == "gate":
                    self.assertTrue(c["island_reaches_ground_via_via"])
                    self.assertFalse(c["disc_reaches_ground_without_diode"])
                    self.assertTrue(c["disc_reaches_ground_with_diode"])
                else:
                    self.assertFalse(c["disc_reaches_ground_with_diode"])

    def test_dc_gate_detects_short(self):
        plain, _ = g.dc_graph(1)
        plain["disc_3"].add("ground")
        plain["ground"].add("disc_3")
        self.assertTrue(g._reach(plain, "disc_3", "ground"))

    def test_clearances(self):
        for l in g.ISLAND_LENGTHS:
            for k, v in g.clearances(l).items():
                self.assertGreater(v, 0.15 - 1e-9, (k, v))

    def test_variables(self):
        for r in g.cases():
            v = g.variables_for(r)
            self.assertEqual(len(v), 16)
            self.assertEqual(v["lid_h"][0], r["lid_mm"])

    def test_new_rows(self):
        nr, sd = g.nr_cases(), g.sd_cases()
        self.assertEqual([r["run_id"] for r in nr], ["NR_D1_slab_H20", "NR_REF_H20", "NR_D2_n2_H20", "NR_D2_n4_H20", "NR_D2_n8_H20"])
        self.assertEqual(len(sd), 24)
        self.assertTrue(all(r["renorm"] is False and r["min_passes"] == 4 for r in nr + sd))
        ids = [r["run_id"] for r in g.all_cases()]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(sorted({r["n_gate"] for r in sd}), [1, 2, 3, 4, 5])
        self.assertEqual({(r["variant"], r["state"]) for r in sd if r["n_gate"] == 3},
                         {(v, st) for v in ("noline", "ylines") for st in ("OFF", "ON", "OFF_C013", "OFF_C018")})
        self.assertEqual(g.row_layout(nr[0])[0], [])
        self.assertAlmostEqual(g.row_layout(nr[0])[1], 2.0)
        self.assertEqual(len(g.variables_for(sd[0])), 14)
        self.assertEqual(g.variables_for(sd[-1])["line_w"][0], 0.1)

    def test_sd_dc(self):
        for v in ("noline", "ylines"):
            for n in (1, 2, 3, 4, 5):
                r = g.sd_dc_gate(n, v)
                self.assertTrue(r["ok"], r)
                self.assertEqual(len(r["pairs"]), n)
        self.assertAlmostEqual(g.line_y0(), 4.382022, places=5)
        self.assertGreater(g.sd_dc_gate(1, "ylines")["clearances_mm"]["cut_to_yline_edge_mm"], 0.15)


if __name__ == "__main__":
    unittest.main()
