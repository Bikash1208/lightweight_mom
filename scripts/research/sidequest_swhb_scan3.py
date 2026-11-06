"""Board designs on the FINAL refined cell D2 (RO4003C 1.524, p 6.4, gap 0.64, L_ser 1.0 nH).  MODEL_ONLY.   python sidequest_swhb_scan3.py"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sidequest_swhb_scan2 as S2
import sidequest_swhb_run as RUN

if __name__ == "__main__":
    f = "cellB_ypatch_d2a.json"
    jobs = [(0.15, "swap", 3.0, sd, f, 60, 0.64, 6.4e-3, 1.0e-9, f"D2_seed{sd}") for sd in (0, 1, 2)]
    jobs += [(0.15, "add", 3.0, sd, f, 60, 0.64, 6.4e-3, 1.0e-9, f"D2_seed{sd}") for sd in (0, 1)]
    jobs += [(0.125, "swap", 3.0, 1, f, 60, 0.64, 6.4e-3, 1.0e-9, "D2_seed1")]
    with ProcessPoolExecutor(6) as ex:
        res = list(ex.map(S2.one, jobs))
    (RUN.RES / "swhb_scan_D2.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
    for r in res:
        if not r["ok"]:
            print(r)
            continue
        s = r["states"]
        print("%s R %.0f %s diodes %d/%d | S1 G %.1f SLL %.1f edge %.2f Prad %.2f | S2 G %.1f SLL %.1f edge %.2f Prad %.2f | other-beam G %.1f / %.1f | XPD %.0f/%.0f" % (r["tag"], r["R"] * 1e3, r["variant"], r["diodes"], r["n_cells"], s[0]["G_peak_dBi"], s[0]["SLL_db"], s[0]["P_edge"], s[0]["P_rad"], s[1]["G_peak_dBi"], s[1]["SLL_db"], s[1]["P_edge"], s[1]["P_rad"], s[0]["G_at_beam2_dBi"], s[1]["G_at_beam1_dBi"], s[0]["XPD_db"], s[1]["XPD_db"]), flush=True)
