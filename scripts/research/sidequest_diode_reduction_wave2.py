"""Second wave of the fewer-diodes study: deeper lambda for the best combinations (to reach 300-500 diodes), appended to diode_reduction.json."""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sidequest_diode_reduction as D
from sidequest_diode_reduction_run import run_spec

OUT = Path("results/sidequest_hybrid/diode_reduction" + D.K.TAG + ".json")

if __name__ == "__main__":
    d = json.loads(OUT.read_text(encoding="utf-8"))
    rows, d5 = d["rows"], d["d5"]
    specs = []
    for name in ("K2", "K4b"):
        dp, dt, ep, _ = d5[name]
        lam0 = D.LAM_BASE[name]
        for m in ("1x2x", "1x2y", "2x2"):
            for f in (1.2, 1.6, 2.0, 2.6, 3.4, 4.5):
                specs.append({"K": name, "id": "C_D2_D3_D5", "D2": m, "D3": "one", "lam": lam0 * f, "dphi": dp, "dtheta": dt, "eps": ep})
                specs.append({"K": name, "id": "C_D2_D5", "D2": m, "lam": lam0 * f, "dphi": dp, "dtheta": dt, "eps": ep})
        for rc in (80.0, 95.0, 115.0):
            for f in (0.8, 1.2, 1.8, 2.6):
                for m in ("1x2y", "2x2"):
                    specs.append({"K": name, "id": "C_D6_D2_D3_D5", "D6": rc, "D3": "one", "D2": m, "lam": lam0 * f, "dphi": dp, "dtheta": dt, "eps": ep})
        for f in (0.15, 0.1, 0.07):
            specs.append({"K": name, "id": "C_D1_D2_D3_D5", "D1": f, "D2": "1x2y", "D3": "one", "lam1": lam0 * 0.8, "dphi": dp, "dtheta": dt, "eps": ep})
            specs.append({"K": name, "id": "C_D1_D2_D3_D5", "D1": f * 2, "D2": "1x2y", "D3": "one", "lam1": lam0 * 0.8, "dphi": dp, "dtheta": dt, "eps": ep})
    t0 = time.time()
    with Pool(4) as pool:
        for r in pool.imap_unordered(run_spec, specs):
            r["id"] = r["spec"].get("id")
            rows.append(r)
    print("wave2 done", len(specs), round(time.time() - t0), flush=True)
    OUT.write_text(json.dumps({"rows": rows, "d5": d5}, indent=1, default=float), encoding="utf-8", newline="\n")
