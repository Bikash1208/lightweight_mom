"""Four-state tables (OFF/OFF, ON/OFF, OFF/ON, ON/ON) of the two chosen switched cells: Bloch Y tensor (propagation along x for the x row, along y for the y row) at the cell's own beta, dY per axis, alpha."""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sidequest_cell_search as C

CASES = [("ro3010", dict(kind="notch", g=0.2, window=0.4, wd=0.4)), ("rt5880", dict(diam=6.0, ls=2.8, pos=0.8))]

if __name__ == "__main__":
    with ProcessPoolExecutor(2) as ex:
        res = list(ex.map(C.job, [(s, kw, 64, None, "full") for s, kw in CASES]))
    for r in res:
        print(r["stack"], r["kw"], r["analysis"], flush=True)
        for k, v in r["Y_mS"].items():
            print("  ", k, [f"{a:.2f}{b:+.2f}j" for a, b in v], flush=True)
    Path("results/sidequest_hybrid/cell_search_full.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
