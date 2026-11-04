"""Random search of the dog-bone gate cell (mesh 32, screening): maximise ON isolation of 3 rows s.t. a matched pass state (OFF S11 <= -18 dB, IL <= 0.3 dB)."""
import json, sys, warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import sidequest_gate_cell as gc

YREF = 1j * 7.82e-3


def job(p):
    warnings.simplefilter("ignore")
    e, s, wd, h, lt, g = p
    out = {"e": e, "s": s, "wd": wd, "h": h, "lt": lt, "g": g}
    for st in ("OFF", "ON"):
        try:
            y, yx = gc.y_xx(gc.make_dogbone(e * 1e-3, s * 1e-3, h * 1e-3, g * 1e-3, st, lt=lt * 1e-3, wd=wd * 1e-3), 32)
            out[st] = {"y": [yx.real * 1e3, yx.imag * 1e3], "r3": gc.strip_transmission(YREF, yx, 3)}
        except Exception as exc:  # noqa: BLE001
            out[st] = {"error": str(exc)[:80]}
    return out


if __name__ == "__main__":
    rng = np.random.default_rng(11)
    pts = []
    q = lambda v: round(float(v), 1)  # noqa: E731  (mm; ALL edges must land on mesh-32 lines = multiples of 0.2 mm)
    for _ in range(120):
        s = q(rng.choice(np.arange(0.8, 4.5, 0.4)))
        wd = q(min(rng.choice([0.8, 1.2, 1.6, 2.0, 2.8, 3.6]), s))
        pts.append((q(rng.choice([0.8, 1.2, 1.6, 2.0, 2.4])), s, wd, q(rng.choice([5.2, 5.6, 6.0])), q(rng.choice([5.6, 6.0])), q(rng.choice([0.4, 0.8]))))
    with ProcessPoolExecutor(4) as ex:
        res = list(ex.map(job, pts))
    Path("results/sidequest_hybrid/gate_search2_mesh32.json").write_text(json.dumps(res, indent=1), encoding="utf-8", newline="\n")
    ok = [r for r in res if "error" not in r["OFF"] and "error" not in r["ON"] and r["OFF"]["r3"] and r["ON"]["r3"]]
    ok.sort(key=lambda r: (r["OFF"]["r3"]["S11_db"] > -18 or r["OFF"]["r3"]["S21_db"] < -0.3, r["ON"]["r3"]["S21_db"]))
    for r in ok[:12]:
        print({k: r[k] for k in ("e", "s", "wd", "h", "lt", "g")}, "OFF b=%.1f S11=%.1f S21=%.2f" % (r["OFF"]["y"][1], r["OFF"]["r3"]["S11_db"], r["OFF"]["r3"]["S21_db"]),
              "| ON y=%.1f%+.1fj S21=%.1f" % (r["ON"]["y"][0], r["ON"]["y"][1], r["ON"]["r3"]["S21_db"]))
