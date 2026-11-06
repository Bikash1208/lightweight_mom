"""Refined cell D2 (RO4003C 1.524 mm, p 6.4, gap 0.64 = 0.1 p): patch-grid solves at meshes 40/60/80, gap +-2 px and eps_r +-0.05 at mesh 60.   python sidequest_cellB_d2.py a|b"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sidequest_cellB_ypatch as Y

if __name__ == "__main__":
    lam, h, g0 = "RO4003C", 1.524, 0.64
    gd = 6.4 / 60
    KR = Y.KR
    if sys.argv[1] == "a":
        jobs = [(lam, h, g0, 40, psi, kr, 6.4) for psi in (0.0, 45.0) for kr in KR]
        jobs += [(lam, h, g0, 60, psi, kr, 6.4) for psi in (0.0, 45.0) for kr in KR]
        jobs += [(lam, h, g0 + s * 2 * gd, 60, 0.0, kr, 6.4) for s in (-1, 1) for kr in KR]
        jobs += [(l2, h, g0, 60, 0.0, kr, 6.4) for l2 in ("RO4003C_lo", "RO4003C_hi") for kr in KR]
        Y.run(jobs, ROOT / "results/sidequest_hybrid/cellB_ypatch_d2a.json", workers=8)
    else:
        jobs = [(lam, h, g0, 80, 0.0, kr, 6.4) for kr in KR]
        Y.run(jobs, ROOT / "results/sidequest_hybrid/cellB_ypatch_d2b.json", workers=4)
