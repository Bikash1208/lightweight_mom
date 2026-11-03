"""Phase 30 S1 (MODEL_ONLY): direct MoM leaky-mode solve of the 7-cell rotating-axis supercell (locality test).

    python scripts/research/phase30_s1_supercell_mom.py <mesh_per_cell> [out_dir]

RO3010 0.635 mm, eps_r 11.2, E1 cell (3.2 mm pitch, super-ellipse n = 4 major 3.05 mm, slot 1.89 x 0.58 mm, outline fixed, slot turned by
180*i/7 deg), grounded (PEC at z = -h), Floquet order 2 x mesh per cell. Stages (all at ky0 = 0):

  A  f-mode, tan d = 0,      kx0 = KX_A (HFSS 0.4 deg/supercell -> 0.31 rad/m): complex f -> Q_rad
  B  f-mode, tan d = 0.0022, kx0 = KX_A: Q_total
  C  f-mode, tan d = 0,      kx0 = KX_C (HFSS 6.74 deg/supercell): f(kx0) -> v_g = 2 pi df / dkx
  D  k-mode, tan d = 0,      f = Re f_A: complex kx0 -> alpha_rad directly
  E  k-mode, tan d = 0.0022, f = Re f_A: alpha_total

Writes <out_dir>/s1_mom_m<mesh>.json (resumable: finished stages are kept) and a .log.
"""

import json
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


from src.dispersion import leaky_supercell_mom as L  # noqa: E402
from src.frequency import free_space_wavenumber  # noqa: E402
from src.materials import Substrate  # noqa: E402

EPS, H, TAND, N = 11.2, 0.635e-3, 0.0022, 7
LX = N * 3.2e-3
KX_A = math.radians(0.4) / LX
KX_C = math.radians(6.74) / LX
C0 = 299792458.0


def main(mesh, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    jpath, lpath = out / f"s1_mom_m{mesh}.json", out / f"s1_mom_m{mesh}.log"
    res = json.loads(jpath.read_text()) if jpath.exists() else {}

    def log(s):
        line = time.strftime("%H:%M:%S ") + s
        print(line, flush=True)
        with open(lpath, "a", encoding="utf-8", newline="\n") as fh:
            fh.write(line + "\n")

    subs = {t: Substrate("ro3010", EPS, t, H) for t in (0.0, TAND)}
    cell, basis = L.build_supercell(L.s1_angles(N), mesh, subs[0.0])
    ox, oy = L.orders_for(basis, N, 2)
    probe = L.harmonic_probe(basis, KX_A, 1)
    log(f"mesh {mesh}: {basis.n_unknowns} unknowns, orders {ox} x {oy}")

    def fmode(tag, tand, kx0, guess):
        if tag in res:
            return complex(*res[tag]["f_ghz"])
        sub = subs[tand]
        t0 = time.time()
        p, _ = L.find_mode(
            lambda p: L.assemble_z_complex(basis, sub, p * 1e9, kx0, 0.0, ox, oy),
            guess,
            guess * 1.01 + 0.01j,
            basis.n_unknowns,
            log=log,
            probe=L.harmonic_probe(basis, kx0, 1),
        )
        res[tag] = {
            "f_ghz": [p.real, p.imag],
            "Q": p.real / (2 * p.imag),
            "tand": tand,
            "kx0": kx0,
            "seconds": time.time() - t0,
        }
        jpath.write_text(json.dumps(res, indent=1))
        log(f"{tag}: f = {p} GHz, Q = {p.real / (2 * p.imag):.2f}")
        return p

    def kmode(tag, tand, f_ghz, guess):
        if tag in res:
            return complex(*res[tag]["kx0"])
        sub = subs[tand]
        t0 = time.time()
        p, _ = L.find_mode(
            lambda p: L.assemble_z_complex(basis, sub, f_ghz * 1e9, p, 0.0, ox, oy),
            guess,
            guess + 0.5 - 0.5j,
            basis.n_unknowns,
            log=log,
            probe=probe,
        )
        k0 = free_space_wavenumber(f_ghz * 1e9)
        res[tag] = {
            "kx0": [p.real, p.imag],
            "alpha_np_m": -p.imag,
            "alpha_over_k0": -p.imag / k0,
            "f_ghz": f_ghz,
            "tand": tand,
            "seconds": time.time() - t0,
        }
        jpath.write_text(json.dumps(res, indent=1))
        log(f"{tag}: kx0 = {p} rad/m, alpha/k0 = {-p.imag / k0:.5f}")
        return p

    guess = complex(sys.argv[3]) if len(sys.argv) > 3 else 10.2 + 0.14j
    fa = fmode("A_rad_f", 0.0, KX_A, guess)
    fmode("B_tot_f", TAND, KX_A, fa)
    fc = fmode("C_rad_f_kx6p74", 0.0, KX_C, fa)
    vg = abs(
        2 * math.pi * (fc.real - fa.real) * 1e9 / (KX_C - KX_A)
    )  # sign: the mode runs toward -x (df/dkx0 < 0); magnitude only
    k0 = free_space_wavenumber(fa.real * 1e9)
    res["derived"] = {
        "vg_over_c": vg / C0,
        "alpha_over_k0_from_Q_rad": C0 / (2 * res["A_rad_f"]["Q"] * vg) if vg else None,
        "alpha_over_k0_from_Q_tot": C0 / (2 * res["B_tot_f"]["Q"] * vg) if vg else None,
        "Q_diel": 1.0 / (1.0 / res["B_tot_f"]["Q"] - 1.0 / res["A_rad_f"]["Q"]),
        "k0_at_f_rad": k0,
    }
    jpath.write_text(json.dumps(res, indent=1))
    log(f"derived: {res['derived']}")
    kmode("D_rad_k", 0.0, fa.real, complex(KX_A, -res["derived"]["alpha_over_k0_from_Q_rad"] * k0))
    kmode("E_tot_k", TAND, fa.real, complex(KX_A, -res["derived"]["alpha_over_k0_from_Q_tot"] * k0))
    log("done")


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2] if len(sys.argv) > 2 else str(ROOT / "results" / "phase30" / "s1_mom"))
