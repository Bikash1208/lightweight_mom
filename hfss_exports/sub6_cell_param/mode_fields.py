"""Locate an eigenmode: solve OFF-along 90 deg on the BAR50 nominal cell, then export |E|, Ez (volume grid) and |Jsurf|
(metal surface) for chosen modes, and reduce them to region fractions + PNG slices (pure numpy PNG writer, no plotting libs).

python mode_fields.py --project <scratch project copy> --output <dir> --modes 1 2 [--no-solve]
The mode selector for Eigenmode field export is probed at run time (candidates below) and the working one is recorded.
"""
import argparse
import json
import struct
import sys
import zlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cell_param_pyaedt as cp  # noqa: E402


def write_png(path, rgb):
    """Minimal PNG writer: rgb uint8 array (h, w, 3)."""
    h, w, _ = rgb.shape
    raw = b"".join(b"\x00" + rgb[y].tobytes() for y in range(h))
    def chunk(t, d):
        c = struct.pack(">I", len(d)) + t + d
        return c + struct.pack(">I", zlib.crc32(t + d) & 0xFFFFFFFF)
    path.write_bytes(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
                     + chunk(b"IDAT", zlib.compress(raw, 6)) + chunk(b"IEND", b""))


def colormap(a):
    """Simple perceptual-ish map: dark blue -> cyan -> yellow -> red."""
    a = np.clip(a, 0, 1)
    r = np.clip(4 * a - 2, 0, 1); g = np.clip(np.where(a < .5, 4 * a, 4 - 4 * a), 0, 1); b = np.clip(2 - 4 * a, 0, 1)
    return (np.stack([r, g, b], -1) * 255).astype(np.uint8)


def read_fld(path):
    """Read an HFSS .fld grid export: rows 'x y z value' after the header lines."""
    pts = []
    with open(path, errors="replace") as f:
        for line in f:
            t = line.split()
            if len(t) >= 4:
                try:
                    pts.append([float(x) for x in t[:4]])
                except ValueError:
                    pass
    return np.array(pts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--modes", type=int, nargs="+", default=[1, 2])
    ap.add_argument("--no-solve", action="store_true")
    a = ap.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    cp.set_profile("bar50"); cp.ROW_TIMEOUT_S = 30 * 60
    h = cp.connect(a.project)
    rep = {"modes_requested": a.modes}
    try:
        v = cp.nominal(); v["Rs"], v["Xs"] = cp.state_values("OFF", v); v["$ph_x"], v["$ph_y"] = 0, 90
        if not a.no_solve:
            r = cp.solve_row(h, a.output, "FIELDS_OFF_along_90", v, 16)
            rep["modes"] = [(m["mode"], round(m["f_ghz"], 4), round(m["q"], 1)) for m in r["modes"]]
        rep["design_solution_intrinsics"] = list(h.design_solutions.intrinsics)
        try:
            rep["setup_default_intrinsics"] = h.setups[0].default_intrinsics
        except Exception as exc:
            rep["setup_default_intrinsics"] = f"ERR {exc}"
        freqs = {m["mode"]: m["f_ghz"] for m in r["modes"]} if not a.no_solve else {}
        # Grid covering the whole cell: x,y in +-5.8, z from FR4/pad plane to just above the copper
        grid = dict(grid_type="Cartesian", grid_start=[-5.8, -5.8, -4.1], grid_stop=[5.8, 5.8, 0.2], grid_step=[0.2, 0.2, 0.1])
        tried = []
        for mode in a.modes:
            f = freqs.get(mode)
            cands = [{"Freq": f"{f}GHz", "Phase": "0deg"}, {"Freq": f"{f:.6f}GHz", "Phase": "0deg"}, {"Phase": "0deg"},
                     {"Mode": str(mode), "Phase": "0deg"}] if f else [{"Phase": "0deg"}]
            for cand in cands:
                fn = a.output / f"Mag_E_mode{mode}.fld"
                try:
                    ok = h.post.export_field_file_on_grid("Mag_E", solution=f"{cp.SETUP} : LastAdaptive", file_name=str(fn),
                                                          intrinsics=cand, **grid)
                    tried.append({"mode": mode, "intrinsics": cand, "ok": bool(ok) and fn.exists()})
                    if ok and fn.exists():
                        break
                except Exception as exc:
                    tried.append({"mode": mode, "intrinsics": cand, "error": str(exc)[:200]})
        rep["export_attempts"] = tried
        # Reduce whatever exported: energy-like |E|^2 fraction per region and PNG slices
        regions = {}
        for mode in a.modes:
            fn = a.output / f"Mag_E_mode{mode}.fld"
            if not fn.exists():
                continue
            d = read_fld(fn)
            if len(d) == 0:
                continue
            x, y, z, e = (d[:, 0] * 1e3 if abs(d[:, 0]).max() < 1 else d[:, 0], d[:, 1] * 1e3 if abs(d[:, 1]).max() < 1 else d[:, 1],
                          d[:, 2] * 1e3 if abs(d[:, 2]).max() < 1 else d[:, 2], d[:, 3])
            e2 = e ** 2; tot = e2.sum()
            il = v["island_l"]; left = -v["slot_w"] / 2 + v["diode_gap"]
            ix, iy = v["island_via_x"], v["island_via_y"]; dx, dy = v["disc_via_x"], v["disc_via_y"]
            zfr4_top, zfr4_bot = -3.21, -4.01
            sel = {
                "FR4 (z -4.01..-3.21)": (z <= zfr4_top) & (z >= zfr4_bot),
                "RT5880 (z -3.175..0)": (z > -3.175) & (z <= 0),
                "above copper (z>0.035)": z > 0.035,
                "island via column r<0.6": np.hypot(x - ix, y - iy) < 0.6,
                "disc via column r<0.6": np.hypot(x - dx, y - dy) < 0.6,
                "island stub zone (x+-0.3, y from via to via+2)": (abs(x - ix) < 0.3) & (y >= iy) & (y <= iy + v["stub_l"]) & (z < -3.2),
                "diode slot gap (island left edge .. slot edge)": (x >= -v["slot_w"] / 2 - 0.2) & (x <= left + 0.2) & (abs(y - v["island_offset"]) < il / 2),
                "island footprint": (x >= left) & (x <= left + v["island_w"]) & (abs(y - v["island_offset"]) < il / 2),
            }
            regions[f"mode{mode}"] = {k: round(float(e2[m].sum() / tot), 4) for k, m in sel.items()}
            # peak location + PNG slice of |E| at the z-plane nearest the copper and mid-FR4
            k = int(np.argmax(e)); regions[f"mode{mode}"]["peak_xyz_mm"] = [round(float(x[k]), 2), round(float(y[k]), 2), round(float(z[k]), 2)]
            for tag, zc in (("copper", 0.0), ("fr4mid", -3.6)):
                zs = np.unique(np.round(z, 3)); zsel = zs[np.argmin(abs(zs - zc))]
                m = np.isclose(np.round(z, 3), zsel)
                xs = np.unique(np.round(x[m], 3)); ys = np.unique(np.round(y[m], 3))
                img = np.zeros((len(ys), len(xs)))
                xi = np.searchsorted(xs, np.round(x[m], 3)); yi = np.searchsorted(ys, np.round(y[m], 3))
                img[yi, xi] = e[m]
                img = np.log10(img + 1e-3 * img.max() + 1e-30); img = (img - img.min()) / max(np.ptp(img), 1e-30)
                big = np.kron(colormap(img[::-1]), np.ones((3, 3, 1), dtype=np.uint8))
                write_png(a.output / f"mode{mode}_logE_{tag}_z{zsel:g}.png", big)
        rep["region_energy_fraction"] = regions
        rep["status"] = "DONE"
    except Exception as exc:
        rep["status"] = "ERROR"; rep["error"] = repr(exc)
        raise
    finally:
        (a.output / "mode_fields.json").write_text(json.dumps(rep, indent=2, default=str) + "\n", encoding="utf-8")
        try:
            h.save_project()
        except Exception:
            pass
        try:
            h.release_desktop()
        except Exception:
            pass
    print(json.dumps(rep, indent=2, default=str)[:3000])


if __name__ == "__main__":
    main()
