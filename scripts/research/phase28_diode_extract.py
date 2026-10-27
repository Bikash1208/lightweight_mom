"""Phase 28: series impedance of the MACOM MA4AGP907 (= MADP-000907-14020, AlGaAs flip-chip PIN) from the vendor-measured 2-port S-parameters (HP8510C, 2004; series-mounted fixture).

    python scripts/research/phase28_diode_extract.py [--refs C:/hmsa_sync/refs/diodes] [--out results/phase28/diodes/ma4agp907_series_z.json]

Series element between two identical fixture lines: S11 = S22 = Z/(Z + 2 Z0), S21 = 2 Z0/(Z + 2 Z0), hence Z = Z0 (S11 + S22)/S21, which cancels symmetric fixture lines (advisor extraction, reproduced here).
ON files give R = Re Z and X = Im Z (R in ohm, X in ohm); OFF files give C = -1/(omega Im Z). The .s2p files are vendor data and are NEVER committed; the extracted values (Z(f), 4 - 20 GHz) and the provenance (file names, sha256,
bias read from the file header, file date) are. Labels: MEASURED_VENDOR_SPARAM.
"""

import hashlib
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

Z0 = 50.0
FILES = {
    "ON_5mA": "MA4AGP907_5mA_Forward_Bias.s2p",
    "ON_10mA": "MA4AGP907_10mA_Forward_Bias.s2p",
    "ON_15mA": "MA4AGP907_15mA_Forward_Bias.s2p",
    "OFF_0V": "MA4AGP907_0V_Reverse_Bias.s2p",
    "OFF_m5V": "MA4AGP907_-5V_Reverse_Bias.s2p",
}
F_BAND_GHZ = (4.0, 20.0)


def read_s2p(path):
    """`(freq_ghz, s)` with `s` of shape (n, 2, 2) from a Touchstone `# GHZ S RI R 50` file (column order S11 S21 S12 S22), plus the header text."""
    header, rows = [], []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("!") or line.startswith("#"):
            header.append(line)
        elif line.strip():
            rows.append([float(v) for v in line.split()])
    a = np.array(rows)
    s = np.empty((len(a), 2, 2), dtype=complex)
    s[:, 0, 0] = a[:, 1] + 1j * a[:, 2]
    s[:, 1, 0] = a[:, 3] + 1j * a[:, 4]
    s[:, 0, 1] = a[:, 5] + 1j * a[:, 6]
    s[:, 1, 1] = a[:, 7] + 1j * a[:, 8]
    return a[:, 0], s, header


def series_z(s):
    return Z0 * (s[:, 0, 0] + s[:, 1, 1]) / s[:, 1, 0]


def bias_from_header(header):
    m = re.search(r"CH2 I=\s*([-+0-9.Ee]+)\s*V=\s*([-+0-9.Ee]+)", " ".join(header))
    return {"i_a": float(m.group(1)), "v_v": float(m.group(2))} if m else None


def extract(refs):
    out = {"labels": ["MEASURED_VENDOR_SPARAM"], "source": "MACOM MA4AGP907 measured S-parameters (HP8510C.07.14, 20 Apr 2004), series-mounted; Z = Z0 (S11 + S22)/S21", "z0_ohm": Z0, "states": {}}
    for key, name in FILES.items():
        p = Path(refs) / name
        f, s, header = read_s2p(p)
        z = series_z(s)
        band = (f >= F_BAND_GHZ[0] - 1e-9) & (f <= F_BAND_GHZ[1] + 1e-9)
        w = 2 * math.pi * f * 1e9
        entry = {
            "file": name,
            "sha256": hashlib.sha256(Path(p).read_bytes()).hexdigest(),
            "file_date": header[0].lstrip("! ").strip(),
            "bias_header": bias_from_header(header),
            "f_ghz": [float(v) for v in f[band]],
            "z_re_ohm": [float(v) for v in z.real[band]],
            "z_im_ohm": [float(v) for v in z.imag[band]],
        }
        if key.startswith("ON"):
            entry["R_mean_4_20_ghz"] = float(np.mean(z.real[band]))
            entry["X_min_max_4_20_ghz"] = [float(np.min(z.imag[band])), float(np.max(z.imag[band]))]
        else:
            c = -1.0 / (w * z.imag)
            entry["c_ff"] = [float(v * 1e15) for v in c[band]]
            entry["c_ff_mean_4_20_ghz"] = float(np.mean(c[band]) * 1e15)
        i10 = int(np.argmin(np.abs(f - 10.0)))
        entry["at_10_ghz"] = {"z_ohm": [float(z[i10].real), float(z[i10].imag)], **({"c_ff": float(-1.0 / (w[i10] * z[i10].imag) * 1e15)} if key.startswith("OFF") else {})}
        out["states"][key] = entry
    return out


def main(argv):
    refs = argv[argv.index("--refs") + 1] if "--refs" in argv else "C:/hmsa_sync/refs/diodes"
    dst = Path(argv[argv.index("--out") + 1] if "--out" in argv else Path(__file__).resolve().parents[2] / "results/phase28/diodes/ma4agp907_series_z.json")
    doc = extract(refs)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(doc, indent=1), encoding="utf-8")
    for k, v in doc["states"].items():
        print(k, v["at_10_ghz"], v.get("R_mean_4_20_ghz"), v.get("X_min_max_4_20_ghz"), v.get("c_ff_mean_4_20_ghz"))


if __name__ == "__main__":
    main(sys.argv[1:])
