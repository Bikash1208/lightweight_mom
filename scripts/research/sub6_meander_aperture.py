"""Ideal upper-hemisphere aperture/array-factor estimate, MODEL_ONLY.

Two aperture sizes: D=150 mm and R=150 mm. Aperture samples are uniform 5 mm
square pixels within the circle. Each half/quadrant has unit-amplitude,
perfect circularly polarised current with a phase ramp toward its assigned
30-degree beam. Field element factor cos(theta); no feed, coupling, loss,
edge taper, or real CP element is included. These are feasibility bounds.
"""

import csv
import json
import math
from pathlib import Path

import numpy as np

C0 = 299792458.0
F = 5.8e9
K0 = 2 * math.pi * F / C0
THETA_BEAM = math.radians(30)


def angular_distance(theta, phi, target_phi):
    dot = np.cos(theta) * math.cos(THETA_BEAM) + np.sin(theta) * math.sin(THETA_BEAM) * np.cos(phi - target_phi)
    return np.arccos(np.clip(dot, -1, 1))


def estimate(radius_m, quadrants, step_m=0.005, angle_step_deg=2):
    grid = np.arange(-radius_m + step_m / 2, radius_m, step_m)
    xx, yy = np.meshgrid(grid, grid)
    keep = xx**2 + yy**2 <= radius_m**2
    x, y = xx[keep], yy[keep]
    pos = np.column_stack([x, y])
    sectors = 4 if quadrants else 2
    labels = np.mod(np.floor((np.arctan2(y, x) + math.pi / sectors) * sectors / (2 * math.pi)).astype(int), sectors)
    phis = np.arange(sectors) * 2 * math.pi / sectors
    hands = ["RHCP" if i % 2 == 0 else "LHCP" for i in range(sectors)]
    weights = {"RHCP": np.zeros(len(x), dtype=complex), "LHCP": np.zeros(len(x), dtype=complex)}
    for i, target_phi in enumerate(phis):
        sel = labels == i
        phase = np.sin(THETA_BEAM) * (x[sel] * np.cos(target_phi) + y[sel] * np.sin(target_phi))
        weights[hands[i]][sel] = np.exp(-1j * K0 * phase)

    theta_deg = np.arange(0, 90 + angle_step_deg, angle_step_deg)
    phi_deg = np.arange(0, 360, angle_step_deg)
    tt, pp = np.meshgrid(np.radians(theta_deg), np.radians(phi_deg), indexing="ij")
    t, p = tt.ravel(), pp.ravel()
    direction = np.column_stack([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p)])
    fields = {hand: np.empty(len(t), complex) for hand in weights}
    for start in range(0, len(t), 256):
        stop = min(start + 256, len(t))
        phase = np.exp(1j * K0 * (direction[start:stop] @ pos.T))
        for hand in weights:
            fields[hand][start:stop] = (phase @ weights[hand]) * np.cos(t[start:stop])
    intensity = {hand: abs(fields[hand]) ** 2 for hand in weights}
    total = intensity["RHCP"] + intensity["LHCP"]
    d_omega = math.radians(angle_step_deg) ** 2 * np.sin(t)
    power = float(np.sum(total * d_omega))
    first_null_rad = math.asin(min(1.0, 1.22 * (C0 / F) / (2 * radius_m / math.sqrt(sectors))))
    beams = []
    for i, target_phi in enumerate(phis):
        hand = hands[i]
        other = "LHCP" if hand == "RHCP" else "RHCP"
        distance = angular_distance(t, p, target_phi)
        local = distance <= math.radians(20)
        idx = np.flatnonzero(local)[int(np.argmax(intensity[hand][local]))]
        co, cross = float(intensity[hand][idx]), float(intensity[other][idx])
        ratio = math.sqrt(co / cross) if cross > 0 else math.inf
        ar = 20 * math.log10((ratio + 1) / (ratio - 1)) if ratio > 1 and math.isfinite(ratio) else (0.0 if math.isinf(ratio) else math.inf)
        outside = np.ones(len(t), dtype=bool)
        for j, phi_j in enumerate(phis):
            if hands[j] == hand:
                outside &= angular_distance(t, p, phi_j) > first_null_rad
        sll = 10 * math.log10(float(np.max(intensity[hand][outside])) / co) if np.any(outside) else None
        beams.append({"target_theta_deg": 30.0, "target_phi_deg": round(math.degrees(target_phi)), "hand": hand,
                      "peak_theta_deg": float(np.degrees(t[idx])), "peak_phi_deg": float(np.degrees(p[idx])),
                      "directivity_dbi": 10 * math.log10(4 * math.pi * co / power),
                      "sll_db": sll, "ar_db": ar,
                      "xpd_db": 10 * math.log10(co / cross) if cross > 0 else None,
                      "targets_met": ar <= 3 and (cross == 0 or co / cross >= 10**1.5) and sll is not None and sll <= -15})
    return {"aperture": "D150" if radius_m == 0.075 else "R150", "radius_mm": radius_m * 1000,
            "partition": f"{sectors}-sector", "sample_count": len(x), "step_mm": step_m * 1000,
            "angle_step_deg": angle_step_deg, "first_null_exclusion_deg": math.degrees(first_null_rad),
            "first_null_method": "1.22 lambda / equal-area circular diameter; approximate for sectors",
            "beams": beams}


def main(out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = [estimate(r, q) for r in (0.075, 0.150) for q in (False, True)]
    (out / "aperture.json").write_text(json.dumps({"label": "MODEL_ONLY / IDEAL_APERTURE", "frequency_hz": F,
        "assumptions": "uniform amplitude; perfect circular polarization per sector; independent coherent sectors; cos(theta) field factor; upper hemisphere only; no feed, loss, coupling or physical cell", "cases": rows}, indent=2) + "\n")
    with (out / "aperture.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, lineterminator="\n", fieldnames=["aperture", "partition", "target_theta_deg", "target_phi_deg", "hand", "peak_theta_deg",
            "peak_phi_deg", "directivity_dbi", "sll_db", "ar_db", "xpd_db", "targets_met"])
        writer.writeheader()
        for row in rows:
            for beam in row["beams"]:
                writer.writerow({"aperture": row["aperture"], "partition": row["partition"], **beam})
    print(json.dumps(rows, indent=2), flush=True)


if __name__ == "__main__":
    import sys
    main(sys.argv[1])
