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
from scipy.ndimage import maximum_filter

C0 = 299792458.0
F = 5.8e9
K0 = 2 * math.pi * F / C0
THETA_BEAM = math.radians(30)


def angular_distance(theta, phi, target_phi):
    dot = np.cos(theta) * math.cos(THETA_BEAM) + np.sin(theta) * math.sin(THETA_BEAM) * np.cos(phi - target_phi)
    return np.arccos(np.clip(dot, -1, 1))


def estimate(radius_m, quadrants, step_m=0.005, angle_step_deg=2, alpha_np_m=None, soft_boundary_mm=0.0,
             feed_floor_m=0.0058):
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
    if alpha_np_m is None:
        radial = np.ones(len(x))
    else:
        rho = np.hypot(x, y)
        radial = np.maximum(rho, feed_floor_m) ** -0.5 * np.exp(-alpha_np_m * rho)
        radial /= np.sqrt(np.mean(radial**2))
    if soft_boundary_mm and quadrants:
        raise ValueError("soft boundary is defined here only for the two-sector case")
    if soft_boundary_mm:
        width = soft_boundary_mm * 1e-3
        right = 0.5 * (1.0 + np.tanh(x / width))
        sector_amplitude = [np.sqrt(right), np.sqrt(1.0 - right)]
    else:
        sector_amplitude = [(labels == i).astype(float) for i in range(sectors)]
    for i, target_phi in enumerate(phis):
        phase = np.sin(THETA_BEAM) * (x * np.cos(target_phi) + y * np.sin(target_phi))
        weights[hands[i]] += radial * sector_amplitude[i] * np.exp(-1j * K0 * phase)

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
        surface = intensity[hand].reshape(tt.shape)
        local_max = surface >= maximum_filter(surface, size=(3, 3), mode=("nearest", "wrap")) * (1 - 1e-12)
        side_candidates = outside & local_max.ravel()
        if np.any(side_candidates):
            side_idx = np.flatnonzero(side_candidates)[int(np.argmax(intensity[hand][side_candidates]))]
            sll = 10 * math.log10(float(intensity[hand][side_idx]) / co)
            side_lobe = {"theta_deg": float(np.degrees(t[side_idx])), "phi_deg": float(np.degrees(p[side_idx])),
                         "hand": hand, "level_db": sll}
        else:
            sll, side_lobe = None, None
        beams.append({"target_theta_deg": 30.0, "target_phi_deg": round(math.degrees(target_phi)), "hand": hand,
                      "peak_theta_deg": float(np.degrees(t[idx])), "peak_phi_deg": float(np.degrees(p[idx])),
                      "directivity_dbi": 10 * math.log10(4 * math.pi * co / power),
                      "sll_db": sll, "worst_sll_lobe": side_lobe, "ar_db": ar,
                      "xpd_db": 10 * math.log10(co / cross) if cross > 0 else None,
                      "targets_met": ar <= 3 and (cross == 0 or co / cross >= 10**1.5) and sll is not None and sll <= -15})
    return {"aperture": "D150" if radius_m == 0.075 else "R150", "radius_mm": radius_m * 1000,
            "partition": f"{sectors}-sector", "sample_count": len(x), "step_mm": step_m * 1000,
            "illumination": "uniform" if alpha_np_m is None else "rho^-1/2 exp(-alpha rho)",
            "alpha_np_m": alpha_np_m, "soft_boundary_mm": soft_boundary_mm,
            "feed_floor_mm": feed_floor_m * 1000 if alpha_np_m is not None else None,
            "angle_step_deg": angle_step_deg, "first_null_exclusion_deg": math.degrees(first_null_rad),
            "first_null_method": "local maxima outside 1.22 lambda / equal-area circular-diameter cone; approximate for sectors",
            "beams": beams}


def main(out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = [estimate(r, q) for r in (0.075, 0.150) for q in (False, True)]
    (out / "aperture.json").write_text(json.dumps({"label": "MODEL_ONLY / IDEAL_APERTURE", "frequency_hz": F,
        "assumptions": "uniform amplitude; perfect circular polarization per sector; independent coherent sectors; cos(theta) field factor; upper hemisphere only; no feed, loss, coupling or physical cell", "cases": rows}, indent=2) + "\n")
    with (out / "aperture.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, lineterminator="\n", fieldnames=["aperture", "partition", "target_theta_deg", "target_phi_deg", "hand", "peak_theta_deg",
            "peak_phi_deg", "directivity_dbi", "sll_db", "sll_theta_deg", "sll_phi_deg", "sll_hand", "ar_db", "xpd_db", "targets_met"])
        writer.writeheader()
        for row in rows:
            for beam in row["beams"]:
                side_lobe = beam.get("worst_sll_lobe") or {}
                writer.writerow({"aperture": row["aperture"], "partition": row["partition"],
                                 **{k: v for k, v in beam.items() if k != "worst_sll_lobe"},
                                 "sll_theta_deg": side_lobe.get("theta_deg"),
                                 "sll_phi_deg": side_lobe.get("phi_deg"), "sll_hand": side_lobe.get("hand")})
    print(json.dumps(rows, indent=2), flush=True)


if __name__ == "__main__":
    import sys
    main(sys.argv[1])
