"""Phase 12A required outputs: deterministic CSV/JSON (no timestamps, no host names), one writer
per file, matching the project's per-phase convention (``phase11.specification``, etc.).
"""

import csv
import json
from dataclasses import asdict


def write_sweep_csv(points, path):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["length_width_mm", "frequency_hz", "s11_xx_mag", "s11_xx_phase_deg", "s11_yx_mag"])
        for p in sorted(points, key=lambda p: (p.frequency_hz, p.length_width_mm)):
            writer.writerow(
                [
                    repr(p.length_width_mm),
                    repr(p.frequency_hz),
                    repr(p.s11_xx_mag),
                    repr(p.s11_xx_phase_deg),
                    repr(p.s11_yx_mag),
                ]
            )
    return path


def _write_json(doc, path):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(json.dumps(doc, indent=1, allow_nan=False) + "\n")
    return path


def write_sweep_json(points, coverage_by_frequency, path):
    doc = {
        "points": [asdict(p) for p in points],
        "coverage_by_frequency": {
            str(f): {
                "coverage_deg": c.coverage_deg,
                "largest_gap_deg": c.largest_gap_deg,
                "phase_min_deg": c.phase_min_deg,
                "phase_max_deg": c.phase_max_deg,
                "n_points": c.n_points,
            }
            for f, c in coverage_by_frequency.items()
        },
    }
    return _write_json(doc, path)


def write_ungrounded_comparison_json(ungrounded_result, path):
    return _write_json(dict(ungrounded_result), path)


def write_analytic_checks_json(slab_points, thin_limit, path):
    doc = {
        "bare_grounded_slab": [
            {
                "frequency_hz": p.frequency_hz,
                "eps_r": p.eps_r,
                "tan_delta": p.tan_delta,
                "thickness_m": p.thickness_m,
                "gamma_mine": {"re": p.gamma_mine.real, "im": p.gamma_mine.imag},
                "gamma_textbook": {"re": p.gamma_textbook.real, "im": p.gamma_textbook.imag},
                "max_error": p.max_error,
            }
            for p in slab_points
        ],
        "vanishing_thickness_limit": {"re": thin_limit.real, "im": thin_limit.imag},
    }
    return _write_json(doc, path)


def write_surface_wave_json(points, path):
    return _write_json([asdict(p) for p in points], path)


def write_convergence_json(points, floquet_error, mesh_error, path):
    doc = {
        "points": [asdict(p) for p in points],
        "refinement_study": {"floquet_error": floquet_error, "mesh_error": mesh_error},
    }
    return _write_json(doc, path)


def write_gate_status_json(decision, path):
    return _write_json(asdict(decision), path)
