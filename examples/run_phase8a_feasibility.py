#!/usr/bin/env python3
"""Phase 8A CLI: run the scientific feasibility gate against the Phase 7 production library.

    python examples/run_phase8a_feasibility.py [--library-dir DIR] [--reference-field FILE]
                                                [--frequency-ghz F] [--polarization x|y]
                                                [--tolerance-deg T] [--output-dir DIR]

With no arguments, this runs Mode A (library-only feasibility) against
``results/unit_cell_library`` and writes ``results/phase8a_feasibility/``. Exit code:
`feasibility_schema.EXIT_CODES[status]`.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.feasibility_gate import run_and_write_phase8a
from src.feasibility_schema import DEFAULT_LIBRARY_DIR, FeasibilityConfig


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library-dir", default=str(DEFAULT_LIBRARY_DIR))
    parser.add_argument("--reference-field", default=None, help="Mode B: reference-field CSV (see reference_field.py)")
    parser.add_argument("--frequency-ghz", type=float, default=None)
    parser.add_argument("--polarization", default="x", choices=["x", "y"])
    parser.add_argument("--tolerance-deg", type=float, default=None)
    parser.add_argument("--expected-geometries", type=int, default=324)
    parser.add_argument("--expected-frequency-points", type=int, default=201)
    parser.add_argument("--output-dir", default="results/phase8a_feasibility")
    args = parser.parse_args(argv)

    kwargs = {
        "library_dir": args.library_dir, "polarization": args.polarization,
        "expected_n_geometries": args.expected_geometries,
        "expected_n_frequency_points": args.expected_frequency_points,
    }
    if args.frequency_ghz is not None:
        kwargs["frequency_hz"] = args.frequency_ghz * 1e9
    if args.tolerance_deg is not None:
        kwargs["phase_match_tolerance_deg"] = args.tolerance_deg
    config = FeasibilityConfig(**kwargs)

    result, output_dir = run_and_write_phase8a(
        config, args.output_dir, reference_field_path=args.reference_field,
    )

    def pf(label, value):
        print(f"{label:<32} {value}")

    print("=" * 50)
    print("PHASE 8A -- SCIENTIFIC FEASIBILITY GATE")
    print("=" * 50)
    print()
    pf("Library integrity:", result.checks.get("library_integrity", "N/A"))
    pf("Phase source validation:", result.checks.get("phase_source_validation", "N/A"))
    pf("Frequency compatibility:", result.checks.get("frequency_compatibility", "N/A"))
    pf("Polarization compatibility:", result.checks.get("polarization_compatibility", "N/A"))
    pf("Reference-field validity:", result.checks.get("reference_field_validity", "N/A"))
    pf("Grid mapping validity:", result.checks.get("grid_mapping_validity", "N/A"))
    print()
    cov = result.library_phase_coverage or {}
    pf("Library phase coverage:", f"{cov.get('span_deg', float('nan')):.1f} deg" if cov else "N/A")
    pf("Largest phase gap:", f"{cov.get('largest_gap_deg', float('nan')):.1f} deg" if cov else "N/A")
    rp = result.required_phase or {}
    pf("Required phase coverage:", f"{rp.get('coverage_deg', float('nan')):.1f} deg" if rp.get("coverage_deg") is not None else "N/A")
    pf("Phase matching:", result.checks.get("phase_matching", "N/A"))
    m = result.matching or {}
    pf("Mean phase error:", f"{m.get('mean_error_deg'):.1f} deg" if m.get("mean_error_deg") is not None else "N/A")
    pf("Maximum phase error:", f"{m.get('max_error_deg'):.1f} deg" if m.get("max_error_deg") is not None else "N/A")
    pf("Grating check:", result.checks.get("grating_check", "N/A"))
    pf("HFSS validation:", "NOT PERFORMED")
    print()
    print("=" * 50)
    print("FINAL STATUS:")
    print(result.status)
    print("=" * 50)
    print()
    print(f"Output: {output_dir}")

    return result.exit_code


if __name__ == "__main__":
    raise SystemExit(main())
