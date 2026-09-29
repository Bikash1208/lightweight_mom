"""Phase 8A orchestrator: runs every gate in `feasibility_schema.STATES`'s order of priority,
assembles `feasibility.json` (Sec. 41) and `FEASIBILITY_REPORT.md` (Sec. 42), and writes every
output file Sec. 40 lists. Nothing here re-solves, interpolates, extrapolates, or invents a
result; every gate module it calls (`feasibility_phase_source`, `feasibility_library_integrity`,
`feasibility_coverage`, `feasibility_matching`, `feasibility_grating`, `feasibility_spatial`,
`reference_field`) already says so in its own docstring.

Priority order (Sec. 36 requires *all* checks to pass for DESIGN_FEASIBLE; this module checks
them in the order below and stops at, and reports, the first one that fails):

    1. library integrity           -> LIBRARY_INVALID
    2. frequency support            -> FREQUENCY_UNSUPPORTED
    3. polarization support         -> POLARIZATION_UNSUPPORTED
    4. reference field (Mode B only) -> REFERENCE_INVALID
    5. phase-source validation      -> PHASE_SOURCE_UNVALIDATED   (always true for this project today)
    6. spatial mapping (Mode B only) -> GRID_MAPPING_UNVALIDATED
    7. phase coverage sufficiency (Mode B only, needs a required-phase array) -> PHASE_COVERAGE_INSUFFICIENT
    8. phase matching (Mode B only) -> PHASE_MATCHING_INSUFFICIENT
    9. grating risk                 -> GRATING_RISK
   10. otherwise                    -> DESIGN_FEASIBLE / DESIGN_FEASIBLE_WITH_WARNINGS

Since step 5 always fails today (`feasibility_phase_source.py`'s module docstring), this is the
status every real run of this gate against the current Phase 7 library reports -- the honest,
expected scientific result (brief Sec. 48), not a software defect.
"""

import json
from dataclasses import dataclass
from pathlib import Path

from . import __version__ as SOLVER_VERSION
from .feasibility_coverage import (
    describe_coverage_vs_resonance,
    library_table_at_frequency,
    phase_coverage_vs_frequency,
    phase_density,
    resonance_summary,
    usable_geometry_count,
)
from .feasibility_grating import grating_check
from .feasibility_library_integrity import check_library_integrity, load_library_bundle
from .feasibility_matching import is_phase_coverage_sufficient, match_required_phase
from .feasibility_phase_source import validate_holographic_phase_source
from .feasibility_plots import (
    plot_circular_phase_coverage,
    plot_coverage_vs_frequency,
    plot_matching_error,
    plot_phase_vs_lw,
)
from .feasibility_required_phase import compute_required_phase, required_phase_statistics
from .feasibility_schema import EXIT_CODES, FeasibilityConfig, PROVIDED_POLARIZATION
from .feasibility_spatial import check_mapping_method_defined
from .phase_extraction import phase_coverage
from .reference_field import validate_beam_definition, validate_reference_field


@dataclass(frozen=True)
class FeasibilityResult:
    status: str
    exit_code: int
    checks: dict
    warnings: tuple
    errors: tuple
    library: dict
    frequency: dict
    phase_source: dict
    library_phase_coverage: dict
    frequency_dependence: dict
    resonance: dict
    grating: dict
    reference_field: dict | None
    required_phase: dict | None
    matching: dict | None
    hfss_validation: dict
    table: list                  # phase-vs-geometry table @ the requested frequency, for reports/plots
    freq_coverage: list          # coverage vs frequency, for reports/plots
    match_results: list | None


def _fail(status, errors, warnings=(), **partial):
    return FeasibilityResult(
        status=status, exit_code=EXIT_CODES[status], checks=partial.pop("checks", {}),
        warnings=tuple(warnings), errors=tuple(errors),
        library=partial.pop("library", {}), frequency=partial.pop("frequency", {}),
        phase_source=partial.pop("phase_source", {}),
        library_phase_coverage=partial.pop("library_phase_coverage", {}),
        frequency_dependence=partial.pop("frequency_dependence", {}),
        resonance=partial.pop("resonance", {}), grating=partial.pop("grating", {}),
        reference_field=partial.pop("reference_field", None),
        required_phase=partial.pop("required_phase", None), matching=partial.pop("matching", None),
        hfss_validation={"performed": False},
        table=partial.pop("table", []), freq_coverage=partial.pop("freq_coverage", []),
        match_results=partial.pop("match_results", None),
    )


def run_feasibility_gate(
    config: FeasibilityConfig, *, reference_field_path=None, beams=None, required_phase_deg=None,
):
    """Runs Mode A (library-only) always; Mode B (design-specific) additionally whenever
    ``reference_field_path``/``beams`` are supplied. ``required_phase_deg``
    (and, for the map plot, ``required_phase_xy_mm``) let a caller who has *already* computed a
    required-phase distribution elsewhere (or a test exercising synthetic data, brief Sec. 44
    Tests 10-11) drive the matching/sufficiency checks without this gate inventing one -- see
    `feasibility_required_phase.py`'s module docstring for why this gate itself never computes
    one from a reference field."""
    checks = {}

    # ---------------------------------------------------------------- 1. library integrity
    bundle = load_library_bundle(config.library_dir)
    integrity = check_library_integrity(bundle, config.expected_n_geometries, config.expected_n_frequency_points)
    checks["library_integrity"] = "PASS" if integrity.ok else "FAIL"
    library_info = {
        "id": bundle.manifest.get("configuration_hash", "")[:16],
        "geometry_count": integrity.n_cells, "successful_cells": integrity.n_ok_cells,
        "failed_cells": integrity.n_failed_cells, "frequency_points": config.expected_n_frequency_points,
        "period_x_mm": bundle.manifest["period_x_mm"], "period_y_mm": bundle.manifest["period_y_mm"],
        "step_files": integrity.n_step_files,
    }
    if not integrity.ok:
        return _fail("LIBRARY_INVALID", integrity.errors, checks=checks, library=library_info)

    # ------------------------------------------------------------------- 2. frequency support
    f_design_ghz = bundle.manifest["design_frequency_GHz"]
    f_request_ghz = (config.frequency_hz / 1e9) if config.frequency_hz is not None else f_design_ghz
    f_min_ghz, f_max_ghz = bundle.manifest["frequency_min_GHz"], bundle.manifest["frequency_max_GHz"]
    ok_cell = next(c for c in bundle.library_doc["cells"] if c.get("solver_status") == "OK")
    on_grid = any(abs(f - f_request_ghz) <= 1e-6 * max(abs(f_request_ghz), 1.0) for f in ok_cell["frequency_GHz"])
    freq_supported = (f_min_ghz <= f_request_ghz <= f_max_ghz) and on_grid
    checks["frequency_compatibility"] = "PASS" if freq_supported else "FAIL"
    frequency_info = {
        "requested_GHz": f_request_ghz, "supported": bool(freq_supported), "exact_grid_point": bool(on_grid),
        "library_min_GHz": f_min_ghz, "library_max_GHz": f_max_ghz,
    }
    if not freq_supported:
        return _fail(
            "FREQUENCY_UNSUPPORTED",
            [f"frequency {f_request_ghz} GHz is outside [{f_min_ghz}, {f_max_ghz}] GHz or not an "
             "exact grid point; no interpolation is performed"],
            checks=checks, library=library_info, frequency=frequency_info,
        )

    # --------------------------------------------------------------- 3. polarization support
    pol_supported = config.polarization == PROVIDED_POLARIZATION
    checks["polarization_compatibility"] = "PASS" if pol_supported else "FAIL"
    if not pol_supported:
        return _fail(
            "POLARIZATION_UNSUPPORTED",
            [f"requested polarization {config.polarization!r} is not provided by this library "
             f"(only {PROVIDED_POLARIZATION!r}, the co-polarized S11xx/S21xx component pinned by "
             "Phase 7's phase-vs-geometry pipeline); no silent substitution is performed"],
            checks=checks, library=library_info, frequency=frequency_info,
        )

    # ----------------------------------------------------------- 4. reference field (Mode B)
    reference_field_info = None
    ref_validation = None
    if reference_field_path is not None:
        ref_validation = validate_reference_field(reference_field_path, expected_frequency_ghz=f_request_ghz)
        checks["reference_field_validity"] = "PASS" if ref_validation.ok else "FAIL"
        reference_field_info = {
            "ok": ref_validation.ok, "n_samples": ref_validation.n_samples,
            "errors": list(ref_validation.errors), "dx_mm": ref_validation.dx_mm, "dy_mm": ref_validation.dy_mm,
        }
        if not ref_validation.ok:
            return _fail(
                "REFERENCE_INVALID", list(ref_validation.errors), checks=checks,
                library=library_info, frequency=frequency_info, reference_field=reference_field_info,
            )
    else:
        checks["reference_field_validity"] = "N/A (Mode A)"

    if beams is not None:
        for i, beam in enumerate(beams):
            bv = validate_beam_definition(beam)
            if not bv.ok:
                return _fail(
                    "REFERENCE_INVALID", [f"beam[{i}]: {e}" for e in bv.errors], checks=checks,
                    library=library_info, frequency=frequency_info, reference_field=reference_field_info,
                )

    # ---------------------------------------------------------- 5. phase-source validation
    phase_source = validate_holographic_phase_source()
    phase_source_info = {
        "quantity": phase_source.quantity_provided,
        "quantity_required": phase_source.quantity_required,
        "validated_for_holographic_use": phase_source.status == "PHYSICALLY_VALIDATED",
        "status": phase_source.status, "note": phase_source.note,
    }
    checks["phase_source_validation"] = "PASS" if phase_source.status == "PHYSICALLY_VALIDATED" else "FAIL"

    # --------------------------------------------------------- library-only analysis (Mode A)
    # Computed regardless of the phase-source verdict above -- it is useful, honestly-labelled
    # data even while the phase-source question remains open (same posture as Phase 7's report).
    table = library_table_at_frequency(bundle.library_doc, f_request_ghz)
    coverage = phase_coverage(table)
    density = phase_density(table)
    n_usable = usable_geometry_count(table)
    library_coverage_info = {
        "span_deg": coverage.coverage_deg, "largest_gap_deg": coverage.largest_gap_deg,
        "naive_span_deg": coverage.naive_span_deg, "n_points": coverage.n_points,
        "n_usable_geometries": n_usable, "n_distinct_phase_states": density.n_distinct_phase_states,
        "smallest_gap_deg": density.smallest_gap_deg, "median_gap_deg": density.median_gap_deg,
        "phase_density_per_deg": density.phase_density_per_deg,
    }

    freq_coverage = phase_coverage_vs_frequency(bundle.library_doc)
    n_resonant, n_non_resonant, resonant_freqs_hz = resonance_summary(bundle.library_doc)
    coverage_vs_resonance = describe_coverage_vs_resonance(freq_coverage, n_resonant, resonant_freqs_hz)
    resonance_info = {
        "n_resonant": n_resonant, "n_non_resonant": n_non_resonant,
        "resonance_frequencies_GHz": [f / 1e9 for f in resonant_freqs_hz],
        "coverage_vs_resonance": coverage_vs_resonance,
    }
    checks["resonance_analysis"] = "COMPUTED"

    grating = grating_check(
        bundle.manifest["period_x_mm"] * 1e-3, bundle.manifest["period_y_mm"] * 1e-3, f_request_ghz * 1e9,
    )
    checks["grating_check"] = "WARNING" if grating.is_grating_risk else "PASS"
    grating_info = {
        "period_m": grating.period_x_m, "frequency_hz": grating.frequency_hz,
        "wavelength_m": grating.wavelength_m, "period_over_wavelength": grating.period_over_wavelength,
        "grating_onset_hz": grating.grating_onset_hz, "is_grating_risk": grating.is_grating_risk,
    }

    if phase_source.status != "PHYSICALLY_VALIDATED":
        return _fail(
            "PHASE_SOURCE_UNVALIDATED", [
                "The current unit-cell scattering quantity has not yet been validated as the "
                "physical phase-control quantity required by the holographic synthesis "
                "formulation.",
            ],
            warnings=(
                [f"even under the tentative {phase_source.quantity_provided} assumption, library "
                 f"phase coverage at {f_request_ghz:g} GHz is only {coverage.coverage_deg:.2f} deg "
                 "of 360 deg (largest gap "
                 f"{coverage.largest_gap_deg:.2f} deg)"]
            ),
            checks=checks, library=library_info, frequency=frequency_info,
            phase_source=phase_source_info, library_phase_coverage=library_coverage_info,
            frequency_dependence={"points": freq_coverage}, resonance=resonance_info,
            grating=grating_info, reference_field=reference_field_info, table=table, freq_coverage=freq_coverage,
        )

    # ------------------------------------------------------- 6. spatial mapping (Mode B only)
    mapping = None
    if reference_field_path is not None:
        mapping = check_mapping_method_defined()
        checks["grid_mapping_validity"] = "PASS" if mapping.defined else "FAIL"
        if not mapping.defined:
            return _fail(
                "GRID_MAPPING_UNVALIDATED", [mapping.note], checks=checks,
                library=library_info, frequency=frequency_info, phase_source=phase_source_info,
                library_phase_coverage=library_coverage_info, frequency_dependence={"points": freq_coverage},
                resonance=resonance_info, grating=grating_info, reference_field=reference_field_info,
                table=table, freq_coverage=freq_coverage,
            )
    else:
        checks["grid_mapping_validity"] = "N/A (Mode A)"

    # ------------------------------------------------- 7-8. required phase / matching (Mode B)
    required_phase_info, matching_info, match_results = None, None, None
    if reference_field_path is not None or beams is not None:
        rp = compute_required_phase(ref_validation, beams, config.polarization)
        required_phase_info = {"available": rp.available, "status": rp.status, "note": rp.note}
        if not rp.available:
            # Same underlying cause as the phase-source gate (no holographic formulation exists
            # in this project, feasibility_required_phase.py's docstring) -- reported under the
            # same state rather than inventing a new one outside Sec. 35's fixed list.
            checks["required_phase_available"] = "FAIL"
            return _fail(
                "PHASE_SOURCE_UNVALIDATED",
                [rp.note], checks=checks, library=library_info, frequency=frequency_info,
                phase_source=phase_source_info, library_phase_coverage=library_coverage_info,
                frequency_dependence={"points": freq_coverage}, resonance=resonance_info,
                grating=grating_info, reference_field=reference_field_info,
                required_phase=required_phase_info, table=table, freq_coverage=freq_coverage,
            )
    elif required_phase_deg is not None:
        # A caller (or a test) supplied an already-computed required-phase array directly,
        # bypassing the (currently nonexistent) generator -- see run_feasibility_gate's docstring.
        stats = required_phase_statistics(required_phase_deg)
        required_phase_info = {
            "available": True, "coverage_deg": stats.circular_coverage_deg,
            "largest_gap_deg": stats.largest_gap_deg, "n_samples": stats.n_samples,
        }
        sufficiency = is_phase_coverage_sufficient(required_phase_deg, table, config.phase_match_tolerance_deg)
        match_results = match_required_phase(required_phase_deg, table, config.phase_match_tolerance_deg)
        matching_info = {
            "matched": sufficiency.n_matched, "unmatched": sufficiency.n_unmatched,
            "mean_error_deg": sufficiency.mean_error_deg, "max_error_deg": sufficiency.max_error_deg,
        }
        checks["phase_matching"] = "PASS" if sufficiency.sufficient else "FAIL"
        if not sufficiency.sufficient:
            return _fail(
                "PHASE_COVERAGE_INSUFFICIENT" if sufficiency.n_matched == 0 else "PHASE_MATCHING_INSUFFICIENT",
                [f"{sufficiency.n_unmatched} of {sufficiency.n_required} required phase sample(s) "
                 f"have no library candidate within {config.phase_match_tolerance_deg:g} deg"],
                checks=checks, library=library_info, frequency=frequency_info,
                phase_source=phase_source_info, library_phase_coverage=library_coverage_info,
                frequency_dependence={"points": freq_coverage}, resonance=resonance_info,
                grating=grating_info, reference_field=reference_field_info,
                required_phase=required_phase_info, matching=matching_info,
                table=table, freq_coverage=freq_coverage, match_results=match_results,
            )

    # ------------------------------------------------------------------------- 9. grating risk
    if grating.is_grating_risk:
        warnings = ["grating-order onset is at or below the requested frequency"]
        status = "DESIGN_FEASIBLE_WITH_WARNINGS"
    else:
        warnings = []
        status = "DESIGN_FEASIBLE"

    return FeasibilityResult(
        status=status, exit_code=EXIT_CODES[status], checks=checks, warnings=tuple(warnings), errors=(),
        library=library_info, frequency=frequency_info, phase_source=phase_source_info,
        library_phase_coverage=library_coverage_info, frequency_dependence={"points": freq_coverage},
        resonance=resonance_info, grating=grating_info, reference_field=reference_field_info,
        required_phase=required_phase_info, matching=matching_info, hfss_validation={"performed": False},
        table=table, freq_coverage=freq_coverage, match_results=match_results,
    )


# --------------------------------------------------------------------------- machine-readable output

def build_feasibility_json(result: FeasibilityResult, config: FeasibilityConfig):
    """Sec. 41's schema, populated entirely from ``result`` -- no example values hard-coded."""
    return {
        "status": result.status,
        "library": result.library,
        "frequency": result.frequency,
        "phase_source": result.phase_source,
        "phase_coverage": {
            "span_deg": result.library_phase_coverage.get("span_deg"),
            "largest_gap_deg": result.library_phase_coverage.get("largest_gap_deg"),
        },
        "required_phase": (
            {"available": result.required_phase.get("available"),
             "coverage_deg": result.required_phase.get("coverage_deg"),
             "largest_gap_deg": result.required_phase.get("largest_gap_deg")}
            if result.required_phase else {"available": False}
        ),
        "matching": result.matching or {"matched": 0, "unmatched": 0, "mean_error_deg": None, "max_error_deg": None},
        "hfss_validation": {"performed": False},
        "checks": result.checks,
        "warnings": list(result.warnings),
        "errors": list(result.errors),
        "solver_version": SOLVER_VERSION,
    }


# ------------------------------------------------------------------------------- human-readable report

def render_feasibility_report_md(result: FeasibilityResult, config: FeasibilityConfig):
    lines = [
        "# FEASIBILITY_REPORT.md",
        "",
        "Phase 8A scientific feasibility gate, run against the Phase 7 production unit-cell "
        "library. This report makes no claim beyond what is computed below; a FAIL here is a "
        "library/formulation limitation, not necessarily a software defect (see Sec. 15).",
        "",
        "## 1. Executive summary",
        "",
        f"**Final status: `{result.status}`** (exit code {result.exit_code}).",
        "",
        "## 2. Library configuration",
        "",
        f"- Geometries: {result.library.get('geometry_count')} ({result.library.get('successful_cells')} "
        f"successful, {result.library.get('failed_cells')} failed).",
        f"- Period: {result.library.get('period_x_mm')} x {result.library.get('period_y_mm')} mm.",
        f"- Frequency points: {result.library.get('frequency_points')}.",
        f"- STEP files: {result.library.get('step_files')}.",
        "",
        "## 3. Library integrity",
        "",
        f"{result.checks.get('library_integrity', 'NOT RUN')}.",
        "",
        "## 4. Phase-control quantity",
        "",
        (f"Provided by the library: `{result.phase_source.get('quantity')}` "
         f"(status: `{result.phase_source.get('status')}`)." if result.phase_source else "Not evaluated."),
        "",
        result.phase_source.get("note", "") if result.phase_source else "",
        "",
        "## 5. Phase coverage",
        "",
    ]
    if result.library_phase_coverage:
        c = result.library_phase_coverage
        lines += [
            f"- Circular coverage at {result.frequency.get('requested_GHz')} GHz: "
            f"**{c.get('span_deg'):.2f} deg of 360 deg**.",
            f"- Largest gap: **{c.get('largest_gap_deg'):.2f} deg**.",
            f"- {c.get('n_points')} successful cells; {c.get('n_usable_geometries')} PHASE_USABLE; "
            f"{c.get('n_distinct_phase_states')} distinct phase states.",
            f"- Smallest/median gap: {c.get('smallest_gap_deg'):.2f} / {c.get('median_gap_deg'):.2f} deg.",
        ]
    else:
        lines.append("Not evaluated (gate stopped at an earlier check).")
    lines += ["", "## 6. Frequency dependence", ""]
    if result.freq_coverage:
        cov_vals = [p["coverage_deg"] for p in result.freq_coverage]
        lines.append(
            f"Circular coverage across the {len(result.freq_coverage)} stored frequency points: "
            f"min {min(cov_vals):.2f} deg, max {max(cov_vals):.2f} deg. See "
            "`plots/phase_coverage_vs_frequency.png` and `phase_coverage_vs_frequency.csv`."
        )
    else:
        lines.append("Not evaluated.")
    lines += ["", "## 7. Resonance analysis", ""]
    if result.resonance:
        lines.append(
            f"{result.resonance.get('n_resonant')} cell(s) RESONANT, "
            f"{result.resonance.get('n_non_resonant')} NON_RESONANT (criterion: "
            "`library_classifier.detect_resonance`, unmodified). "
            f"{result.resonance.get('coverage_vs_resonance', {}).get('note', '')}"
        )
    else:
        lines.append("Not evaluated.")
    lines += ["", "## 8. Reference-field compatibility", ""]
    lines.append(
        "Design-specific mode was not invoked for this run (no reference field supplied)."
        if result.reference_field is None else
        f"ok={result.reference_field.get('ok')}, n_samples={result.reference_field.get('n_samples')}."
    )
    lines += ["", "## 9. Required holographic phase", ""]
    if result.required_phase is None:
        lines.append(
            "Not computed: this project has no implemented holographic-synthesis formulation "
            "(`feasibility_required_phase.py`). See Sec. 14/15."
        )
    else:
        lines.append(str(result.required_phase))
    lines += ["", "## 10. Phase matching", ""]
    lines.append(str(result.matching) if result.matching else "Not evaluated.")
    lines += ["", "## 11. Spatial sampling / mapping", ""]
    lines.append(
        "No validated reference-grid -> MTS-grid mapping method exists in this project "
        "(`feasibility_spatial.check_mapping_method_defined`)."
    )
    lines += ["", "## 12. Grating-lobe / period check", ""]
    if result.grating:
        g = result.grating
        lines.append(
            f"period/wavelength = {g.get('period_over_wavelength'):.4f}; onset = "
            f"{g.get('grating_onset_hz', 0) / 1e9:.2f} GHz; is_grating_risk = "
            f"{g.get('is_grating_risk')}."
        )
    else:
        lines.append("Not evaluated.")
    lines += [
        "", "## 13. HFSS validation status", "",
        "**HFSS validation: NOT PERFORMED.** MoM library internal self-consistency was checked "
        "(Sec. 3); this is not HFSS electromagnetic validation.",
        "", "## 14. Final feasibility decision", "",
        f"`{result.status}`.",
        "", "## 15. Limitations", "",
    ]
    lines += [f"- {e}" for e in result.errors] or ["- None beyond what is stated above."]
    lines += ["", "## 16. Recommended next scientific step", ""]
    if result.status == "PHASE_SOURCE_UNVALIDATED":
        lines.append(
            "Resolve, as an explicit architecture decision, whether this library's holographic "
            "phase-control quantity should be S11xx (reflective) or S21xx (transmissive), and "
            "whether a ground plane will be added; only then can this gate's phase-source check "
            "pass. Independently, the measured library phase coverage "
            f"({result.library_phase_coverage.get('span_deg', 0):.2f} deg) is far short of a full "
            "360 deg holographic phase requirement and would itself need library expansion."
        )
    elif result.status == "LIBRARY_INVALID":
        lines.append("Fix the library integrity failures listed in Sec. 3/15 before proceeding.")
    else:
        lines.append("None required by this gate for the current request.")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------------- orchestration

def run_and_write_phase8a(
    config: FeasibilityConfig, output_dir, *, reference_field_path=None, beams=None,
    required_phase_deg=None,
):
    """Sec. 40's full output layout, written from one gate run."""
    output_dir = Path(output_dir)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    result = run_feasibility_gate(
        config, reference_field_path=reference_field_path, beams=beams,
        required_phase_deg=required_phase_deg,
    )

    (output_dir / "feasibility.json").write_text(
        json.dumps(build_feasibility_json(result, config), indent=1) + "\n"
    )
    (output_dir / "feasibility_report.md").write_text(render_feasibility_report_md(result, config))

    if result.table:
        _write_csv(output_dir / "library_phase_coverage.csv",
                   ["L_mm", "W_mm", "phase_deg", "magnitude_db", "S11_real", "S11_imag"], result.table)
        _write_csv(output_dir / "phase_vs_geometry.csv",
                   ["L_mm", "W_mm", "phase_deg", "magnitude_db", "S11_real", "S11_imag"], result.table)
    if result.freq_coverage:
        _write_csv(output_dir / "phase_coverage_vs_frequency.csv",
                   ["frequency_GHz", "coverage_deg", "largest_gap_deg", "phase_min_deg", "phase_max_deg"],
                   result.freq_coverage)

    try:
        if result.table:
            plot_phase_vs_lw(
                result.table,
                sorted({r["L_mm"] for r in result.table}), sorted({r["W_mm"] for r in result.table}),
                result.frequency.get("requested_GHz", 0.0), plots_dir / "library_phase_vs_geometry.png",
            )
            from .phase_extraction import phase_coverage as _cov
            plot_circular_phase_coverage(
                result.table, _cov(result.table), result.frequency.get("requested_GHz", 0.0),
                plots_dir / "library_phase_coverage.png",
            )
        if result.freq_coverage:
            plot_coverage_vs_frequency(result.freq_coverage, plots_dir / "phase_coverage_vs_frequency.png")
        if result.match_results:
            plot_matching_error(result.match_results, plots_dir / "phase_matching_error.png")
    except ImportError:
        pass   # matplotlib not installed; JSON/CSV/report are still written (req: never claim a plot exists)

    return result, output_dir


def _write_csv(path, columns, rows):
    import csv

    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([row.get(c) if isinstance(row, dict) else getattr(row, c) for c in columns])
