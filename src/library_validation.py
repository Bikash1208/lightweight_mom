"""Phase 7: the deterministic HFSS validation subset (req. 22) and internal consistency checks
(req. 20's "CSV and JSON internally consistent") -- this module makes **no** claim that HFSS
validation has occurred; `LIBRARY_VALIDATION.md` says so explicitly (req. 22's "do not claim HFSS
agreement until actual HFSS data are supplied", req. 31 item 18).
"""

import csv
import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SubsetEntry:
    cell_id: str
    l_mm: float
    w_mm: float
    reasons: tuple


def select_hfss_validation_subset(config, table, classifications):
    """Deterministic selection covering the categories req. 22 lists, from data already computed
    (``table``: `phase_extraction.phase_vs_geometry_table`; ``classifications``:
    `library_classifier.classify_cell` results) -- no random sampling, no manual picks. Each
    reason is evaluated independently and a cell selected by more than one keeps every reason it
    earned (e.g. the smallest patch may also be the non-resonant example)."""
    l_values = sorted(config.l_values_mm)
    w_values = sorted(config.w_values_mm)
    l_min, l_max = l_values[0], l_values[-1]
    l_mid = l_values[len(l_values) // 2]
    w_mid = w_values[len(w_values) // 2]

    reasons_by_lw = {}

    def add(l_mm, w_mm, reason):
        reasons_by_lw.setdefault((l_mm, w_mm), []).append(reason)

    add(l_min, l_min, "small patch (smallest L and W in the sweep)")
    add(l_max, l_max, "large patch (largest L and W in the sweep)")
    add(l_mid, w_mid, "approximately square, medium-sized patch (sweep midpoint)")
    add(l_min, l_max, "small/large aspect ratio (narrowest L, widest W)")
    add(l_max, l_min, "small/large aspect ratio (widest L, narrowest W)")

    ok = [c for c in classifications if c.solver_status == "OK"]
    resonant = [c for c in ok if c.resonance is not None and c.resonance.is_resonant]
    non_resonant = [c for c in ok if c.resonance is not None and not c.resonance.is_resonant]
    if resonant:
        best = max(resonant, key=lambda c: (c.resonance.resonance_metric, c.l_mm, c.w_mm))
        add(best.l_mm, best.w_mm, "resonant cell (highest resonance prominence in the library)")
    if non_resonant:
        first = min(non_resonant, key=lambda c: (c.l_mm, c.w_mm))
        add(first.l_mm, first.w_mm, "non-resonant cell (deterministic: smallest (L, W) with no detected resonance)")

    if table:
        by_phase = sorted(table, key=lambda r: (r["phase_deg"], r["L_mm"], r["W_mm"]))
        n = len(by_phase)
        quartile_indices = sorted({0, n // 4, n // 2, (3 * n) // 4, n - 1})
        for q in quartile_indices:
            row = by_phase[q]
            add(
                row["L_mm"],
                row["W_mm"],
                f"phase-region representative (S11xx phase quartile sample, phase={row['phase_deg']:.1f} deg)",
            )

    cell_by_lw = {(c.l_mm, c.w_mm): c.cell_id for c in classifications}
    subset = [
        SubsetEntry(cell_by_lw[(l_mm, w_mm)], l_mm, w_mm, tuple(reasons))
        for (l_mm, w_mm), reasons in sorted(reasons_by_lw.items())
        if (l_mm, w_mm) in cell_by_lw
    ]
    return subset


def write_hfss_validation_subset_csv(subset, path):
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["cell_id", "L_mm", "W_mm", "reason_for_selection"])
        for entry in subset:
            writer.writerow([entry.cell_id, repr(entry.l_mm), repr(entry.w_mm), "; ".join(entry.reasons)])
    return path


# ------------------------------------------------------------------------ internal consistency


@dataclass(frozen=True)
class ConsistencyReport:
    ok: bool
    n_json_records: int
    n_csv_rows: int
    mismatches: tuple


def check_json_csv_consistency(json_path, csv_path, rtol=1e-9):
    """Req. 20: `unit_cell_library.csv` (flat) and `unit_cell_library.json` (nested by cell,
    `unit_cell_library.write_library_json`) must agree numerically wherever they overlap
    (real/imag/mag/mag_db/phase_deg of S11/S21/S12/S22 per (cell, frequency)) -- both read
    independently from disk, never compared against in-memory state that produced either."""
    doc = json.loads(Path(json_path).read_text(encoding="utf-8"))
    with open(csv_path, newline="", encoding="utf-8") as fh:
        csv_rows = list(csv.DictReader(fh))
    by_key = {(row["cell_id"], row["f_GHz"]): row for row in csv_rows}

    mismatches = []
    n_json_records = 0
    for cell in doc["cells"]:
        if cell["solver_status"] != "OK":
            continue
        for i, f_ghz in enumerate(cell["frequency_GHz"]):
            n_json_records += 1
            row = by_key.get((cell["cell_id"], repr(f_ghz)))
            if row is None:
                candidates = [
                    r for r in csv_rows if r["cell_id"] == cell["cell_id"] and abs(float(r["f_GHz"]) - f_ghz) < 1e-6
                ]
                row = candidates[0] if candidates else None
            if row is None:
                mismatches.append(f"{cell['cell_id']} @ {f_ghz:g} GHz: missing from CSV")
                continue
            for family in ("S11", "S21", "S12", "S22"):
                entry = cell[family][i]
                for field_ in ("real", "imag", "mag", "mag_db", "phase_deg"):
                    expected = entry[field_]
                    got = float(row[f"{family}_{field_}"])
                    if abs(got - expected) > rtol * max(abs(expected), 1.0):
                        mismatches.append(
                            f"{cell['cell_id']} @ {f_ghz:g} GHz: {family}_{field_} JSON={expected!r} CSV={got!r}"
                        )
    return ConsistencyReport(len(mismatches) == 0, n_json_records, len(csv_rows), tuple(mismatches))


# ------------------------------------------------------------------------------------ report


def render_library_validation_md(config, subset, consistency, build):
    lines = [
        "# LIBRARY_VALIDATION.md",
        "",
        "## HFSS validation status",
        "",
        "**No HFSS comparison has been performed.** Every number in this library is produced by "
        "this project's own lightweight periodic MoM solver (Phases 1-6A, unmodified). The "
        "subset below is a deterministic, documented *recommendation* of which cells to validate "
        "in HFSS first; it is not evidence that any of them agree with HFSS, because no HFSS "
        "simulation has been run or imported for this library (req. 22).",
        "",
        "## Recommended HFSS validation subset",
        "",
        f"{len(subset)} of {config.n_geometries} cells, selected deterministically (no random "
        "sampling, no manual picks) to cover: a small patch, a large patch, an approximately "
        "square medium patch, both extreme aspect ratios, the most resonant cell found, one "
        "non-resonant cell, and phase-quartile representatives at the design frequency. Full "
        "table: `HFSS_VALIDATION_SUBSET.csv`.",
        "",
        "| cell_id | L_mm | W_mm | reason(s) |",
        "|---|---|---|---|",
    ]
    for entry in subset:
        lines.append(f"| {entry.cell_id} | {entry.l_mm:g} | {entry.w_mm:g} | {'; '.join(entry.reasons)} |")

    lines += [
        "",
        "## Internal consistency (JSON vs. CSV)",
        "",
        f"`unit_cell_library.json` ({consistency.n_json_records} records) vs. "
        f"`unit_cell_library.csv` ({consistency.n_csv_rows} rows): "
        f"{'**consistent**' if consistency.ok else f'**{len(consistency.mismatches)} mismatch(es) found**'}.",
    ]
    if not consistency.ok:
        lines.append("")
        for m in consistency.mismatches[:20]:
            lines.append(f"- {m}")
        if len(consistency.mismatches) > 20:
            lines.append(f"- (+{len(consistency.mismatches) - 20} more)")

    lines += [
        "",
        "## Build completeness",
        "",
        f"- geometries: {config.n_geometries} expected, {build.n_ok_cells} solved OK, "
        f"{build.n_invalid_cells} invalid/skipped.",
        f"- computed (this run): {build.computed}; reused from checkpoint: {build.reused}.",
    ]
    return "\n".join(lines) + "\n"
