"""Phase 8B Sec. 30-32: the 324-geometries-vs-169-distinct-phase-states finding (first surfaced by
Phase 8A's `feasibility_coverage.phase_density`), investigated -- not merely re-asserted -- here.
No solver code is touched; this module only reads already-solved data and groups it.
"""

from dataclasses import dataclass

DEDUP_EPS_DEG = 1e-6   # bit-level dedup, same tolerance as feasibility_coverage.phase_density


@dataclass(frozen=True)
class DistinctStateRow:
    l_mm: float
    w_mm: float
    phase_deg: float
    equivalence_group: int


@dataclass(frozen=True)
class DistinctStateAnalysis:
    n_geometries: int
    n_distinct_phase_states: int
    n_duplicate_records: int
    rows: tuple
    mesh_hypothesis_note: str


def analyze_distinct_states(table, mesh_dx_mm=None):
    """Sec. 31: assigns every row of ``table`` (`feasibility_coverage.library_table_at_frequency`'s
    schema, or any list of dicts with L_mm/W_mm/phase_deg) an ``equivalence_group`` id -- rows
    whose phase agrees to within `DEDUP_EPS_DEG` share a group. Sec. 32: this is explicitly a
    *phase*-diversity count, not a claim about *geometric* diversity (``n_geometries`` is reported
    alongside, never conflated)."""
    if not table:
        raise ValueError("analyze_distinct_states: empty table")
    ordered = sorted(table, key=lambda r: r["phase_deg"])
    rows = []
    group = 0
    prev_phase = None
    for r in ordered:
        if prev_phase is None or (r["phase_deg"] - prev_phase) > DEDUP_EPS_DEG:
            group += 1
        rows.append(DistinctStateRow(r["L_mm"], r["W_mm"], r["phase_deg"], group))
        prev_phase = r["phase_deg"]

    n_distinct = group
    n_dup = len(table) - n_distinct

    note = (
        f"{n_dup} of {len(table)} geometries share a bit-identical (within {DEDUP_EPS_DEG:g} deg) "
        f"stored phase with at least one other geometry ({n_distinct} distinct phase states). "
        "This is consistent with, not contradictory to, the mesh-coarseness limitation "
        "`LIBRARY_REPORT.md` Sec. 16 already documents (32x32 cells over the 4.8 mm period gives "
        "dx = dy = 0.15 mm, `library_schema.DEFAULT_MESH`): a sufficiently large L, when "
        "discretized to whole 0.15 mm mesh cells, can round the *same* number of metal cells for "
        "two nominally different W values whenever the physically relevant (x-polarized) response "
        "is dominated by L rather than W -- an ordinary staircase-quantization effect on the "
        "L axis's neighbourhood, not an independent bug. This module does not modify the mesh or "
        "solver to test this hypothesis further (out of Phase 8B's scope, Sec. 3); it is reported "
        "descriptively, as Sec. 30 requires, not asserted as certain."
    )
    return DistinctStateAnalysis(
        n_geometries=len(table), n_distinct_phase_states=n_distinct, n_duplicate_records=n_dup,
        rows=tuple(rows), mesh_hypothesis_note=note,
    )


def write_distinct_state_csv(analysis, path):
    import csv

    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["L_mm", "W_mm", "phase_deg", "equivalence_group"])
        for row in analysis.rows:
            writer.writerow([repr(row.l_mm), repr(row.w_mm), repr(row.phase_deg), row.equivalence_group])
    return path
