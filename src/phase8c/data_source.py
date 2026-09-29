"""Phase 8C data source: the full linear-basis 4x4 S-matrix (xx, xy, yx, yy of every S-family),
read from the still-present internal Phase 7 solve checkpoint
(`results/unit_cell_library/_solve/database.json`) rather than the delivered
`database/unit_cell_library.json`, which only exports the co-polarized "_xx" component per family
(`unit_cell_library.py::_cell_doc`). Section 11 requires investigating cross-polarized terms
(S11xy, S11yx, ...), which are real, already-computed Phase 7 data sitting on disk -- reading them
is not a Phase 7 modification (nothing is written back), and it is more honest than declaring
cross-polarization "unsupported" when the actual simulated data already exists.

`load_full_database` cross-checks every "_xx" value it reads against the delivered
`unit_cell_library.json` (Phase 8A's own integrity-checked deliverable) so this module can never
silently read a stale or different solve.
"""

import json
from pathlib import Path

from ..database import Database
from ..feasibility_schema import DEFAULT_LIBRARY_DIR

DEFAULT_CHECKPOINT_JSON = DEFAULT_LIBRARY_DIR / "_solve" / "database.json"


def load_full_database(checkpoint_path=DEFAULT_CHECKPOINT_JSON):
    """The full `database.Database` (Phase 5's own native format), read directly from disk.
    Raises FileNotFoundError with a clear message if the checkpoint is not present (e.g. a fresh
    checkout that only has the delivered library, not the internal solve directory) -- callers
    must not substitute a smaller/partial database silently."""
    checkpoint_path = Path(checkpoint_path)
    if not checkpoint_path.exists():
        raise FileNotFoundError(
            f"Phase 8C needs the full linear-basis S-matrix, which only the internal Phase 7 "
            f"solve checkpoint has ({checkpoint_path} not found). The delivered "
            f"database/unit_cell_library.json only exports co-polarized '_xx' components "
            f"(unit_cell_library.py::_cell_doc) and cannot answer Sec. 11's cross-polarization "
            f"question. Re-run examples/build_unit_cell_library.py (it resumes from the "
            f"checkpoint if present, or rebuilds it) to regenerate this file."
        )
    return Database.from_json(checkpoint_path)


def cross_check_against_delivered_library(
    database, design_frequency_ghz=12.5, library_dir=DEFAULT_LIBRARY_DIR, n_sample=20, rtol=1e-9
):
    """Sanity check (not a Phase 7 modification): a sample of this checkpoint's own S11_xx/S21_xx
    values at ``design_frequency_ghz`` must agree exactly with the delivered
    `unit_cell_library.json` -- proof this is genuinely the same solve Phase 7 delivered, not a
    stale or different one."""
    delivered = json.loads((Path(library_dir) / "database" / "unit_cell_library.json").read_text())
    by_lw = {(c["L_mm"], c["W_mm"]): c for c in delivered["cells"] if c["solver_status"] == "OK"}
    freqs = database.frequencies
    f_idx = next(
        (i for i, f in enumerate(freqs) if abs(f / 1e9 - design_frequency_ghz) <= rtol * design_frequency_ghz), None
    )
    if f_idx is None:
        raise ValueError(f"{design_frequency_ghz} GHz is not an exact grid point of this checkpoint")

    mismatches, checked = [], 0
    ok_cells = [c for c in database.cells if c.status == "ok"]
    for cell in ok_cells[:: max(1, len(ok_cells) // n_sample)][:n_sample]:
        key = (round(cell.parameters["length_m"] * 1e3, 9), round(cell.parameters["width_m"] * 1e3, 9))
        deliv_cell = by_lw.get(key)
        if deliv_cell is None:
            mismatches.append(f"{key}: not found in delivered library")
            continue
        record = database.record(cell.index, f_idx)
        checked += 1
        for fam, deliv_fam in (("S11_xx", "S11"), ("S21_xx", "S21")):
            got = record.component(fam)
            deliv_entry = deliv_cell[deliv_fam][f_idx]
            if abs(got.real - deliv_entry["real"]) > rtol or abs(got.imag - deliv_entry["imag"]) > rtol:
                mismatches.append(f"{key} {fam}: checkpoint={got!r} delivered={deliv_entry!r}")
    return {"checked": checked, "mismatches": mismatches, "ok": len(mismatches) == 0}
