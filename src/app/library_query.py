"""Phase 13 Part 3: turns a registered library's stored responses at one frequency into the
``(record_id, L_mm, W_mm, library_id, complex_value)`` candidate tuples
`phase8d_rev.response_matching.nearest_response_match` (and, through it,
`reflectarray.designer.assign_cells_from_library`) expects -- shared by the reflectarray and HMSA
routers so there is one place this conversion happens.
"""


def library_candidates_at_frequency(lib, library_id, frequency_hz, family="S11"):
    """``lib``: a loaded `unit_cell_library.UnitCellLibrary`. Only cells with an exact stored
    record at ``frequency_hz`` are included (no interpolation -- `UnitCellLibrary.response`'s own
    rule); an empty result is a real, reportable outcome, not silently padded."""
    frequency_ghz = frequency_hz / 1e9
    candidates = []
    for cell in lib.all_cells():
        if cell["solver_status"] != "OK":
            continue
        try:
            entry = lib.response(cell["L_mm"], cell["W_mm"], frequency_ghz, family=family)
        except ValueError:
            continue
        candidates.append(
            (cell["cell_id"], cell["L_mm"], cell["W_mm"], library_id, complex(entry["real"], entry["imag"]))
        )
    return candidates
