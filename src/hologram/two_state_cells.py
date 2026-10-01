"""Phase 15 Part E: two-state (1-bit) unit-cell selection for the Coded design method (Mishra,
Barad, Sarkar, IEEE TVT 2026). A 1-bit coded hologram needs two reflective cell states with
reflection phases ~180deg apart and comparable magnitude -- Mishra's own pair is a full square
patch vs. a patch with an annular slot (`src.geometry.AnnularSlotPatch`), but the achieved
separation is a MEASUREMENT of whatever two real candidates are chosen, never assumed to be
exactly 180deg (`separation_warning` reports the shortfall honestly, per this project's own
scientific-honesty convention -- see CLAUDE.md)."""

from dataclasses import dataclass

from ..phase import principal_phase_deg, wrapped_difference_deg

IDEAL_SEPARATION_DEG = 180.0
SEPARATION_WARNING_THRESHOLD_DEG = 30.0  # how far from 180deg before this is flagged, not silently used


@dataclass(frozen=True)
class TwoStateCell:
    """One candidate reflective state: a human label, its own record id (traceable back to a
    library entry or an uploaded row), the achieved complex reflection coefficient, and whatever
    geometry parameters identify it for export (`{"kind": "rectangle", "L_mm", "W_mm"}` or
    `{"kind": "annular_slot", "a_mm", "b_mm", "R_mm", "r_mm"}`)."""

    label: str
    record_id: str
    value: complex
    geometry: dict


@dataclass(frozen=True)
class TwoStateChoice:
    state_0: TwoStateCell
    state_1: TwoStateCell
    achieved_separation_deg: float


def candidates_from_square_library(candidates):
    """`candidates`: `(record_id, L_mm, W_mm, library_id, complex_value)` tuples (the format
    `phase8d_rev.response_matching.nearest_response_match`/`hmsa_synthesis.map_control_to_geometry_library`
    already use) -- wraps each as a `TwoStateCell` labelled by its own square side length."""
    return [
        TwoStateCell(
            label=f"square patch L={l_mm:g}mm", record_id=record_id, value=complex(value),
            geometry={"kind": "rectangle", "L_mm": l_mm, "W_mm": w_mm},
        )
        for record_id, l_mm, w_mm, _library_id, value in candidates
    ]


def candidates_from_annular_slot_library(cell_docs, frequency_ghz):
    """`cell_docs`: `library_sweep.build_annular_slot_library_cells`'s own return list -- wraps
    each successfully-solved cell, at the frequency sample nearest `frequency_ghz`, as a
    `TwoStateCell`."""
    out = []
    for doc in cell_docs:
        if doc["solver_status"] != "OK" or not doc["frequency_GHz"]:
            continue
        idx = min(range(len(doc["frequency_GHz"])), key=lambda i: abs(doc["frequency_GHz"][i] - frequency_ghz))
        entry = doc["S11"][idx]
        value = complex(entry["real"], entry["imag"])
        out.append(
            TwoStateCell(
                label=f"annular slot R={doc['R_mm']:g}mm", record_id=doc["cell_id"], value=value,
                geometry={
                    "kind": "annular_slot", "a_mm": doc["a_mm"], "b_mm": doc["b_mm"],
                    "R_mm": doc["R_mm"], "r_mm": doc["r_mm"],
                },
            )
        )
    return out


def best_separated_pair(candidates):
    """Brute-force O(n^2) search (library sizes here are at most a few dozen to ~100 candidates --
    this project never guesses at a smarter search without measuring that it is needed) for the
    pair of `candidates` whose reflection PHASES are closest to 180deg apart -- real, measured
    values only. `candidates` may mix cell families (e.g. square-patch AND annular-slot
    `TwoStateCell`s together), which is exactly how Mishra's own two-state pair (a plain patch vs.
    one WITH an annular slot) is found in practice."""
    candidates = list(candidates)
    if len(candidates) < 2:
        raise ValueError(f"best_separated_pair needs at least 2 candidates, got {len(candidates)}")
    best = None
    for i in range(len(candidates)):
        phase_i = principal_phase_deg(candidates[i].value)
        for j in range(i + 1, len(candidates)):
            phase_j = principal_phase_deg(candidates[j].value)
            sep = abs(wrapped_difference_deg(phase_i, phase_j))
            if best is None or sep > best[0]:
                best = (sep, candidates[i], candidates[j])
    sep, a, b = best
    return TwoStateChoice(state_0=a, state_1=b, achieved_separation_deg=float(sep))


def two_state_from_values(label_0, value_0, label_1, value_1, geometry_0=None, geometry_1=None):
    """Part E option (c): an HFSS-measured two-state table supplied directly by the user, rather
    than searched from one of this project's own libraries."""
    a = TwoStateCell(label_0, "uploaded_state_0", complex(value_0), geometry_0 or {})
    b = TwoStateCell(label_1, "uploaded_state_1", complex(value_1), geometry_1 or {})
    sep = abs(wrapped_difference_deg(principal_phase_deg(a.value), principal_phase_deg(b.value)))
    return TwoStateChoice(state_0=a, state_1=b, achieved_separation_deg=float(sep))


def separation_warning(choice: TwoStateChoice, threshold_deg=SEPARATION_WARNING_THRESHOLD_DEG):
    """`None` if the achieved separation is within `threshold_deg` of the ideal 180deg; otherwise a
    plain-language warning string naming the actual shortfall -- never silently accepted as "close
    enough" (this project's own scientific-honesty convention, CLAUDE.md)."""
    shortfall = abs(IDEAL_SEPARATION_DEG - choice.achieved_separation_deg)
    if shortfall <= threshold_deg:
        return None
    return (
        f"These two cell states achieve only {choice.achieved_separation_deg:.1f}deg of phase "
        f"separation (ideal: 180deg) -- {shortfall:.1f}deg short. 1-bit coding will lose more "
        f"beam efficiency than an ideal pair; consider a different library or geometry pair."
    )
