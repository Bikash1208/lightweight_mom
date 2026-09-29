"""Sec. 15: complex-response nearest-neighbour search. Amplitude, phase and complex error are
always reported separately (Sec. 15's explicit "do not collapse these into a single score unless
a physically justified objective function is explicitly defined" -- none is defined anywhere in
this project, so none is invented here); `complex_error` (the plain |difference|) is used only to
*order* candidates, never to declare a single "winner" without also showing the breakdown.
"""

from dataclasses import dataclass

from ..phase import principal_phase_deg, wrapped_difference_deg


@dataclass(frozen=True)
class ResponseMatch:
    record_id: str
    L_mm: float
    W_mm: float
    library_id: str
    complex_value: complex
    complex_error: float  # |required - candidate|
    amplitude_error: float  # |  |required| - |candidate|  |
    phase_error_deg: float  # circular wrapped difference, degrees


def complex_distance(required, candidate):
    return abs(complex(required) - complex(candidate))


def nearest_response_match(required_complex, candidates):
    """``candidates``: an iterable of (record_id, L_mm, W_mm, library_id, complex_value). Returns
    every candidate's `ResponseMatch` sorted by ``complex_error`` ascending -- the *whole* list,
    not a single chosen entry (mirrors `unit_cell_library.UnitCellLibrary.find_phase_candidates`'s
    "all matches, never a single chosen winner" policy, extended to complex-plane distance)."""
    required = complex(required_complex)
    req_mag, req_phase = abs(required), float(principal_phase_deg(required))
    matches = []
    for record_id, l_mm, w_mm, library_id, value in candidates:
        value = complex(value)
        matches.append(
            ResponseMatch(
                record_id=record_id,
                L_mm=l_mm,
                W_mm=w_mm,
                library_id=library_id,
                complex_value=value,
                complex_error=complex_distance(required, value),
                amplitude_error=abs(abs(value) - req_mag),
                phase_error_deg=float(abs(wrapped_difference_deg(float(principal_phase_deg(value)), req_phase))),
            )
        )
    return sorted(matches, key=lambda m: m.complex_error)
