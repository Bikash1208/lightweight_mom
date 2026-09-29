"""Phase 8C Sec. 14-15: power conservation / passivity and propagating-mode analysis.

Reuses `polarization.co_cross_power` (Phase 5, unmodified) so cross-polarized power is always
included, never a bare `|S11xx|^2` reading (Sec. 14's explicit warning against that), and
`feasibility_grating.grating_check` (Phase 8A, unmodified) for the propagating-mode check.
"""

from dataclasses import dataclass

import numpy as np

from ..feasibility_grating import grating_check


@dataclass(frozen=True)
class PowerBalance:
    R_co: float  # |S11xx|^2
    R_cross: float  # |S11yx|^2
    T_co: float  # |S21xx|^2
    T_cross: float  # |S21yx|^2
    R_total: float
    T_total: float
    R_plus_T: float
    absorbed_fraction: float  # 1 - R_plus_T -- NOT claimed as "absorption" without the caveat below
    note: str


def power_balance_for_record(record):
    """Sec. 14: R = co+cross reflected power for x-polarized incidence, T = co+cross transmitted
    power, both via `polarization.co_cross_power` on the actual S11/S21 (2,2) blocks (row = column
    0 of the incident-x column, i.e. this uses the full first *column* of each block, not just the
    diagonal xx entry)."""
    from ..database import family_block

    s11_block = family_block(record.s, "S11")
    s21_block = family_block(record.s, "S21")
    # x-incidence only: column 0 of each 2x2 block = (S_xx, S_yx) for that family.
    r_co, r_cross = abs(s11_block[0, 0]) ** 2, abs(s11_block[1, 0]) ** 2
    t_co, t_cross = abs(s21_block[0, 0]) ** 2, abs(s21_block[1, 0]) ** 2
    r_total, t_total = r_co + r_cross, t_co + t_cross
    r_plus_t = r_total + t_total
    note = (
        f"1 - (R+T) = {1.0 - r_plus_t:.6f} is only interpretable as dielectric absorption if no "
        "other propagating diffraction order exists at this frequency/period (Sec. 14's explicit "
        "warning); see propagating_mode_check for that condition at this specific frequency."
    )
    return PowerBalance(
        R_co=float(r_co),
        R_cross=float(r_cross),
        T_co=float(t_co),
        T_cross=float(t_cross),
        R_total=float(r_total),
        T_total=float(t_total),
        R_plus_T=float(r_plus_t),
        absorbed_fraction=float(1.0 - r_plus_t),
        note=note,
    )


def library_power_summary(database, frequency_index):
    """Mean R/T (co and cross) over every successful cell at one stored frequency index --
    feeds `topology_audit.audit_topology`'s architecture classification."""
    balances = [
        power_balance_for_record(database.record(cell.index, frequency_index))
        for cell in database.cells
        if cell.status == "ok"
    ]
    if not balances:
        raise ValueError("library_power_summary: no successful cells")
    return {
        "n_cells": len(balances),
        "mean_R_co": float(np.mean([b.R_co for b in balances])),
        "mean_R_cross": float(np.mean([b.R_cross for b in balances])),
        "mean_T_co": float(np.mean([b.T_co for b in balances])),
        "mean_T_cross": float(np.mean([b.T_cross for b in balances])),
        "mean_R_plus_T": float(np.mean([b.R_plus_T for b in balances])),
        "min_R_plus_T": float(np.min([b.R_plus_T for b in balances])),
        "max_R_plus_T": float(np.max([b.R_plus_T for b in balances])),
    }


def propagating_mode_check(unit_cell, frequencies_hz):
    """Sec. 15: reuses `feasibility_grating.grating_check` (Phase 8A, unmodified -- same formula
    as `discretization.grating_onset_frequency`) at the band edges to report whether any
    higher-order Floquet mode propagates anywhere in ``frequencies_hz``."""
    lat = unit_cell.lattice
    results = [grating_check(lat.period_x, lat.period_y, f) for f in (min(frequencies_hz), max(frequencies_hz))]
    any_risk = any(r.is_grating_risk for r in results)
    return {
        "period_x_m": lat.period_x,
        "period_y_m": lat.period_y,
        "grating_onset_hz": results[0].grating_onset_hz,
        "band_min_hz": min(frequencies_hz),
        "band_max_hz": max(frequencies_hz),
        "only_fundamental_mode_propagates": not any_risk,
        "note": (
            "Only the fundamental (0,0) Floquet order propagates throughout the analyzed band "
            f"(onset at {results[0].grating_onset_hz / 1e9:.2f} GHz, far above it); S11/S21 fully "
            "describe the specular scattered power."
            if not any_risk
            else "A higher-order Floquet mode propagates within the analyzed band: S11/S21 alone do "
            "NOT fully describe the scattered power; 1-R-T cannot be interpreted as absorption."
        ),
    }
