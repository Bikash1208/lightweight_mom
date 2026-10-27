"""Phase 28 item 3: evaluation of a diode-switched CP rotating-axis cell from its normal-incidence TENSOR SHEET (the Phase 26 B screening method, additive).

For one state of a cell (a `UnitCell` with lumped diode loads) the 2 x 2 patch-layer sheet admittance `Y_s` is extracted exactly as in `dispersion.tensor_sheet`
(`reflection_matrix` -> `sheet_admittance_tensor_from_gamma`; the lumped loads enter through the MoM). From it:

* mean susceptance `b_bar`, anisotropy `delta` (half the split of the principal susceptances, S) and axis angle (`tensor_sheet.decompose`);
* the TM0 phase constant `beta_bar` of the mean sheet (`sheet_model.find_beta_from_sheet_admittance`);
* `dX_aniso` [ohm]: the difference of the TM surface reactances `X = Im(1 / (Y + Y_slab,TM(beta_bar)))` along the two principal axes
  (`leakage_model.state_surface_reactances`) -- the contrast measure of the item-2 feasibility report (tensor row);
* the dissipation `g_max` (largest principal conductance, S) and the ratio `g_max / |delta|`.

Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET (diode values), VIA_OMITTED_APPROXIMATION. A number from this module is a screening value of the
sheet model, not a dispersion (alpha, beta) of a row of cells: that is the state library (item 4).
"""

import warnings

import numpy as np

from ..basis_functions import MeshSettings
from ..dispersion.leakage_model import leakage_alpha_first_order_tensor, state_surface_reactances
from ..dispersion.sheet_model import find_beta_from_sheet_admittance
from ..dispersion.tensor_sheet import decompose, reflection_matrix, sheet_admittance_tensor_from_gamma
from ..frequency import free_space_wavenumber


def sheet_tensor(unit_cell, frequency_hz, mesh_n, floquet_factor=1):
    """`Y_s` (2 x 2 complex, S) of a grounded cell at `frequency_hz` on a `mesh_n` x `mesh_n` rooftop mesh; the Floquet cutoff is `floquet_factor * mesh_n`
    (1: the Phase 26 screening setting; 2: the converged cutoff of the Phase 20 R4 audit for cells with sub-wavelength lumped loads)."""
    mesh = MeshSettings(mesh_n, mesh_n, int(floquet_factor * mesh_n))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = reflection_matrix(unit_cell, frequency_hz, mesh)
    return sheet_admittance_tensor_from_gamma(g, unit_cell.substrate, frequency_hz)


def design_loading(frequency_hz, substrate, beta_over_k0=1.06):
    """The mean sheet susceptance `b` [S] whose TM0 constant is `beta_over_k0 k0` (the mean loading of the item-2 contrast derivation, root-solved)."""
    from scipy.optimize import brentq

    k0 = free_space_wavenumber(frequency_hz)

    def f(b):
        return (
            find_beta_from_sheet_admittance(1j * b * 1e-3, substrate, frequency_hz).beta_rad_per_m / k0 - beta_over_k0
        )

    return brentq(f, 0.5, 20.0, xtol=1e-9) * 1e-3


def principal_metrics(y, substrate, frequency_hz, beam_theta_deg=30.0, beta_design=1.06):
    """Screening metrics of a sheet tensor (module docstring). The susceptance tensor is `Im(Y)`; `Y` need not be exactly symmetric (reported)."""
    yi, delta, axis, asym = decompose(y)
    b_bar = float(np.imag(yi))
    out = {
        "b_bar_ms": b_bar * 1e3,
        "delta_ms": delta * 1e3,
        "axis_deg": float(np.degrees(axis)),
        "asymmetry_ms": float(asym) * 1e3,
        "g_mean_ms": float(np.real(yi)) * 1e3,
    }
    # principal conductances (real part of the tensor): the largest eigenvalue
    gre = np.real(0.5 * (np.asarray(y) + np.asarray(y).T))
    out["g_max_ms"] = float(np.max(np.linalg.eigvalsh(gre))) * 1e3
    out["g_over_delta"] = out["g_max_ms"] / max(out["delta_ms"], 1e-12)
    k0 = free_space_wavenumber(frequency_hz)
    try:
        mode = find_beta_from_sheet_admittance(1j * b_bar, substrate, frequency_hz)
        beta = mode.beta_rad_per_m
        xs = state_surface_reactances([1j * (b_bar + delta), 1j * (b_bar - delta)], beta, substrate, frequency_hz)
        out.update(
            {
                "beta_over_k0": beta / k0,
                "dx_aniso_ohm": float(abs(xs[0] - xs[1])),
                "x_principal_ohm": [float(xs[0]), float(xs[1])],
            }
        )
        # the leakage this tensor sheet delivers to a beam at `beam_theta_deg` (Phase 26 closed form) and the same anisotropy on the DESIGN mean loading
        out["alpha_beam_np_per_m"] = float(
            leakage_alpha_first_order_tensor(
                beta, 1j * b_bar, delta, np.radians(beam_theta_deg), substrate, frequency_hz
            )
        )
        b_des = design_loading(frequency_hz, substrate, beta_design)
        beta_des = find_beta_from_sheet_admittance(1j * b_des, substrate, frequency_hz).beta_rad_per_m
        xd = state_surface_reactances([1j * (b_des + delta), 1j * (b_des - delta)], beta_des, substrate, frequency_hz)
        out["dx_aniso_at_design_loading_ohm"] = float(abs(xd[0] - xd[1]))
    except Exception as exc:  # noqa: BLE001 - no bound mode at this loading: reported, not hidden
        out.update(
            {
                "beta_over_k0": None,
                "dx_aniso_ohm": None,
                "x_principal_ohm": None,
                "note": f"NO_BOUND_MODE: {type(exc).__name__}",
            }
        )
    return out


def evaluate_states(cell, kit, frequency_hz, mesh_n, states=None):
    """`{state_index: metrics}` for the cell's states (`cell.states()`, `cell.unit_cell(state, kit)`), plus the tensors."""
    states = cell.states() if states is None else states
    out = {}
    for idx, st in enumerate(states):
        uc = cell.unit_cell(st, kit)
        y = sheet_tensor(uc, frequency_hz, mesh_n)
        m = principal_metrics(y, uc.substrate, frequency_hz)
        m["y_s"] = [[[float(v.real), float(v.imag)] for v in row] for row in y]
        out[idx] = m
    return out
