"""Phase 28 item 3: CELL_CANDIDATES.md generated from the characterisation JSONs only (`phase28_cell_characterise.py`); every number printed is a value of those files
(`tests/test_phase28_cell_report.py`). Prose carries no numbers of its own.

    python scripts/research/phase28_cell_report.py <out.md> <cell.json> [<cell.json> ...] [--search-summary <summary.json>]
"""

import json
import sys
from pathlib import Path


def f(v, nd=2):
    return "n/a" if v is None else f"{v:.{nd}f}"


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out


def _row(rows, **kw):
    kw.setdefault("floquet_factor", 1)
    return next(r for r in rows if all(r.get(k) == v for k, v in kw.items()))


def candidate(doc):
    L = []
    a = L.append
    sp, d = doc["spec"], doc["derived"]
    p = sp["params"]
    a(f"## Candidate {sp['name']}")
    a("")
    a(
        "Hub-and-sector cell (`coded_cell_candidates`, the family is in the table): a hub and eight sectors, each joined to the hub by a neck that carries one diode. Mode `bowtie`: state k = the two necks on axis k conduct (the hub and the two opposite sectors form a bow-tie along that axis), the other six are open. "
        "Mode `notch`: every neck conducts except the two of the perpendicular axis (the metal fill stays high, two opposite sectors float). The four states are the same cell turned in 45 degree steps, so the sheet tensor rotates with the state."
    )
    a("")
    keys = [k for k in p if k != "mode"]
    L += table(
        ["family", "period [mm]", "diode kit"] + [k + " [mm]" for k in keys],
        [(sp.get("family", "sector_hub"), f(sp["period_mm"], 1), doc["kit"]) + tuple(f(p[k] * 1e3, 2) for k in keys)],
    )
    a("")
    a("### Mesh convergence (normal-incidence tensor sheet, 10.0 GHz)")
    a("")
    rows = []
    for mesh in sp["meshes"]:
        for s in range(4):
            r = _row(doc["rows"], state=s, mesh=mesh, f_ghz=10.0, scale=None)
            rows.append(
                (
                    mesh,
                    s,
                    f(r["b_bar_ms"], 3),
                    f(r["delta_ms"], 3),
                    f(r["beta_over_k0"], 4),
                    f(r["dx_aniso_ohm"], 1),
                    f(r["dx_aniso_at_design_loading_ohm"], 1),
                    f(r["alpha_beam_np_per_m"], 3),
                    f(r["g_max_ms"], 3),
                    f(r["axis_deg"], 1),
                )
            )
    L += table(
        [
            "mesh",
            "state",
            "mean loading b [mS]",
            "anisotropy delta [mS]",
            "beta/k0",
            "dX_aniso own loading [ohm]",
            "dX_aniso design loading [ohm]",
            "alpha at the beam [Np/m]",
            "largest conductance [mS]",
            "principal axis [deg]",
        ],
        rows,
    )
    a("")
    a(
        f"Convergence is judged on the relative change between mesh {d['mesh_prev']} and mesh {d['mesh_last']}, largest over the four states, against the tolerance of the table (an ASSUMED criterion: {f(d['convergence_tolerance'], 2)} on the anisotropy, the loading and the reactance difference; the leakage at the beam grows as the square of the anisotropy, so its tolerance is twice that):"
    )
    a("")
    L += table(
        ["metric", "largest relative change", "tolerance", "converged"],
        [
            (k, f(v["max_relative_change"], 4), f(v["tolerance"], 2), "yes" if v["converged"] else "NO")
            for k, v in d["convergence"].items()
        ],
    )
    a("")
    a(
        "All metrics converged: "
        + ("yes." if d["convergence_all"] else "NO - the numbers of this candidate are not to be used until they are.")
    )
    a("")
    fl = d["floquet"]
    a(
        f"Floquet cutoff: the same four states at mesh {fl['mesh']} with the cutoff {fl['factor']} times the mesh number (the converged setting of the Phase 20 R4 audit for lumped loads) against the cutoff equal to the mesh number: largest relative change "
        + ", ".join(f"{k} {f(v, 4)}" for k, v in fl["max_relative_change"].items())
        + "; cutoff converged: "
        + ("yes." if fl["converged"] else "NO.")
    )
    a("")
    a(f"### Targets and state balance (mesh {d.get('mesh_main', d['mesh_last'])})")
    a("")
    t = doc["targets_ohm"]
    L += table(
        [
            "dX_aniso definition",
            "smallest over the four states [ohm]",
            f"meets {f(t['target'], 0)} ohm",
            f"meets {f(t['stretch'], 0)} ohm (stretch)",
        ],
        [
            (
                "own mean loading (the cell's beta)",
                f(d["target"]["dx_aniso_ohm"]["min_over_states"], 1),
                "yes" if d["target"]["dx_aniso_ohm"]["meets_target"] else "NO",
                "yes" if d["target"]["dx_aniso_ohm"]["meets_stretch"] else "NO",
            ),
            (
                "design mean loading (beta/k0 of the item-2 contrast derivation)",
                f(d["target"]["dx_aniso_at_design_loading_ohm"]["min_over_states"], 1),
                "yes" if d["target"]["dx_aniso_at_design_loading_ohm"]["meets_target"] else "NO",
                "yes" if d["target"]["dx_aniso_at_design_loading_ohm"]["meets_stretch"] else "NO",
            ),
        ],
    )
    a("")
    a(
        "The same anisotropy gives a different reactance difference on a different mean loading (the surface wave is more or less bound), so both definitions are shown. The leakage the cell delivers to a beam at the design polar angle (`alpha at the beam`) is the quantity the hologram needs: compare it with the contrast requirement of the item-2 feasibility sweep at this frequency."
    )
    a("")
    L += table(
        ["quantity over the four states", "smallest", "largest", "relative spread"],
        [(k, f(v["min"], 3), f(v["max"], 3), f(v["spread_relative"], 3)) for k, v in d["states"].items()],
    )
    a("")
    a(
        f"### Diode tolerance (plus and minus twenty percent of the diode model parameters, states 0 and 1, mesh {d.get('mesh_main', d['mesh_last'])})"
    )
    a("")
    a(
        "ASSUMED_SENSITIVITY_NOT_DEVICE_DATA: the scale factors multiply the series resistance, the package inductance and the junction capacitance of the kit; the MADP package inductance is NOT_AVAILABLE_FROM_CURRENT_DATA and is zero in the nominal model, so its scale has no effect there."
    )
    a("")
    rows = [
        (
            json.dumps(x["scale"]),
            x["state"],
            f(x["delta_ms"], 3),
            f(100 * x["change_delta_ms"], 1),
            f(x["dx_aniso_ohm"], 1),
            f(100 * x["change_dx_aniso_ohm"], 1),
            f(x["alpha_beam_np_per_m"], 3),
            f(100 * x["change_alpha_beam_np_per_m"], 1),
        )
        for x in d["tolerance"]
    ]
    L += table(
        [
            "scale",
            "state",
            "delta [mS]",
            "change in delta [%]",
            "dX_aniso own [ohm]",
            "change [%]",
            "alpha at the beam [Np/m]",
            "change [%]",
        ],
        rows,
    )
    a("")
    w = d["tolerance_worst_change"]
    a(
        f"Worst relative change over all tolerance cases: delta {f(100 * w['delta_ms'], 1)} %, dX_aniso (own loading) {f(100 * w['dx_aniso_ohm'], 1)} %, dX_aniso (design loading) {f(100 * w['dx_aniso_at_design_loading_ohm'], 1)} %, alpha at the beam {f(100 * w['alpha_beam_np_per_m'], 1)} %."
    )
    a("")
    if d.get("ls_sensitivity"):
        a(
            "Package inductance of the MADP model (NOT_AVAILABLE_FROM_CURRENT_DATA, zero in the nominal model): ASSUMED values as a sensitivity, not device data."
        )
        a("")
        L += table(
            [
                "package inductance [nH]",
                "state",
                "delta [mS]",
                "change [%]",
                "dX_aniso own [ohm]",
                "change [%]",
                "alpha at the beam [Np/m]",
                "change [%]",
            ],
            [
                (
                    f(x["madp_Ls_nH"], 1),
                    x["state"],
                    f(x["delta_ms"], 3),
                    f(100 * x["change_delta_ms"], 1),
                    f(x["dx_aniso_ohm"], 1),
                    f(100 * x["change_dx_aniso_ohm"], 1),
                    f(x["alpha_beam_np_per_m"], 3),
                    f(100 * x["change_alpha_beam_np_per_m"], 1),
                )
                for x in d["ls_sensitivity"]
            ],
        )
        a("")
    a(f"### Contrast bandwidth (four states, mesh {d.get('mesh_main', d['mesh_last'])})")
    a("")
    L += table(
        [
            "frequency [GHz]",
            "smallest dX_aniso own [ohm]",
            "smallest dX_aniso design [ohm]",
            "smallest alpha at the beam [Np/m]",
        ],
        [
            (f(x["f_ghz"], 2), f(x["min_dx_aniso_ohm"], 1), f(x["min_dx_design_ohm"], 1), f(x["min_alpha_np_per_m"], 3))
            for x in d["bandwidth"]["rows"]
        ],
    )
    a("")
    for tag, label in (
        ("above_target_design_loading", "smallest design-loading dX_aniso at or above the target"),
        (
            "within_70pct_of_nominal_own_loading",
            "smallest own-loading dX_aniso within 70 percent of its 10.0 GHz value",
        ),
    ):
        b = d["bandwidth"][tag]
        if b["f_low_ghz"] is None:
            a(f"- {label}: threshold {f(b['threshold_ohm'], 1)} ohm is not met at the design frequency.")
        else:
            a(
                f"- {label} (threshold {f(b['threshold_ohm'], 1)} ohm): from {f(b['f_low_ghz'], 2)} to {f(b['f_high_ghz'], 2)} GHz, width {f(b['width_ghz'], 2)} GHz (grid step 0.25 GHz).".replace(
                    "0.25", f(0.25, 2)
                )
            )
    a("")
    a("### Control hardware (150 mm hexagon, same counting rules as the item-2 report)")
    a("")
    hw = d["hardware"]
    L += table(
        [
            "diodes per cell",
            "forward diodes per state",
            "cells",
            "diodes",
            "forward diodes",
            "forward current [mA]",
            "DC power at that current [W]",
            "DC power at 1 mA [W]",
        ],
        [
            (
                hw["diodes_per_cell"],
                hw["forward_diodes_per_state"],
                hw["n_cells"],
                hw["diodes_total"],
                hw["forward_diodes_total"],
                f(doc["i_forward_ma"], 0),
                f(hw["dc_power_w_at_i_forward"], 1),
                f(hw["dc_power_w_at_1ma"], 2),
            )
        ],
    )
    a("")
    a(
        f"DC power = forward diodes x forward current x supply rail ({f(doc['supply_rail_v'], 1)} V). Not modelled: the bias network (a connection to the hub and to every sector pad through chokes) and the diode package footprint; both are `NOT_AVAILABLE_FROM_CURRENT_DATA` here."
    )
    a("")
    return L


def ls_section(ls):
    L = ["## Sensitivity to the package inductance of the diode (candidate " + ls["spec"]["name"] + ")", ""]
    L.append(
        "ASSUMED_SENSITIVITY_NOT_DEVICE_DATA. The MADP model has no package inductance in the repo (NOT_AVAILABLE_FROM_CURRENT_DATA, zero in the nominal model); the ON state of a diode is a series resistance plus a reactance that grows with the inductance (about six ohm per tenth of a nanohenry at 10.0 GHz). The documented SMP1345-079LF model (package inductance 0.7 nH) is shown for reference. Mesh "
        + str(ls["mesh"])
        + ", states 0 and 1."
    )
    L.append("")
    L += table(
        [
            "diode model",
            "package inductance [nH]",
            "state",
            "delta [mS]",
            "dX_aniso own [ohm]",
            "dX_aniso design [ohm]",
            "alpha at the beam [Np/m]",
        ],
        [
            (
                r["kit"],
                "documented 0.7 nH" if r["kit"] == "smp" else f(r["madp_Ls_nH"], 2),
                r["state"],
                f(r["delta_ms"], 2),
                f(r["dx_aniso_ohm"], 1),
                f(r["dx_aniso_at_design_loading_ohm"], 1),
                f(r["alpha_beam_np_per_m"], 3),
            )
            for r in ls["rows"]
        ],
    )
    L.append("")
    return L


def summary(docs):
    L = ["## Summary of the candidates", ""]
    rows = []
    for doc in docs:
        d, sp = doc["derived"], doc["spec"]
        al = d["states"]["alpha_beam_np_per_m"]
        rows.append(
            (
                sp["name"],
                f(sp["period_mm"], 1),
                sp["params"].get("mode", "bowtie"),
                d["hardware"]["diodes_per_cell"],
                d["hardware"]["forward_diodes_per_state"],
                "yes" if d["convergence_all"] else "NO",
                f(d["target"]["dx_aniso_ohm"]["min_over_states"], 1),
                f(d["target"]["dx_aniso_at_design_loading_ohm"]["min_over_states"], 1),
                "yes"
                if d["target"]["dx_aniso_at_design_loading_ohm"]["meets_target"]
                or d["target"]["dx_aniso_ohm"]["meets_target"]
                else "NO",
                f(al["min"], 2) + " to " + f(al["max"], 2),
                f(d["states"]["delta_ms"]["min"], 2) + " to " + f(d["states"]["delta_ms"]["max"], 2),
                f(100 * d["tolerance_worst_change"]["alpha_beam_np_per_m"], 1),
                f(d["hardware"]["dc_power_w_at_i_forward"], 0),
                f(d["hardware"]["dc_power_w_at_1ma"], 1),
            )
        )
    L += table(
        [
            "candidate",
            "pitch [mm]",
            "mode",
            "diodes per cell",
            "forward diodes per state",
            "converged",
            "smallest dX_aniso own [ohm]",
            "smallest dX_aniso design [ohm]",
            "meets the target in either definition",
            "alpha at the beam over the states [Np/m]",
            "delta over the states [mS]",
            "worst diode-tolerance change in alpha [%]",
            "DC power at the forward current [W]",
            "DC power at 1 mA [W]",
        ],
        rows,
    )
    L += [
        "",
        "`converged` is the mesh criterion below; a candidate that is not converged is listed for the record and its numbers are not to be used. `mode`: `bowtie` = the two necks of the state's axis conduct, `notch` = every neck conducts except the two of the perpendicular axis. The 103 ohm target and 132 ohm stretch goal are the values approved from the 10.54 GHz item-2 sweep (combined 1.5 dB and 1 dB totals); at 10.0 GHz the same item-2 totals need lower anisotropies, and the leakage the surface must reach is `alpha_cap`; the 10.0 GHz values are in `FEASIBILITY_REPORT.md` and `FREQUENCY_COMPARISON.md`. `meets the target` refers to the anisotropy only; compare the leakage column with `alpha_cap` separately.",
        "",
    ]
    return L


def build(docs, search=None, ls=None):
    L = ["# Phase 28 item 3 - diode-switched CP cell candidates at 10.0 GHz (stop point)", ""]
    a = L.append
    a(
        "**Labels:** `MODEL_ONLY`, `NOT_HFSS_VALIDATED`, `PROVISIONAL_DATASHEET`, `VIA_OMITTED_APPROXIMATION`. Generated from the characterisation JSON files only."
    )
    a("")
    a(
        "Method: the normal-incidence tensor sheet of each state of the cell (`coded_cell_eval`: the Phase 26 B extraction `Y_s` from the 2 x 2 reflection matrix, the lumped diode loads of the Phase 20 R1 / R2 solver extension acting inside the MoM), "
        "its principal susceptances (mean loading and anisotropy), the TM0 phase constant of the mean sheet and the reactance difference of the two principal axes (the contrast measure of the item-2 report). "
        "This is a screening of the sheet model, not a dispersion (alpha, beta) of a row of cells; that is the state library of item 4."
    )
    a("")
    if search:
        a("## What was tried")
        a("")
        L += table(
            [
                "family",
                "pitch [mm]",
                "samples",
                "valid",
                "best smallest-state dX_aniso own [ohm]",
                "best smallest-state alpha at the beam [Np/m]",
            ],
            [
                (
                    s["family"],
                    f(s["period_mm"], 1),
                    s["n_samples"],
                    s["n_valid"],
                    f(s["best_dx_own_ohm"], 2),
                    f(s["best_alpha_np_per_m"], 3),
                )
                for s in search["families"]
            ],
        )
        a("")
        for note in search["notes"]:
            a("- " + note)
        a("")
    L += summary(docs)
    for doc in docs:
        L += candidate(doc)
    if ls:
        L += ls_section(ls)
    return "\n".join(L) + "\n"


def main(argv):
    out = Path(argv[0])
    search = None
    ls = None
    paths = []
    i = 1
    while i < len(argv):
        if argv[i] == "--search-summary":
            search = json.loads(Path(argv[i + 1]).read_text(encoding="utf-8"))
            i += 2
        elif argv[i] == "--ls-threshold":
            ls = json.loads(Path(argv[i + 1]).read_text(encoding="utf-8"))
            i += 2
        else:
            paths.append(argv[i])
            i += 1
    docs = [json.loads(Path(p).read_text(encoding="utf-8")) for p in paths]
    out.write_text(build(docs, search, ls), encoding="utf-8", newline="\n")
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1:])
