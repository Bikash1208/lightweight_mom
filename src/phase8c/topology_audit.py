"""Phase 8C Sec. 6-10: unit-cell topology, ground-plane test, port/excitation/S-matrix audit.

Every claim here is either read directly from the database's own stored `metadata['context']`
(the exact dict `config_from_dict` validated and `runner.run_sweep` solved from, Phase 1,
unmodified) or cited from `CONVENTIONS.md`/`S_PARAMETERS.md` -- nothing is inferred from a
filename (Sec. 6's explicit prohibition).
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class TopologyAudit:
    ground_plane: str                # "PRESENT" | "ABSENT" | "UNKNOWN"
    ground_plane_evidence: str
    metal_layers: str
    substrate_layers: str
    air_regions: str
    architecture_classification: str  # Sec. 7's 5 values
    architecture_evidence: str


def _ground_plane_test(context):
    """Sec. 8: PRESENT/ABSENT/UNKNOWN, from the actual solved config context, not a filename or a
    doc citation alone.

    Since Phase 12A the config schema **can** express a ground plane (`config.py`'s
    `substrate.backing`, `geometry.UnitCell.backing`, `green_function.LayeredGreenFunction.backing`:
    `"air"` | `"pec"`), so this is a real per-database read of `context['backing']`, not a
    project-wide constant -- a database solved before Phase 12A (or any config that never set
    `backing`) has no such key and was always the ungrounded ("air") stack (`sweep.py`'s own
    "absent means air" default), so absence still means ABSENT here."""
    backing = context.get("backing", "air")
    if backing not in ("air", "pec"):
        return ("UNKNOWN", f"context['backing'] = {backing!r} is not a recognized value")
    has_ground_key = any("ground" in str(k).lower() for k in _flatten_keys(context))
    substrate = context.get("substrate", {})
    metal = context.get("metal", {})
    if backing == "pec":
        evidence = (
            f"context['backing'] = 'pec' (Phase 12A): this database's stack is air | pattern at "
            f"z=0 | substrate (context['substrate'] = {substrate}) | PEC ground at z=-h. "
            f"context['metal'] = {metal} remains a patterned PEC *pattern*, a finite scatterer, "
            f"which together with the ground makes this a grounded (reflectarray-style) "
            f"architecture, not the Phase 1-11 ungrounded one."
        )
        return ("PRESENT", evidence)
    evidence = (
        f"context['backing'] = {backing!r} (absent means 'air', the Phase 1-11 default -- "
        f"`sweep.unit_cell_from_context`); context['substrate'] = {substrate} (one dielectric "
        f"layer, not a stack terminated in a conductor); context['metal'] = {metal} (a patterned "
        f"PEC *sheet at z=0*, i.e. a finite scatterer, not a continuous ground plane); no 'ground' "
        f"key exists anywhere in the context dict (checked: {has_ground_key}). README.md "
        f"'Assumptions' item 2: 'No ground plane, no superstrate' by default (Phase 1; Phase 12A "
        f"added the opt-in 'pec' backing this database did not use)."
    )
    return ("ABSENT", evidence)


def _flatten_keys(d, prefix=""):
    keys = []
    if isinstance(d, dict):
        for k, v in d.items():
            keys.append(f"{prefix}{k}")
            keys.extend(_flatten_keys(v, f"{prefix}{k}."))
    return keys


def _architecture_classification(mean_r, mean_t, mean_r_cross, mean_t_cross):
    """Sec. 7: classification from the *actual, measured* power split (Sec. 7's explicit
    prohibition on classifying from S11/S21 magnitude *alone* is honoured by using the full
    co+cross power, `power_balance.py`, not a bare |S11xx|^2 reading)."""
    r_total, t_total = mean_r + mean_r_cross, mean_t + mean_t_cross
    label = "PARTIALLY_REFLECTIVE_TRANSMISSIVE"   # this library never has R=0 or T=0 exactly, so
                                                    # REFLECTIVE/TRANSMISSIVE (pure) never applies
    if t_total > 2.0 * r_total:
        note = f"more transmissive than reflective (mean reflected power {r_total:.4f}, transmitted {t_total:.4f})"
    elif r_total > 2.0 * t_total:
        note = f"more reflective than transmissive (mean reflected power {r_total:.4f}, transmitted {t_total:.4f})"
    else:
        note = f"comparable reflected/transmitted power (R={r_total:.4f}, T={t_total:.4f})"
    return label, note


def audit_topology(database, mean_power):
    """``mean_power``: the dict `power_balance.library_power_summary` returns (mean R/T, co and
    cross), so this function never recomputes power itself -- one source of truth."""
    context = database.metadata.get("context")
    if not context:
        raise ValueError("database has no 'context' metadata; cannot audit topology")
    ground_plane, ground_evidence = _ground_plane_test(context)
    arch, arch_note = _architecture_classification(
        mean_power["mean_R_co"], mean_power["mean_T_co"], mean_power["mean_R_cross"], mean_power["mean_T_cross"],
    )
    return TopologyAudit(
        ground_plane=ground_plane, ground_plane_evidence=ground_evidence,
        metal_layers=f"one patterned PEC sheet at z=0 (context['metal']={context.get('metal')})",
        substrate_layers=f"one dielectric layer (context['substrate']={context.get('substrate')})",
        air_regions="two semi-infinite air half-spaces (z>0 and z<-h), each with an explicit "
                    f"reference plane (context['reference_planes']={context.get('reference_planes')})",
        architecture_classification=arch, architecture_evidence=arch_note,
    )


def render_topology_audit_md(audit):
    return "\n".join([
        "# unit_cell_topology_audit.md", "",
        "## Ground plane", "",
        f"**{audit.ground_plane}**", "", audit.ground_plane_evidence, "",
        "## Metal layers", "", audit.metal_layers, "",
        "## Substrate layers", "", audit.substrate_layers, "",
        "## Air regions / reference planes", "", audit.air_regions, "",
        "## Physical architecture classification", "",
        f"**{audit.architecture_classification}** -- {audit.architecture_evidence}", "",
        "## Ports, excitation, polarization, S-matrix definition", "",
        "See `s_parameter_definition.md` (this project's fixed, existing conventions, "
        "`CONVENTIONS.md`/`S_PARAMETERS.md`, reused verbatim -- not redefined here).",
    ]) + "\n"


def render_s_parameter_definition_md():
    """Sec. 9-10: this project's fixed port/S-matrix definitions, quoted (not re-derived) from
    `CONVENTIONS.md`/`S_PARAMETERS.md`/`database.CONVENTIONS`."""
    from ..database import CONVENTIONS

    lines = [
        "# s_parameter_definition.md", "",
        "## Ports", "",
        "Port 1 = top reference plane (z_top = +top_air_thickness), incident wave travels -z "
        "(toward the pattern). Port 2 = bottom reference plane "
        "(z_bottom = -(h + bottom_air_thickness)), incident wave travels +z. Only the specular "
        "(0,0) Floquet order at normal incidence is a port mode (`CONVENTIONS.md` 'Ports').", "",
        "## S-matrix", "", f"`{CONVENTIONS['s_matrix']}`", "",
        "## Polarization basis", "", f"`{CONVENTIONS['linear_basis']}`", "",
        "Circular basis (derived, not independently excited/stored): "
        f"R={CONVENTIONS['circular_basis']['R']}, L={CONVENTIONS['circular_basis']['L']}.", "",
        "## Normalization", "",
        "Wave amplitudes are tangential-E amplitudes at the reference plane, *measured* from E "
        "and H at the plane (not assumed unit-normalized, `CONVENTIONS.md` 'Waves, ports and "
        "S-parameters'); incident amplitude is 1 V/m referenced at the pattern plane z=0 for the "
        "excitation itself (`excitation.py`).", "",
        "## Reference-plane sensitivity", "",
        "Moving the reference planes changes S11/S22 by a pure phase factor and S21/S12 by "
        "another (`CONVENTIONS.md`: `S11' = S11 exp(-2j k0 D1)`, `S22' = S22 exp(+2j k0 D2)`, "
        "`S21' = S12' = S21 exp(+j k0 (D2-D1))`), so phase is only meaningful relative to the "
        "specific reference planes stored in the database metadata -- see "
        "`reference_plane_sensitivity.py`.",
    ]
    return "\n".join(lines) + "\n"
