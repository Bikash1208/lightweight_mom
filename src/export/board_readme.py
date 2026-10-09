"""Phase 21 D3: the README that ships with a multi-port board export -- the port/feed table, the beam
table, the cell-pair phases and how to excite one port or all of them in HFSS."""

from .exact_board import fmt


def _beam_rows(ports):
    rows = []
    for p in ports:
        report = p.get("report")
        if not report:
            continue
        for k, b in enumerate(report["beams"]):
            cross = "-" if b.get("cross_pol_level_db") is None else "%.1f" % b["cross_pol_level_db"]
            rows.append(
                "| %d | %d | (%.1f°, %.1f°) | %.1f°, %.1f° | %.2f | %.1f | %s |"
                % (
                    p["port_id"] + 1,
                    k + 1,
                    b["requested_theta_deg"],
                    b["requested_phi_deg"],
                    b["theta_deg"],
                    b["phi_deg"],
                    b["pointing_error_deg"],
                    b["directivity_dbi"],
                    cross,
                )
            )
    return rows


def build_board_readme(board, result, name, eps_r, tan_delta, frequency_hz, files):
    """The README text. `result` is the `/synthesize_ports` job result; `files` the exported file names."""
    ports = result["ports"]
    f_ghz = frequency_hz / 1e9
    lines = [
        "# %s" % name,
        "",
        "%d-port multi-port HMSA board, %.4g GHz. Substrate eps_r %g, tan d %g, %s mm thick, PEC ground. "
        "Cells %s mm. The board is %s x %s mm (x %s..%s, y %s..%s), one monopole per sub-aperture."
        % (
            len(ports),
            f_ghz,
            eps_r,
            tan_delta,
            fmt(board.thickness_mm),
            fmt(board.period_mm),
            fmt(board.width_mm),
            fmt(board.height_mm),
            fmt(board.x_min),
            fmt(board.x_max),
            fmt(board.y_min),
            fmt(board.y_max),
        ),
        "%d patches in total. Every coordinate is exact (rounded to 0.001 mm; the largest rounding applied was "
        "%.4g mm), and each port's cells sit at that port's own positions, not on a shared lattice." % (
            len(board.patches),
            board.max_rounding_mm,
        ),
        "",
        "## Ports and feeds",
        "",
        "| Port | Feed (x, y) mm | Sub-aperture x mm | Sub-aperture y mm | Cells | Patches |",
        "|---|---|---|---|---|---|",
    ]
    for p, (fx, fy), (x0, x1, y0, y1) in zip(ports, board.feeds_mm, board.port_rects_mm, strict=True):
        lines.append(
            "| %d | (%s, %s) | %s .. %s | %s .. %s | %d | %d |"
            % (p["port_id"] + 1, fmt(fx), fmt(fy), fmt(x0), fmt(x1), fmt(y0), fmt(y1), p["n_cells"], p["patch_count"])
        )
    beam_rows = _beam_rows(ports)
    lines += ["", "## Beams (model prediction, not an HFSS result)", ""]
    if beam_rows:
        lines += [
            "| Port | Beam | Target (θ, φ) | Predicted (θ, φ) | Pointing error ° | Directivity dBi | Cross-pol dB |",
            "|---|---|---|---|---|---|---|",
            *beam_rows,
        ]
    else:
        lines.append("NOT_AVAILABLE_FROM_CURRENT_DATA (the design was run without a loss tangent, so no pattern was computed).")
    two_state = result.get("two_state")
    lines += ["", "## Cell states", ""]
    if two_state and two_state.get("states"):
        lines += ["| State | Cell | abs(Γ) dB | Phase at the patch plane |", "|---|---|---|---|"]
        for k, s in enumerate(two_state["states"]):
            lines.append("| %d | %s | %.3f | %.1f° |" % (k, s["label"], s["magnitude_db"], s["phase_deg_patch_plane"]))
        lines.append("")
        lines.append(
            "Separation %.1f°. Check the pair with one HFSS unit cell first (master/slave boundaries, %s mm period, Floquet "
            "port, PEC ground, %.4g GHz) and nudge the patch size if the separation is far from 180°."
            % (two_state["achieved_separation_deg"], fmt(board.period_mm), f_ghz)
        )
    elif result.get("state_table"):
        lines += ["| Level | Target phase | Cells | Mean achieved | Mean error |", "|---|---|---|---|---|"]
        for r in result["state_table"]:
            ach = "-" if r["mean_achieved_phase_deg"] is None else "%.1f°" % r["mean_achieved_phase_deg"]
            err = "-" if r["mean_phase_error_deg"] is None else "%.1f°" % r["mean_phase_error_deg"]
            lines.append("| %d | %.0f° | %d | %s | %s |" % (r["level"], r["target_phase_deg"], r["n_cells"], ach, err))
    else:
        lines.append("Continuous-phase design: each cell has its own size (see layout.csv).")
    lines += [
        "",
        "## Files",
        "",
        *["- `%s`" % f for f in files],
        "",
        "`hfss_build_board.py` is a native AEDT script (HFSS: Tools > Run Script). It builds the substrate, ground "
        "with a coax opening at each feed, all patches united into one body \"Patches\" (Perfect E), the PEC "
        "monopoles with %d lumped ports (50 ohm), an air box with a radiation boundary, a %.4g GHz setup, a sweep "
        "and a far-field sphere. Design variables: `mono_h` (default 7.1 mm: set it to your monopole), `r_pin`, `r_out`, "
        "`pad_xy`, `pad_top`, `pad_bot`. Set `BUILD_MONOPOLES = False` at the top for geometry only. Creating a few "
        "thousand rectangles takes minutes." % (len(ports), f_ghz),
        "",
        "## Exciting the ports in HFSS",
        "",
        "- **One port at a time**: HFSS > Excitations > Edit Sources, set port N to 1 W and every other port to 0 W "
        "(or simply read the S-parameters: S(N,N) is the match, S(M,N) the isolation to port M).",
        "- **All ports together**: Edit Sources, give every port 1 W at the phase offsets you want (all 0° is the "
        "default); the far-field sphere then shows the superposed beams.",
        "- The board is large (about %.0f wavelengths across): use one port at a time first, an IE/FE-BI region or a "
        "mixed-order mesh, and compare with the model's beam table above."
        % (max(float(board.width_mm), float(board.height_mm)) / (299792458.0 / frequency_hz * 1e3)),
        "",
        "Model-only: none of these numbers has been compared with an HFSS solve of this board.",
        "",
    ]
    return "\n".join(lines)
