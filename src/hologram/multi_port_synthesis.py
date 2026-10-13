"""Phase 15 Part C (gap-fill): runs N independent single-aperture HMSA syntheses
(`hologram.hmsa_synthesis.synthesize_formulation_a/b`, UNMODIFIED), one per `ports.PortRegion` --
see `ports.py`'s module docstring for why a port needs no new synthesis machinery of its own, only
its own cropped native grid, feed and beam list.

**Scope reduction, documented (not hidden)**: every port of one multi-port design shares the SAME
formulation and the same dispersion table / phase library (`per_port_synthesize`, below, is one
callable bound once by the caller) -- the brief's own "each [port] with ... its own method" is not
implemented as a genuinely independent per-port formulation choice in this first version. Each
port DOES get its own feed position, own beam list, and its own polarization compensation (every
other part of Part A/B is fully per-port)."""

from dataclasses import dataclass, replace

import numpy as np

from .ports import IsolationFlag, PortRegion, isolation_estimate


def region_cell_mask(region: PortRegion, x_m, y_m, earlier_regions=()):
    """Phase 23: which cells of a port's design grid belong to the port -- cell centre inside the
    region's polygon and not inside an EARLIER port's region (a centre exactly on a shared edge is
    owned by the lower-numbered port only, so no cell is built twice). Every cell for a rectangle."""
    mask = region.contains_array(x_m, y_m)
    for other in earlier_regions:
        mask &= ~other.contains_array(x_m, y_m)
    return mask


def mask_design_to_region(design, mask):
    """A copy of `design` (`HmsaDesign` or `CodedDesign`) with every cell outside `mask` removed: its
    `aperture_field` is zero there (nothing radiates from outside the sub-aperture) and it is dropped
    from `cells` (so exports place no patch there -- the empty state). `level_map`/`x_m`/`y_m` keep the
    full cropping-box grid."""
    field = None
    if design.aperture_field is not None:
        field = np.where(mask, design.aperture_field, 0.0)
    kept = tuple(c for c in design.cells if bool(mask[c.i, c.j]))
    return replace(design, cells=kept, aperture_field=field)


def port_native_grid(region: PortRegion, resolution_m):
    """A grid spanning `region` at `resolution_m` pitch (finer than the eventual unit-cell period
    -- `synthesize_formulation_a/b`'s own `resample_to_cell_grid` resamples this down to the
    chosen period, exactly as it already does for a real imported `.fld`)."""
    if resolution_m <= 0.0:
        raise ValueError(f"resolution_m must be > 0, got {resolution_m!r}")
    nx = max(4, round(region.width_m / resolution_m))
    ny = max(4, round(region.height_m / resolution_m))
    x = np.linspace(region.x_min_m, region.x_max_m, nx)
    y = np.linspace(region.y_min_m, region.y_max_m, ny)
    return np.meshgrid(x, y, indexing="xy")


def port_reference_field(region: PortRegion, resolution_m, k_swg=None, uploaded_field=None):
    """This port's own reference field, on its own native grid: either the standard analytic
    monopole about `region`'s own feed, or an uploaded field resampled onto this port's grid --
    reusing `hologram.multi_feed.superposed_reference_field` on a single-feed list (a port's own
    field IS exactly that function's "one feed" case; no new interpolation code is written here)."""
    from .multi_feed import FeedSpec, superposed_reference_field

    x_m, y_m = port_native_grid(region, resolution_m)
    feed = FeedSpec(region.feed_x_m, region.feed_y_m, amplitude=1.0, reference_field=uploaded_field)
    ref = superposed_reference_field([feed], x_m, y_m, k_swg)
    if uploaded_field is None:
        # Phase-only, as an uploaded .fld is: the monopole's rho^-1/2 spreading would be read as an
        # amplitude taper that Formulation A's modulation range cannot follow (measured on the
        # full board: 54% unreachable analytic vs 0% phase-only, uploaded field 15%).
        phase_only = np.exp(1j * np.angle(ref.Ex))
        ref = replace(ref, Ex=phase_only, Ey=np.zeros_like(phase_only), phase_deg=np.degrees(np.angle(phase_only)))
    return ref


@dataclass(frozen=True)
class PortSynthesisResult:
    port_id: int
    region: PortRegion
    design: object  # hologram.hmsa_synthesis.HmsaDesign
    feed_xy_m: tuple


@dataclass(frozen=True)
class MultiPortDesign:
    ports: "tuple[PortSynthesisResult, ...]"
    isolation_flags: "tuple[IsolationFlag, ...]"
    period_m: float


def synthesize_ports(
    port_regions,
    port_beams,
    per_port_synthesize,
    period_m,
    resolution_m=None,
    k_swg=None,
    uploaded_fields=None,
    min_lambda_sw=1.0,
):
    """`port_beams`: `{port_id: [Beam, ...]}`. `per_port_synthesize(reference_field, beams,
    feed_xy_m) -> HmsaDesign`: a caller-supplied partial application of `synthesize_formulation_a`
    or `_b` with every formulation-specific argument (dispersion table / library candidates /
    X_avg / modulation index / etc.) already bound -- see module docstring for why this keeps the
    SAME formulation across every port. `uploaded_fields`: optional `{port_id: ReferenceFieldGrid}`.
    `k_swg`: needed both for an analytic (non-uploaded) port feed AND for the isolation estimate's
    own lambda_sw; `None` skips the isolation estimate entirely (returns `()`) rather than guess."""
    if not port_regions:
        raise ValueError("synthesize_ports: at least one port region is required")
    resolution_m = resolution_m if resolution_m is not None else period_m / 4.0
    uploaded_fields = uploaded_fields or {}
    results = []
    for region in port_regions:
        beams = port_beams.get(region.port_id)
        if not beams:
            raise ValueError(f"port {region.port_id}: at least one beam is required")
        ref_field = port_reference_field(region, resolution_m, k_swg, uploaded_fields.get(region.port_id))
        feed_xy_m = (region.feed_x_m, region.feed_y_m)
        if region.is_polygon:
            # Phase 23: the synthesis callable may take `region=` to see only the cells inside the
            # polygon while it measures/refines; the result is then masked here regardless.
            design = per_port_synthesize(ref_field, beams, feed_xy_m, region=region)
            earlier = [r for r in port_regions if r is not region and r.port_id < region.port_id]
            mask = region_cell_mask(region, design.x_m, design.y_m, earlier)
            design = mask_design_to_region(design, mask)
        else:
            design = per_port_synthesize(ref_field, beams, feed_xy_m)
        results.append(PortSynthesisResult(region.port_id, region, design, feed_xy_m))
    isolation_flags = isolation_estimate(port_regions, k_swg, min_lambda_sw) if k_swg else ()
    return MultiPortDesign(tuple(results), isolation_flags, period_m)


def combined_board_lw_grid_m(multi_port: MultiPortDesign, cell_m_fn):
    """Assembles ONE rectangular `(length_m, width_m)` nested-list grid spanning the union of
    every port's own rectangle, on ONE shared `period_m` lattice indexed from the board's own
    corner (the convention `export.array_export.write_hfss_step` assumes) -- for EXPORT/
    VISUALIZATION only; each port's own synthesis result (`MultiPortDesign.ports`) is reported to
    the user separately and is NOT re-gridded.

    `cell_m_fn(GeometryCell) -> (length_m, width_m) | None` converts one solved cell to its
    export geometry (`None` for an unreachable cell -- rendered as "no patch").

    Each combined-grid cell's geometry is the NEAREST solved cell (by centre distance) of
    whichever port's own rectangle contains that grid cell's centre: a port's own grid is centred
    independently within its own rectangle (`ports.py`'s own docstring), so two neighbouring
    ports do not necessarily share one exact lattice phase at the seam between them -- nearest-cell
    snapping absorbs that sub-period mismatch, a documented, bounded (<= one period) limitation of
    the COMBINED DRAWING, not a claim that the physical layout is seamless at every port boundary.
    A grid cell inside no port's rectangle (a `gap_m` guard row) is `(0.0, 0.0)` ("no patch").

    Returns a `CombinedBoard`.
    """
    period_m = multi_port.period_m
    x_min = min(p.region.x_min_m for p in multi_port.ports)
    x_max = max(p.region.x_max_m for p in multi_port.ports)
    y_min = min(p.region.y_min_m for p in multi_port.ports)
    y_max = max(p.region.y_max_m for p in multi_port.ports)
    n_x = max(1, round((x_max - x_min) / period_m))
    n_y = max(1, round((y_max - y_min) / period_m))
    px, py = (x_max - x_min) / n_x, (y_max - y_min) / n_y

    grid = []
    for i in range(n_y):
        cy = y_min + (i + 0.5) * py
        row = []
        for j in range(n_x):
            cx = x_min + (j + 0.5) * px
            owner = next((p for p in multi_port.ports if p.region.contains(cx, cy)), None)
            cell_value = None
            if owner is not None and owner.design.cells:
                cells = owner.design.cells
                idx = int(np.argmin([(c.x_m - cx) ** 2 + (c.y_m - cy) ** 2 for c in cells]))
                cell_value = cell_m_fn(cells[idx])
            row.append(cell_value if cell_value is not None else (0.0, 0.0))
        grid.append(row)
    return CombinedBoard(
        lw_grid_m=grid,
        width_m=(x_max - x_min),
        height_m=(y_max - y_min),
        x_min_m=x_min,
        y_min_m=y_min,
        period_x_m=px,
        period_y_m=py,
    )


@dataclass(frozen=True)
class CombinedBoard:
    lw_grid_m: list
    width_m: float
    height_m: float
    x_min_m: float
    y_min_m: float
    period_x_m: float  # the ACTUAL per-cell pitch used (board width / n_x -- close to but not
    period_y_m: float  # necessarily exactly `period_m`, since n_x/n_y are rounded integers)

    def cell_center_m(self, i, j):
        """The (x_m, y_m) centre this module's own corner-indexed lattice places grid cell
        `(i, j)` at -- the SAME convention `export.array_export.write_hfss_step`/
        `app.plots.patch_layout_png` need to draw each cell at its true board position."""
        cx = self.x_min_m + (j + 0.5) * self.period_x_m
        cy = self.y_min_m + (i + 0.5) * self.period_y_m
        return cx, cy


def export_frame_offset_m(multi_port: MultiPortDesign, board: CombinedBoard):
    """The point `(x, y)` in the board's global frame that the STEP origin sits on. `export.array_export`
    writes a board centred on the origin (cell j at -L/2 + (j + 1/2) P, L = n_x P), whereas the combined
    grid and the ports are indexed from the board's own corner. Feeds and port rectangles must be
    shifted by this offset before they are written next to the cells."""
    n_x = len(board.lw_grid_m[0])
    n_y = len(board.lw_grid_m)
    period = multi_port.period_m
    return board.x_min_m + n_x * period / 2.0, board.y_min_m + n_y * period / 2.0


def export_frame_feeds_and_rects_mm(multi_port: MultiPortDesign, board: CombinedBoard):
    """`(feeds_mm, rects_mm)` in the STEP frame, one entry per port in port order: each feed as
    `(x_mm, y_mm)` and each port rectangle as `(x_min_mm, x_max_mm, y_min_mm, y_max_mm)`."""
    ox, oy = export_frame_offset_m(multi_port, board)
    feeds = [((p.feed_xy_m[0] - ox) * 1e3, (p.feed_xy_m[1] - oy) * 1e3) for p in multi_port.ports]
    rects = [
        (
            (p.region.x_min_m - ox) * 1e3,
            (p.region.x_max_m - ox) * 1e3,
            (p.region.y_min_m - oy) * 1e3,
            (p.region.y_max_m - oy) * 1e3,
        )
        for p in multi_port.ports
    ]
    return feeds, rects
