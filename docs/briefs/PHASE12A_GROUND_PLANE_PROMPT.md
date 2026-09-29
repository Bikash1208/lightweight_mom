# Phase 12A brief — PEC ground plane (grounded-substrate stack)

You are working in `lightweight_mom`. Read `CLAUDE.md`, `CONVENTIONS.md`, `GREEN_FUNCTION.md`,
`MOM_FORMULATION.md` and `S_PARAMETERS.md` before changing anything.

## Goal

Add an optional **PEC ground plane at the bottom of the substrate** (z = −h), so the solver can model
the stack used by reflectarrays and grounded metasurfaces:

```
air (z > 0, incident side) | PEC pattern at z = 0 | substrate (−h < z < 0) | PEC ground at z = −h
```

The existing ungrounded stack (air | pattern | substrate | air) must keep working **bit-identically**.
This phase adds a capability and measures what it gives; it does not build a production library,
does not re-run Phases 8–11, and does not unlock geometry inversion.

Why: the ungrounded rectangular-patch library reaches only 24.8° (S11xx) / 43.9° (S21xx) phase
coverage at 12.5 GHz, and every grounded holographic/reflectarray formulation (Fong 2010, Minatti,
Pozar 1997, Budhu & Grbic 2021) was `ARCHITECTURE_INCOMPATIBLE` for lack of a ground plane. With a
ground, S21 ≡ 0 by construction, so the S11-vs-S21 control-variable question is settled by the
architecture (reflection), not chosen by phase coverage.

## Physics (verify, don't trust)

In the TE/TM transmission-line model of `src/green_function.py`, the ground is a short circuit at
z = −h. Everything else in the pipeline (rooftops, Z assembly, Floquet sums, supercells) is
unchanged, because it only consumes `G` and the unpatterned-stack field.

- Down-looking impedance at z = 0 becomes `Z_down = j·Z_sub·tan(b·h)`, `Z_sub` = substrate modal
  impedance (TE: ωμ0/b, TM: b/(ωε0ε_c)). In the code's existing ρ/`u = exp(−2jbh)` notation this is the
  bottom-interface reflection `ρ → −1`: `Z_down = Z_sub·(1 − u)/(1 + u)`.
- `Z_eq = 1/(Y_air + 1/Z_down)`. Write it in the same overflow-free style as the existing closed
  forms (no bare `tan`, finite at a = 0, b = 0 and `b·h → 0`), e.g.
  `Z_eq = (1 − u) / (Y_air(1 − u) + Y_sub(1 + u))`. Derive both TE and TM and add the derivation to
  `GREEN_FUNCTION.md`.
- Fields: in the substrate `V(−h) = 0`; there is no field below z = −h (return zeros there, or
  reject z < −h — pick one, document it, test it).
- Unpatterned stack, top port: `Γ = (Z_in − Z_air)/(Z_in + Z_air)` with `Z_in = Z_down`; field at
  z = 0 is `1 + Γ`. The **bottom port does not exist** for a grounded stack — raise a clear error.
- Surface-wave poles: a grounded dielectric slab has a TM0 surface-wave pole with no cutoff, at
  `k0 < k_t < Re(k_sub)`. Add a diagnostic (warning + per-record flag, like the existing grating
  check) when any Floquet harmonic's `k_t` lies within a documented tolerance of that interval /
  pole. Verify numerically that for the default 4.8 mm and 8 mm periods at 10–15 GHz no harmonic
  falls there at normal incidence — report the margin, don't assume it.

## Design requirements

1. **Config.** Add the choice explicitly, e.g. a `backing` key (`air` | `pec`) in the `substrate`
   section (or a new `stack` section — choose, justify, document). Existing configs have no such key
   and must load unchanged and produce identical results: treat absence as `air`, and record this as
   the one documented default in `CONVENTIONS.md` (it preserves prior behaviour; it is not a physical
   guess). With `backing: pec`, `reference_planes.bottom_air_thickness` is meaningless — reject it
   (strict config: irrelevant keys are errors).
2. **Kernel.** Extend `LayeredGreenFunction` with an optional backing field defaulting to air, so
   every existing constructor call is untouched. Implement `equivalent_impedances`,
   `_sheet_mode_transfer`, `_stack_mode_transfer`, `sheet_fields`, `stack_fields` for the grounded
   case. Do not alter the ungrounded code path's arithmetic.
3. **Excitation / S-parameters.** `stack_tangential_field`, `ReferencePlanes`,
   `compute_s_parameters`, `scattering_from_fields`, `at_planes`, `stack_s_parameters` must handle a
   one-port (top-only) stack. S11 is the 2×2 [received, incident] block as now. S21/S12/S22 are
   **undefined** (no port 2), not zero: represent that explicitly, and make any attempt to read them
   an error rather than a silent 0.
4. **Database / downstream.** `database.py` records are 4×4 two-port matrices. Audit
   `database.py`, `sweep.py`, `runner.py`, `polarization.py`, `library_schema.py`,
   `phase8d_rev/library_entry.py` (add a stack/backing tag so grounded and ungrounded libraries are
   never mixed in `search_unit_cells`) and anything else touching port 2 or `bottom_air_thickness`
   (grep: `hfss_comparison.py`, `validation_cases.py`, `validation_export.py`, `phase8b/conventions.py`,
   `phase8c/topology_audit.py`, `examples/sweep_config.yaml`). Choose the minimal, explicit one-port
   representation; bump schema versions where the on-disk format changes; old files must still load.
   Frozen earlier-phase modules: extend, don't rewrite. If a module can't support one-port records
   cleanly, make it reject grounded input with a clear error and list it in the report.
5. **CAD.** `cad_export.py`: add the ground as a metal solid at the substrate bottom when
   `backing: pec` (needed later for HFSS). Optional-dependency rules unchanged.

## Validation (all as tests; property-based, not "it runs")

- Existing suite: all current tests pass with no edits to their expectations. Prove the ungrounded
  path is bit-identical (e.g. re-solve a few Phase 7 cells and compare to stored records exactly),
  and that nothing under `results/` changed.
- Analytic grounded slab, no pattern: `|Γ| = 1` for tanδ = 0; phase matches
  `Γ = (j·η_s·tan(k_s·h) − η0)/(j·η_s·tan(k_s·h) + η0)`, `η_s = η0/√ε_r`, over a frequency sweep; lossy:
  `|Γ| < 1` and `1 − |Γ|²` equals absorbed power; `h → 0` gives `Γ → −1`.
- Independent reference: extend `tests/reference_layered.py` (independent field matching, not the
  kernel's code) with a PEC termination and check `sheet_fields` / `stack_fields` at many z,
  including `E_t(z = −h) = 0` and the `H_t` jump `ẑ×ΔH = J` at the sheet.
- Kernel: new closed form equals the textbook `tan` form where the latter is finite; stable at
  a = 0, b = 0, very thick/thin substrates; `Z_eq` Hermitian-part sign consistent with passivity.
- Patterned: full-PEC cell gives S11 = −1 (ground irrelevant); lossless patch gives
  `|S11_xx|² + |S11_yx|² = 1` to solver precision; reciprocity `S11 = S11ᵀ` at normal incidence;
  x/y and L↔W mirror symmetry for rectangles; mesh/Floquet convergence study for one grounded patch.
- Config: `backing: pec` with `bottom_air_thickness` rejected; unknown backing values rejected;
  bottom-port and S21 access on grounded stacks raise.
- Physics sanity (report, don't assert a target): a small grounded rectangular-PEC-patch sweep on the
  **same** Rogers RT/Duroid 5880, 1.575 mm, 4.8 mm period, at 12.5 GHz and a few frequencies in
  10–15 GHz (e.g. L = W over the Phase 7 range, ~10–20 values). Measure reflection-phase coverage,
  largest gap and |S11| (loss), side by side with the ungrounded Phase 7 numbers, and compare the
  curve's shape qualitatively with Pozar et al. 1997's variable-size-patch phase curves (cite; don't
  digitize unless you record the source figure and tolerance). If coverage is still far from 360°,
  say so plainly — thicker/thinner substrate or other topologies are the next question, not this
  phase's to answer.

## Deliverables

- Code + tests as above; `examples/run_phase12a_ground_plane.py` writing
  `results/phase12a_ground_plane/` (sweep CSV/JSON, plots: phase vs L at several frequencies,
  |S11| vs L, grounded vs ungrounded comparison, convergence), deterministic outputs.
- `PHASE12A_GROUND_PLANE.md` (method, derivation pointer, measured numbers, limitations) and updates
  to `CONVENTIONS.md`, `GREEN_FUNCTION.md`, `S_PARAMETERS.md`, `DATABASE.md`, `README.md`
  (layout table, status, "Still open"), `CLAUDE.md` if commands change.
- Explicit final status from a fixed set, e.g. `GROUNDED_STACK_VALIDATED` /
  `GROUNDED_STACK_PARTIAL` (with the reason), plus the measured grounded phase coverage. HFSS status
  stays whatever it is; geometry inversion and STEP export stay LOCKED.

## Rules

- No fabricated data; measured numbers only; label anything synthetic.
- Don't change the ungrounded physics, earlier phases' results, or their gate outcomes.
- Don't pick a control variable by phase coverage; here reflection is fixed by the architecture.
- If a design choice above is genuinely ambiguous (config shape, one-port record format), make the
  minimal explicit choice, document why, and flag it in the final report rather than stopping.
- Run the full test suite at the end and report counts.
