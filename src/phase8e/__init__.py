"""Phase 8E: the physical control-variable / holographic-scattering bridge.

Builds strictly on Phase 8B (`phase8b.object_wave`/`reference_wave`/`holographic_modulation`/
`control_variable`), Phase 8C (`phase8c.scattering_analysis`/`architecture_test`/`topology_audit`),
and Phase 8D-REV (`phase8d_rev.library_entry`/`response_matching`/`build_request`) -- nothing here
re-derives a conclusion those modules already reached; it extends them with the pieces the project
never had: reference/object-wave spectral analysis, reference-plane propagation, a generalized
(vector/matrix) control-variable candidate, a genuine 2-D-complex-plane controllability test, and
continuous-vs-discrete beam reconstruction *before* any geometry inversion.

Frozen, never modified or re-run: the MoM solver, mesh classifier, Green's function, Phase 7
database/geometry/STEP files, the Phase 8D-REV `results/universal_database/` outputs. No GUI, no
final MTS STEP, no geometry inversion (Sec. 31) -- this package only ever decides whether geometry
inversion is *permitted* (`decision_tree.py`'s Sec. 20 gate).

Central modules:

- `provenance.py` -- Sec. 27's formulation-provenance registry (candidates A-D).
- `reference_field_input.py` -- Sec. 3: real-or-synthetic reference field, explicit about which.
- `spectral.py` -- Sec. 7-8: 2-D FFT spectral analysis of the reference field and the object wave.
- `reference_plane.py` -- Sec. 6: reference-plane compatibility and angular-spectrum propagation
  (a real derivation, `periodic_bc.normal_wavenumber` reused, never an arbitrary phase offset).
- `control_response.py` -- Sec. 9-10: C_req(x,y) for each candidate formulation.
- `architecture_compatibility.py` -- Sec. 11-12: does the real (ungrounded) topology satisfy each
  candidate's assumptions.
- `complex_controllability.py` -- Sec. 13-14: 2-D-complex-plane library matching (never phase
  coverage alone).
- `local_periodicity.py` -- Sec. 15: always `LOCAL_PERIODICITY_UNVALIDATED` (no coupled-cell
  solver exists to validate it), with the supporting sampling-density number reported anyway.
- `candidate_library_period.py` -- Sec. 16-17: reuses `phase8d_rev.candidate_period`/`build_request`.
- `beam_reconstruction.py` -- Sec. 18-19: continuous and discrete aperture -> far-field comparison.
- `decision_tree.py` -- Sec. 20, 22, 25, 29: the fixed-order gate and 13-value status taxonomy.
"""
