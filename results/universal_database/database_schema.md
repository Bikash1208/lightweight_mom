# database_schema.md

Generated from the actual `phase8d_rev` dataclasses (`dataclasses.fields`), so this file cannot drift from the code.

## UnitCellRecord (Sec. 7)

- `record_id`: `<class 'str'>`
- `library_id`: `<class 'str'>`
- `topology_id`: `<class 'str'>`
- `simulation_id`: `<class 'str'>`
- `topology`: `<class 'str'>`
- `p_x_mm`: `<class 'float'>`
- `p_y_mm`: `<class 'float'>`
- `L_mm`: `<class 'float'>`
- `W_mm`: `<class 'float'>`
- `extra_geometry_parameters`: `<class 'dict'>`
- `substrate_name`: `<class 'str'>`
- `epsilon_r`: `<class 'float'>`
- `tan_delta`: `<class 'float'>`
- `substrate_thickness_mm`: `<class 'float'>`
- `additional_layers`: `<class 'tuple'>`
- `frequency_hz`: `<class 'float'>`
- `incidence_theta_deg`: `<class 'float'>`
- `incidence_phi_deg`: `<class 'float'>`
- `polarization_basis`: `<class 'str'>`
- `propagation_direction`: `<class 'str'>`
- `periodic_x`: `<class 'bool'>`
- `periodic_y`: `<class 'bool'>`
- `port_type`: `<class 'str'>`
- `floquet_mode_definition`: `<class 'str'>`
- `s_parameters`: `<class 'dict'>`
- `reflected_power_co`: `<class 'float'>`
- `reflected_power_cross`: `<class 'float'>`
- `transmitted_power_co`: `<class 'float'>`
- `transmitted_power_cross`: `<class 'float'>`
- `absorbed_power_fraction`: `<class 'float'>`
- `polarization_conversion_note`: `<class 'str'>`
- `phase_coverage_metadata`: `None`
- `provenance`: `<class 'src.phase8d_rev.validation_levels.Provenance'>`

## UnitCellFamily (Sec. 9)

- `family_id`: `<class 'str'>`
- `topology`: `<class 'str'>`
- `parameter_schema`: `<class 'tuple'>`
- `fixed_parameters`: `<class 'dict'>`
- `material_stacks`: `<class 'tuple'>`
- `supported_periods_mm`: `<class 'tuple'>`
- `supported_frequency_range_hz`: `<class 'tuple'>`
- `supported_polarizations`: `<class 'tuple'>`
- `supported_incidence`: `<class 'tuple'>`
- `solver_method`: `<class 'str'>`
- `validity_constraints`: `<class 'dict'>`
- `geometry_export_method`: `<class 'str'>`
- `library_ids`: `<class 'tuple'>`

## Validation levels (Sec. 23)

`DATA_PRESENT`, `NUMERICALLY_VALIDATED`, `CROSS_SOLVER_VALIDATED`, `HFSS_VALIDATED`, `EXPERIMENTALLY_VALIDATED`

## LibraryEntry roles (Sec. 3, 17)

`reference_validation_library`, `project_library`

## Period/frequency/substrate/polarization/incidence are per-entry, not global

Period (p_x_mm, p_y_mm), frequency, substrate, polarization and incidence are all per-LibraryEntry fields, never global constants; a family aggregates them across every registered entry of one topology (see family.py).

