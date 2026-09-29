# query_spec.md

`query.search_unit_cells(universal_db, query)` applies the following filters in this fixed order (Sec. 14); it never ranks surviving candidates by phase coverage or any other controllability metric.

## UnitCellQuery fields

- `topology`: `str | None` (`None` = unconstrained)
- `substrate_name`: `str | None` (`None` = unconstrained)
- `eps_r`: `float | None` (`None` = unconstrained)
- `tan_delta`: `float | None` (`None` = unconstrained)
- `substrate_thickness_mm`: `float | None` (`None` = unconstrained)
- `p_x_mm`: `float | None` (`None` = unconstrained)
- `p_y_mm`: `float | None` (`None` = unconstrained)
- `frequency_hz`: `float | None` (`None` = unconstrained)
- `polarization`: `str | None` (`None` = unconstrained)
- `incidence`: `str | None` (`None` = unconstrained)
- `min_validation_status`: `str | None` (`None` = unconstrained)

## Filter order

1. topology
2. material_stack (substrate name/eps_r/tan_delta/thickness)
3. frequency (exact stored grid point only, no interpolation)
4. period (p_x_mm, p_y_mm)
5. polarization
6. incidence
7. response_availability (at least one OK-solved cell)
8. validity_status (>= min_validation_status, if given)

## Result

`QueryResult.status` is `"MATCH"` (with the surviving `LibraryEntry` objects) or `"DATABASE_COVERAGE_MISSING"` (with `eliminated_at` naming the first filter stage that reduced candidates to zero, and the full `filter_trail`). A missing result is turned into a `build_request.LibraryBuildRequest` -- never a fabricated, interpolated, or rescaled match (Sec. 21).

