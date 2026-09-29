"""Sec. 26-27: `database_schema.json`/`.md`, `query_spec.md`, `library_build_request_schema.json`.
Every schema below is generated from the actual dataclass fields in this package (via
`dataclasses.fields`), so the written files can never drift out of sync with the code that
produces the data they describe.
"""

import dataclasses
import json

from .build_request import LibraryBuildRequest
from .family import UnitCellFamily
from .library_entry import ROLES
from .query import UnitCellQuery
from .record_schema import UnitCellRecord
from .validation_levels import VALIDATION_LEVELS


def _field_map(cls):
    return {f.name: str(f.type) for f in dataclasses.fields(cls)}


FORMAT_NAME = "lightweight_mom.universal_unit_cell_database"
FORMAT_VERSION = 1


def build_database_schema_dict():
    return {
        "format": FORMAT_NAME,
        "version": FORMAT_VERSION,
        "record_schema": _field_map(UnitCellRecord),
        "family_schema": _field_map(UnitCellFamily),
        "validation_levels": list(VALIDATION_LEVELS),
        "library_entry_roles": list(ROLES),
        "note": (
            "Period (p_x_mm, p_y_mm), frequency, substrate, polarization and incidence are all "
            "per-LibraryEntry fields, never global constants; a family aggregates them across "
            "every registered entry of one topology (see family.py)."
        ),
    }


def write_database_schema_json(path):
    doc = build_database_schema_dict()
    path.write_text(json.dumps(doc, indent=1) + "\n")
    return path


def render_database_schema_md():
    doc = build_database_schema_dict()
    lines = [
        "# database_schema.md", "",
        "Generated from the actual `phase8d_rev` dataclasses (`dataclasses.fields`), so this file "
        "cannot drift from the code.", "",
        "## UnitCellRecord (Sec. 7)", "",
    ]
    for name, typ in doc["record_schema"].items():
        lines.append(f"- `{name}`: `{typ}`")
    lines += ["", "## UnitCellFamily (Sec. 9)", ""]
    for name, typ in doc["family_schema"].items():
        lines.append(f"- `{name}`: `{typ}`")
    lines += [
        "", "## Validation levels (Sec. 23)", "",
        ", ".join(f"`{v}`" for v in doc["validation_levels"]), "",
        "## LibraryEntry roles (Sec. 3, 17)", "",
        ", ".join(f"`{r}`" for r in doc["library_entry_roles"]), "",
        "## Period/frequency/substrate/polarization/incidence are per-entry, not global", "",
        doc["note"], "",
    ]
    return "\n".join(lines) + "\n"


def render_query_spec_md():
    fields = _field_map(UnitCellQuery)
    lines = [
        "# query_spec.md", "",
        "`query.search_unit_cells(universal_db, query)` applies the following filters in this "
        "fixed order (Sec. 14); it never ranks surviving candidates by phase coverage or any other "
        "controllability metric.", "",
        "## UnitCellQuery fields", "",
    ]
    for name, typ in fields.items():
        lines.append(f"- `{name}`: `{typ}` (`None` = unconstrained)")
    lines += [
        "", "## Filter order", "",
        "1. topology", "2. material_stack (substrate name/eps_r/tan_delta/thickness)",
        "3. frequency (exact stored grid point only, no interpolation)",
        "4. period (p_x_mm, p_y_mm)", "5. polarization", "6. incidence",
        "7. response_availability (at least one OK-solved cell)",
        "8. validity_status (>= min_validation_status, if given)", "",
        "## Result", "",
        "`QueryResult.status` is `\"MATCH\"` (with the surviving `LibraryEntry` objects) or "
        "`\"DATABASE_COVERAGE_MISSING\"` (with `eliminated_at` naming the first filter stage that "
        "reduced candidates to zero, and the full `filter_trail`). A missing result is turned into "
        "a `build_request.LibraryBuildRequest` -- never a fabricated, interpolated, or rescaled "
        "match (Sec. 21).", "",
    ]
    return "\n".join(lines) + "\n"


def build_library_build_request_schema_dict():
    return {
        "format": "lightweight_mom.library_build_request",
        "version": 1,
        "fields": _field_map(LibraryBuildRequest),
        "note": "Produced only by build_request.generate_build_request on a DATABASE_COVERAGE_MISSING "
                "query result; never auto-submitted to a solver (Sec. 18: 'must NOT silently generate it').",
    }


def write_library_build_request_schema_json(path):
    path.write_text(json.dumps(build_library_build_request_schema_dict(), indent=1) + "\n")
    return path
