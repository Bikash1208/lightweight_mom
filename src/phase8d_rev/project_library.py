"""Sec. 17: the distinction between the long-lived `UniversalDatabase` and the small
`project_library` a future HMSA design actually uses -- the latter is *extracted* from the former,
never a separate, independently-maintained dataset.
"""

from dataclasses import dataclass

from .query import UnitCellQuery, search_unit_cells


@dataclass(frozen=True)
class ProjectLibrary:
    query: UnitCellQuery
    source_library_ids: tuple  # LibraryEntry.library_id values this subset was extracted from
    entries: tuple  # the matching LibraryEntry objects themselves (read-only refs)


def extract_project_library(universal_db, query: UnitCellQuery):
    """Sec. 17: run ``query`` against ``universal_db`` and, if it matches, return a
    `ProjectLibrary` referencing (not copying) the surviving entries. Raises ValueError if the
    query does not match anything -- callers must handle `DATABASE_COVERAGE_MISSING` via
    `query.search_unit_cells`/`build_request.generate_build_request` themselves before calling
    this (Sec. 18: a missing project library produces a build request, not an empty library)."""
    result = search_unit_cells(universal_db, query)
    if result.status != "MATCH":
        raise ValueError(
            f"extract_project_library: query did not match ({result.status} at "
            f"'{result.eliminated_at}'); use build_request.generate_build_request(result) instead"
        )
    return ProjectLibrary(
        query=query,
        source_library_ids=tuple(e.library_id for e in result.matching_entries),
        entries=result.matching_entries,
    )
