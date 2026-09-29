"""Part 3, app page 3 ("Reflectarray designer"): the required-phase-map formula (standard
reflectarray theory, e.g. Huang & Encinar, *Reflectarray Antennas*, Ch. 1) and cell selection from
a chosen library, reusing `phase8d_rev.response_matching.nearest_response_match` unmodified. Every
design this produces carries `LOCAL_PERIODICITY_APPROXIMATION` and, where the per-cell phase
gradient exceeds Phase 12D's own measured limit, `STEEP_GRADIENT(...)`.
"""
