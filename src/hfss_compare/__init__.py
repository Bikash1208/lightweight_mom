"""Part 3, app page 7 ("HFSS compare"): upload a Touchstone (.sNp) file, pick the matching MoM
case, overlay vs MoM, error metrics -- reusing `phase9.hfss_reference.read_touchstone` and
`phase12d.hfss_ingest.compare_single_cell` unmodified (both already generic over an arbitrary
Touchstone path). This package adds only a small on-disk ledger so ad-hoc app-side comparisons
(beyond the fixed G1-G6 set Phase 12D/13 already recorded) are remembered across app restarts."""
