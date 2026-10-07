# Inputs

Example input files used by the app and the docs/examples that walk through it.

| File | What it is |
|---|---|
| `reference_fields/efl_phase_EzZ_example.fld` | An HFSS field-calculator export of the monopole feed's reference surface-wave field (phase only — the Z-component of E, exported along a plane) at 10.54 GHz. The worked example `docs/USER_GUIDE.md`'s Step 1 and `docs/START_HERE.md` both point to: upload it, pick a substrate and cell size, and walk the 3-step wizard. Carries one documented feed-point NaN, which the app interpolates and reports (`FEED_REGION_MASKED`) rather than dropping silently. |

Upload your own `.fld` export the same way — nothing about this file is special-cased beyond
being the one example checked into the repo.
