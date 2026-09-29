"""Part 3, app page 4 ("Pattern prediction"): far-field pattern, beam pointing, SLL and a
directivity estimate from an aperture field, reusing `phase8e.beam_reconstruction.far_field_pattern`
(the standard planar-aperture radiation integral) unmodified. Always labelled
`LOCAL_PERIODICITY_APPROXIMATION` by the caller: the aperture field itself, for a reflectarray or
HMSA design, comes from a per-cell lookup-table assumption, not a full-wave solve of the whole
array."""
