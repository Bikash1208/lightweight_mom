# unit_cell_topology_audit.md

## Ground plane

**ABSENT**

context['substrate'] = {'name': 'Rogers RT/Duroid 5880', 'eps_r': 2.2, 'tan_delta': 0.0009, 'thickness_m': 0.001575} (one dielectric layer, not a stack terminated in a conductor); context['metal'] = {'model': 'pec'} (a patterned PEC *sheet at z=0*, i.e. a finite scatterer, not a continuous ground plane); no 'ground' key exists anywhere in the context dict (checked: False); `config.py`/`materials.py`/`geometry.py`/`periodic_bc.py` contain zero occurrences of 'ground' (grep-verified) -- the config schema cannot express a ground plane at all. README.md 'Assumptions' item 2: 'No ground plane, no superstrate' (Phase 1, unchanged through Phase 7).

## Metal layers

one patterned PEC sheet at z=0 (context['metal']={'model': 'pec'})

## Substrate layers

one dielectric layer (context['substrate']={'name': 'Rogers RT/Duroid 5880', 'eps_r': 2.2, 'tan_delta': 0.0009, 'thickness_m': 0.001575})

## Air regions / reference planes

two semi-infinite air half-spaces (z>0 and z<-h), each with an explicit reference plane (context['reference_planes']={'top_air_thickness_m': 0.005, 'bottom_air_thickness_m': 0.005})

## Physical architecture classification

**PARTIALLY_REFLECTIVE_TRANSMISSIVE** -- more transmissive than reflective (mean reflected power 0.1922, transmitted 0.8070)

## Ports, excitation, polarization, S-matrix definition

See `s_parameter_definition.md` (this project's fixed, existing conventions, `CONVENTIONS.md`/`S_PARAMETERS.md`, reused verbatim -- not redefined here).
