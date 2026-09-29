# s_parameter_definition.md

## Ports

Port 1 = top reference plane (z_top = +top_air_thickness), incident wave travels -z (toward the pattern). Port 2 = bottom reference plane (z_bottom = -(h + bottom_air_thickness)), incident wave travels +z. Only the specular (0,0) Floquet order at normal incidence is a port mode (`CONVENTIONS.md` 'Ports').

## S-matrix

`4x4 [[S11, S12], [S21, S22]] over (port1 x, port1 y, port2 x, port2 y); entry name S<out><in>_<received><incident>; port 1 = top reference plane (wave incident travelling -z), port 2 = bottom reference plane (wave incident travelling +z)`

## Polarization basis

`global x and y unit vectors, identical at both ports`

Circular basis (derived, not independently excited/stored): R=(x + j y)/sqrt(2), L=(x - j y)/sqrt(2).

## Normalization

Wave amplitudes are tangential-E amplitudes at the reference plane, *measured* from E and H at the plane (not assumed unit-normalized, `CONVENTIONS.md` 'Waves, ports and S-parameters'); incident amplitude is 1 V/m referenced at the pattern plane z=0 for the excitation itself (`excitation.py`).

## Reference-plane sensitivity

Moving the reference planes changes S11/S22 by a pure phase factor and S21/S12 by another (`CONVENTIONS.md`: `S11' = S11 exp(-2j k0 D1)`, `S22' = S22 exp(+2j k0 D2)`, `S21' = S12' = S21 exp(+j k0 (D2-D1))`), so phase is only meaningful relative to the specific reference planes stored in the database metadata -- see `reference_plane_sensitivity.py`.
