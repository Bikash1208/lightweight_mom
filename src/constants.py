"""Physical constants and global conventions.

All quantities are SI (metres, hertz, seconds, radians). See CONVENTIONS.md
for the full list of sign and coordinate conventions; the ones that affect
this module are:

* Time convention: exp(+j*omega*t). A phasor F represents Re{F exp(+j w t)}.
* Consequently a lossy medium has Im(eps) < 0 and Im(k) < 0.
"""

# Speed of light in vacuum [m/s]. Exact by SI definition.
C0 = 299_792_458.0

# Vacuum permeability [H/m]. CODATA 2018 measured value (mu0 is no longer
# exactly 4*pi*1e-7 after the 2019 SI redefinition).
MU0 = 1.25663706212e-6

# Vacuum permittivity [F/m], defined here from MU0 and C0 so that
# EPS0 * MU0 * C0**2 == 1 holds to rounding error.
EPS0 = 1.0 / (MU0 * C0**2)

# Free-space wave impedance [ohm] = sqrt(MU0/EPS0) = MU0*C0.
ETA0 = MU0 * C0

# Documented global time convention. Nothing in the code branches on this
# string; it exists so that reports/exports can state the convention.
TIME_CONVENTION = "exp(+jwt)"
