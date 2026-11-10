# BASELINE / G1 island-cell HFSS packages

MODEL_ONLY / NOT_HFSS_VALIDATED. Four geometry packages, each with one STEP
shared by OFF/ON, one-page HFSS_SETUP.md, states.csv, model references and
an eigenmode results template. No AEDT scripts.

| Package | Physical slot L x W, mm | Physical island L x W, mm | Hardware |
|---|---:|---:|---|
| BASELINE_novia | 8.696 x 2.900 | 4.348 x 1.300 | no vias / no FR4 |
| BASELINE_vias_fr4_stubs | 8.696 x 2.900 | 4.348 x 1.300 | two vias + FR4 + calibration stubs |
| G1_novia | 9.184 x 2.658 | 3.866 x 1.054 | no vias / no FR4 |
| G1_vias_fr4_stubs | 9.184 x 2.658 | 3.866 x 1.054 | two vias + FR4 + calibration stubs |

Disc radius 5.435, period 11.600, RT5880 3.175, copper 0.035 mm. Physical
diode sheet is exactly 1.000 x 0.350 mm. BASELINE copies the physical board
builder dimensions, not the discretized MoM island outline. G1 vertices
are snapped to 0.001 mm (half-slot length 4.592, half-island length 1.933,
slot half-width 1.329, island y=-0.725..+0.329). Differences from the G1
search's continuous input dimensions are <=0.001 mm; MoM realized dimensions
are separately recorded per mesh in results/sub6_island_meshtrend/.

All STEP shapes are rotated +90 degrees to match ISL_vias2: across=x,
along=y. User-authorized via variant retains both 2.0 x 0.3 mm +x-directed
stubs in the island frame (+y after rotation), on the FR4 bottom face.
Via centres in island frame are (-1.780,-0.530), (-3.630,-3.630); barrel
diameter 0.3, bottom pad 0.6, anti-pad 0.6, FR4 thickness 0.8 mm.
G1 island-to-barrel x-edge clearance is only 0.003 mm: exact geometry fits,
but this is not a fabrication-clearance qualification. No via relocation
or top-pad enlargement was introduced.

AIR is a single overlapping six-face box. Top PEC lid is at +150 mm.
Bottom is 30 mm below the lowest copper, with Radiation, retaining the
calibration cavity fix. This user-requested lower extension supersedes
the generic ground-bottom-only AIR limit in EXPORT_STANDARD; all other
export gates apply. RT5880/FR4 contain explicit barrel holes to avoid
metal/dielectric solid overlaps. The two anti-pad voids are air in the
ground-thickness band, as in the archived calibration stack.

## Method and readouts

See CALIBRATION_METHOD.md. Read f, Q and phase of the identified top-surface
TM branch for both directions and both states. From a local f-vs-phase
curve at 5.8 GHz obtain beta and group velocity, then alpha_x=omega/(2Qvg).
At least three HFSS phase points are requested to determine the slope
without importing the mesh-unreliable MoM ON dispersion into the result.

The old calibration used parallel-RLC equivalents; these packages use the
user's specified Impedance boundary values at 5.8 GHz. Record that model
difference. A fixed-Rs/Xs phase sweep measures dispersion of that specified
frozen-boundary HFSS model; it must not silently be labelled a calibration
of the full frequency-dependent series-L diode circuit. Use the same
recorded boundary convention for no-via/via comparisons.

The uniform-cell Q-derived alpha_x is loss attenuation. It does not measure
the first-order radiation leakage of a rotating-axis hologram. Those
leakage-alpha numbers are included only as labelled mesh-trend references.
Directional no-via beta and perturbative sheet-loss estimates, all roots,
and extrapolation flags are in predictions.json/model_reference.csv.
No 3D via/FR4/stub prediction is available from this MoM model; those
entries are explicitly NOT_AVAILABLE_FROM_CURRENT_DATA rather than reused
no-via numbers. No current-diode HFSS results are claimed.

The 150-mm PEC lid matches the calibration setup. Near-light-line ON
roots can have appreciable field at that lid; predictions.json records
the model exp(-kappa*0.150) field-tail estimate per root. Check Ez at the
lid and flag AIR_HEIGHT_NOT_CONVERGED if it is not small; do not describe
a lid-sensitive result as an open-air calibration without a height check.

## Gates

geometry_summary.json and each geometry_checks.json contain pre-export and
re-import results: names, valid BRep, three metal solids in one object,
no solid overlap, six-face AIR, physical top mask/volume, single exact
sheet with coincident copper edges, and EXACT_DIMS. STEP SHA256 is recorded.
No HFSS solution has been run by this session.
