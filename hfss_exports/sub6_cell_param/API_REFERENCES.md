# Documented APIs and implementation status

User decision 2026-11-10 explicitly authorizes PyAEDT scripts and copying them, superseding the
older handoff's STEP/manual-only and no-Python-copy instructions for this package.
No recorded scripts or guessed AEDT command arrays are used.

PyPI distribution: [pyaedt](https://pypi.org/project/pyaedt/), import `ansys.aedt.core`, installed release 1.9.0.
The package is not distributed under the requested descriptive name `ansys-aedt-core`.

Documented interfaces used:

- [Hfss](https://aedt.docs.pyansys.com/version/stable/API/_autosummary/ansys.aedt.core.hfss.Hfss.html):
  non_graphical/version/project/design/solution_type/new_desktop; variable assignment;
  materials.add_material and permittivity/dielectric_loss_tangent; modeler.model_units.
- Modeler3D.create_box, create_cylinder, create_rectangle, subtract, unite, get_faceid_from_position;
  Object3d.name, solve_inside, faces; FacePrimitive.is_planar.
  [3D modeler](https://aedt.docs.pyansys.com/version/stable/API/Modeler.html).
- [Hfss.assign_impedance_to_sheet](https://aedt.docs.pyansys.com/version/stable/API/_autosummary/ansys.aedt.core.hfss.Hfss.assign_impedance_to_sheet.html).
- [Hfss.assign_lattice_pair](https://aedt.docs.pyansys.com/version/stable/API/_autosummary/ansys.aedt.core.hfss.Hfss.assign_lattice_pair.html),
  InputPhaseDelay with phase_delay_param1; Hfss.assign_perfect_e.
- [Hfss.assign_radiation_boundary_to_faces](https://aedt.docs.pyansys.com/version/stable/API/_autosummary/ansys.aedt.core.hfss.Hfss.assign_radiation_boundary_to_faces.html):
  **documented call exists, but the installed AEDT rejects Radiation in Eigenmode.**
- Mesh.assign_length_mesh, maximum_elements=None; Hfss.create_setup, SetupHFSS.props/update;
  Hfss.validate_full_design, save_project, cleanup_solution, export_profile, release_desktop.
- [Hfss.analyze_setup](https://aedt.docs.pyansys.com/version/stable/API/_autosummary/ansys.aedt.core.hfss.Hfss.analyze_setup.html):
  cores, use_auto_settings, blocking. Request every logical CPU; actual license capability unverified.
- PostProcessor3D.available_report_quantities/get_solution_data; SolutionData.get_expression_data/units_data.
  Retrieval follows the [official eigenmode example](https://examples.aedt.docs.pyansys.com/version/dev/examples/high_frequency/emc/eigenmode.html)
  categories `Eigen Q`, `Eigen Modes`, report category `Eigenmode`, without hardcoded quantity names.

Actual execution reached impedance, lattice pairs and PEC lid assignment. It stopped on Radiation;
mesh/setup/solve/export calls after that line have not been exercised on this server.
The same-version user HFSS calibration cannot be claimed reproduced until its actual bottom boundary is known.

The [2026 R1 Radiation documentation](https://ansyshelp.ansys.com/public/views/secured/electronics/v261/en/Subsystems/HFSS/Content/HFSS/AssigningRadiationBoundaries.htm)
mentions Eigenmode, but this does not override the 2025 R1 runtime rejection.
The [2025 R1 PML documentation](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/HFSS/AssigningPMLBoundaries.htm)
also states PML radiation boundaries are not generated in Eigenmode projects; PML is not an assumed replacement.
