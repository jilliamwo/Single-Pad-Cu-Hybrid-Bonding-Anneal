# Native thermoelastic baseline

`thermoelastic_solver.py` solves small-strain, axisymmetric elasticity or optional
J2 plasticity with prescribed thermal loading using FEniCSx 0.10.0. ANSYS parity is not established.
The archived phase-field code under `old/` is preserved for provenance.
For the active coupled model and terminal workflow, see [RUNNING_THE_MODEL.md](../RUNNING_THE_MODEL.md).

## Reference geometry and constraints

All dimensions are parameters of `ReferenceGeometry`. Defaults are:

- SiO2 base: r = 0–500 nm, z = 0–500 nm.
- SiO2 annulus: r = 250–500 nm, z = 500–1000 nm.
- Cu: r = 0–250 nm, z = 500–996.5 nm.
- Recess: r = 0–250 nm, z = 996.5–1000 nm, modeled with residual stiffness.
- Axis ur = 0; bottom ur = uz = 0; outer side ur = 0.
- Remaining exterior components are traction-free.

`MechanicalBCs` exposes these constraints. The default mesh contains rows at
both Cu and dielectric surfaces and at the base interface. `--nz` controls each
of the base and pad layers; two additional rows of cells resolve the recess.

Cu and SiO2 elastic constants are supported by least-squares fits of interior
normal stresses to elastic strains in the SiO2 CSV. Thermal inputs still require
confirmation against the ANSYS input deck:
Cu E=129 GPa, nu=0.34, alpha=17e-6/K; SiO2 E=73 GPa, nu=0.17,
alpha=0.55e-6/K; stress-free temperature 293.15 K, anneal 523.15 K.

## Run

From the repository root with the FEniCSx environment active:

```sh
python new/validate_baseline.py --check \
  --ansys old/ansys-results/250_500nm_3.5_SiO2.csv \
  --output /tmp/bonding-baseline-validation
```

Outputs are a JSON summary, physical Cu/dielectric surface CSVs, and an optional
ANSYS surface comparison. Comparison uses nodes on the reference Cu surface,
without a displacement offset, and interpolates FEniCSx uz onto ANSYS radii.

Checks cover zero thermal load, linear displacement/quadratic energy scaling,
fixed boundary values, finite nonnegative energy, normalization-length
invariance, and homogeneous analytic free expansion with negligible stress.

## API and limitations

`solve` requires an explicit `geometry` and the full mechanics mesh. The
Cu-only phase-field mesh does not contain the dielectric and cannot be substituted.
`static_coupling.py` transfers fields using the explicit coordinate map
z_phase = z_reference - base_thickness/L0.

`u_dimless_field` and `UY_dimless` are normalized by L0; `u_field` is in metres.
`cu_surface_m` and `dielectric_surface_m` contain deformed physical (r,z).
`W_field` is Cu-localized energy; `W_full_field` includes every material.
Energy and stress outputs use DG0 cell fields to avoid smoothing jumps across
material interfaces. `W_ansys_like_field` now uses full Cu-localized elastic
energy: the SiO2 CSV's W exactly equals half the sum of all three normal
stress/elastic-strain products plus shear stress times engineering shear strain.
The older zz-only reconstruction described in the context does not match this CSV.


Dielectric identity is fixed in reference coordinates while its displacement is
solved. Optional scalar phi on the mechanics mesh changes Cu/void constitutive
properties through a clipped smooth interpolation, independently of dielectric.
Surface diagnostics currently describe the flat reference baseline surfaces;
they do not extract evolving phi contours. Dynamic coupling, CMP variants,
creep, transient heating, and contact are not implemented. Unsupported transient
and contact requests raise errors. Residual void stiffness is a numerical
regularization, not a contact law. Mesh and stiffness convergence must be assessed
for each validated case, especially near interface corners.

## Initial validation results (previous Cu defaults: 110 GPa, nu=0.343)

On the 40 x 82 baseline mesh, all `--check` assertions passed in FEniCSx 0.10.0.
Cu-center uz was 2.227950 nm versus 2.821300 nm in the SiO2 ANSYS CSV;
the surface RMSE across 127 ANSYS radii was 0.505199 nm. Cu-edge uz was
1.501476 nm versus 1.251000 nm in ANSYS.

Refining to 80 x 162 cells gave center uz = 2.229539 nm. At the original
resolution, sweeping residual stiffness from 100 to 100000 Pa changed center
uz by less than 4e-9 nm. These checks do not explain the ANSYS discrepancy.
Peak cell W increased from 2.659e6 to 3.303e6 J/m³ with refinement, so peak
energy is not converged. Confirm ANSYS material constants, thermal reference
state, and exact model constraints before claiming replacement parity.

## Elastic constants identified from the ANSYS export

Fit sigma_i = lambda*trace(epsilon_el) + 2*mu*epsilon_el_i jointly to
all three normal components, excluding nodes near material interfaces and
exterior boundaries. The 250 °C SiO2 CSV gives:

| Region | Nodes | E (GPa) | nu | Relative stress fit residual |
|---|---:|---:|---:|---:|
| Cu interior | 2232 | 129.000040 | 0.339999974 | 1.69e-5 |
| SiO2 base interior | 505 | 72.999835 | 0.169999480 | 2.47e-5 |
| SiO2 annulus interior | 1601 | 72.999997 | 0.170000536 | 2.12e-5 |

Masks in nm: Cu r<240, 510<z<986; base r<490, 10<z<490;
annulus 260<r<490, 510<z<990. These fits use normal elastic strains only,
so the shear-strain convention does not affect the inferred constants.

The default Cu constants are now 129 GPa and 0.34. With these values, the
40 x 82 run passes all analytic/consistency checks but gives center uz =
2.190619 nm, edge uz = 1.532893 nm, and ANSYS surface RMSE = 0.520676 nm.
Thus the elastic-constant correction does not resolve the mismatch. Thermal
expansion and stress-free temperature are not determined by this stress/elastic-
strain fit and remain provisional. Those comparisons used the elastic native
solver. The actual ANSYS material-law settings remain unverified; optional
plasticity has since been added as described below. Creep remains unimplemented.

## Calibration audit

See [the mismatch/calibration report](calibration_results/REPORT.md). The
reproducible `calibrate_baseline.py` fits effective thermal strains and checks
surface holdouts plus interior displacement/stress. Its attempted thermal fits
are saved as rejected candidates; thermal solver defaults remain provisional.
The export shows evidence inconsistent with uniform-temperature elastic-only
Cu deformation, and its W column is verified to be full elastic energy.

[Remaining checks](calibration_results/REMAINING_CHECKS.md) cover all 11 source
CSVs, serial/MPI agreement, radial-confinement verification, additional mesh
refinement, zero recess, and boundary sensitivity. The audit identifies distinct
Cu elastic constants in the Polyimide cases and an apparent 1050 Pa dielectric
modulus in the diamond export; review those inputs before comparing the sweep.

## Calibration across cases

[The temperature/material sweep report](calibration_results/SWEEP_REPORT.md)
contains native FEA fits for all 11 exports and held-out-temperature predictions.
SiCOH and Polyimide A/B have provisional fits checked on the 160 x 322 mesh;
[their parameter profiles](calibration_results/provisional_material_profiles.json)
assume a 22 C reference and have only been checked at 250 C. The SiO2 thermal
calibration does not transfer accurately through 350 C. The older unlabeled
250 C export shows evidence of different outer-boundary behavior and is tested
separately. Run `calibrate_sweep.py --help` for the sweep and diagnostic options.

## Optional plasticity

The [plasticity report](calibration_results/plasticity/REPORT.md) documents the
new Cu J2 model and direct comparisons across all 11 ANSYS exports. The default
solver remains elastic; enable the material-sheet Cu law explicitly:

```python
from plasticity import BilinearPlasticity
# Supply all existing mesh/material/geometry arguments as before.
result = solve(..., plasticity_cu=BilinearPlasticity(), load_steps=40)
```

Defaults in `BilinearPlasticity` are yield=172.3 MPa and uniaxial tangent
modulus=517 MPa. Calibration runs explicitly assume a 22 C reference and use a
uniform, quasi-static thermal ramp. This does not reconstruct the ANSYS transient
thermal field or add inertia/creep. Dielectrics remain elastic unless a separate,
justified `plasticity_dielectric` law is supplied.

Plastic strains are stored at quadrature points; nonlinear increments commit
history only after convergence. Each call integrates from a virgin state and
returns `plastic_state` plus a DG0 `equivalent_plastic_strain_field`. Repeated
calls do not continue prior history. Evolving-phi plasticity is not supported yet.

```sh
python new/validate_plasticity.py
python new/calibrate_plasticity.py --resolution 40 --steps 20
```
