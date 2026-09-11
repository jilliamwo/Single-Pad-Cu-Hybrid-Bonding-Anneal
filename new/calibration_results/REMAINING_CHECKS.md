# Remaining checks: results

All locally executable solver assertions passed. The existing ANSYS CSVs still
do not establish elastic-only parity. No thermal calibration or boundary
changes were applied to production defaults.

## Solver and numerical checks

- The existing zero-load, thermal scaling, energy scaling, fixed-boundary,
  finite energy, normalization-length, and analytic free-expansion checks pass
  on two MPI ranks.
- Serial and two-rank runs on the same 12 x 26 cell mesh agree in Cu surface
  displacement to 2.1e-14 nm. Statistics agree within 1e-10 relative tolerance,
  excluding timing. See `mpi_check.json`.
- A homogeneous cylinder with radial confinement matches the exact solution
  ur=0, uz=alpha*dT*(1+nu)/(1-nu)*z. Maximum dimensionless displacement error
  is below 9e-16; radial stress matches -E*alpha*dT/(1-nu), and axial stress
  is zero within 0.1 Pa.
- The zero-recess geometry produces finite displacements and a continuous
  Cu/dielectric edge height.
- Full Cu-localized energy and the corrected ANSYS-comparison energy agree.

With Cu E=129 GPa, nu=0.34, and the provisional 20 C thermal reference:

| Radial cells | Vertical cells | Cu-center uz (nm) | Cu mean W (J/m³) | Peak W (J/m³) |
|---|---:|---:|---:|---:|
| 40 | 82 | 2.190619 | 537571 | 2.792e6 |
| 80 | 162 | 2.191888 | 535501 | 3.512e6 |
| 160 | 322 | 2.192245 | 534795 | 4.393e6 |

Center displacement changes only 0.000357 nm in the final refinement, versus
about 0.629 nm discrepancy from ANSYS. Mean energy changes about 0.13% in the
final refinement. Peak energy is still not converged; corner/interface peaks
should not be used to tune material parameters.

Diagnostic boundary variants at resolution 40 gave Cu-center uz:

| Constraint variant | Center uz (nm) |
|---|---:|
| Confirmed constraints | 2.190619 |
| Bottom radial displacement free | 2.190551 |
| Outer radial displacement free | 1.867688 |
| Both radial constraints free | 1.859006 |

These variants do not fix the mismatch and were not promoted to defaults.
The exported bottom uz is zero in all 11 files. Outer-side uz is nonzero,
consistent with leaving axial displacement free there. Radial constraints
cannot be verified directly because radial displacement was not exported.

## All 11 ANSYS CSVs audited

The older CSVs use different column names/orders. The loader now recognizes
all three observed layouts and UTF-8 BOM headers. All input values parsed as
finite numbers. Elastic constants below are inferred from interior normal
stress/elastic-strain pairs, not copied from material labels:

| Export group | Cu E (GPa) | Cu nu | Dielectric E | Dielectric nu |
|---|---:|---:|---:|---:|
| SiO2 and temperature series | 129 | 0.34 | 73 GPa | 0.17 |
| BCB | 129 | 0.34 | 2.9 GPa | 0.34 |
| SiCN | 129 | 0.34 | 150 GPa | 0.25 |
| SiCOH | 129 | 0.34 | 6 GPa | 0.22 |
| Polyimide A | **117** | **0.30** | 1.7 GPa | 0.34 |
| Polyimide B | **117** | **0.30** | 3.2 GPa | 0.34 |
| diamond | 129 | 0.34 | **1050 Pa** | 0.07 |

Two findings require checking before comparing the dielectric sweep:

1. The Polyimide files also change the Cu elastic constants, so this is not a
   sweep with identical Cu inputs.
2. The diamond stress/strain columns imply E=1050 Pa in both dielectric
   regions, with relative constitutive residuals about 1–2e-5. Check the actual
   material entry and export units. This is evidence of a possible input/unit
   error; it does not by itself identify where that error occurred.

All energy columns agree with full elastic energy using engineering shear.
The older `250C_500nm_3.5recess.csv` differs by roughly 8e-6 relative, consistent
with rounded exported components; other files reconstruct exactly at parsed
precision. Therefore the zz-only energy reconstruction is unsuitable here.

Nine files contain 221 repeated node IDs. At coincident coordinates rounded
to 1e-5 nm, the SiO2 case has 247 repeated-coordinate groups. Along the Cu side,
paired uz values agree; along the base interface the largest difference is
0.004754 nm. These small base differences could reflect interface treatment
or export details; they do not explain a 0.629 nm center-surface discrepancy.
The calibration script fits unique surface radii. Additional sensitivity RMSEs
in `remaining_checks.json` are computed over surface rows, so duplicated rows
receive repeated weight; use the original calibration report for held-out metrics.

## Temperature-reference check

Using only the top-axis intervals (980–997 nm) of the 200/250/300/350 C
exports, a linear fit of duz/dz - elastic_zz versus nominal filename temperature
gives alpha≈16.917e-6/K and an apparent reference temperature≈20.83 C.
Residual strain RMSE is 6.3e-6.

This independently supports the scale of the current Cu CTE and a reference
near room temperature. It does not confirm an exact 20 or 22 C setting: it
assumes the filenames represent solved temperatures and that these top-axis
intervals have negligible additional inelastic strain. Nodal averaging and
finite differences also contribute uncertainty. The result does not support
the much larger Cu thermal strain needed by the rejected surface-only fit.

## Checks requiring unavailable ANSYS information

The repository contains no ANSYS project or result database. The CSVs contain
no total strain, equivalent plastic strain, temperature field, radial
displacement, or result-time metadata. Therefore the following remain unchecked:

- Whether hardening was active in the material actually assigned to Cu.
- Whether total strain decomposes into elastic plus thermal strain only.
- Whether all CSV columns came from the same result time/state.
- The actual temperature field and stress-free temperature.
- Exact contact/interface definitions and radial constraints in the solved model.

The previously reported stress plateau and axial compatibility discrepancy
remain evidence to investigate, not proof that plasticity was active.

## Reproduce

```sh
python new/run_remaining_checks.py
mpiexec -n 2 python new/validate_baseline.py --nr 12 --nz 12 --check --output /tmp/bonding-mpi-check
python new/validate_baseline.py --nr 12 --nz 12 --output /tmp/bonding-serial-check
```

Detailed export, mesh, boundary, and analytic results are in
`remaining_checks.json`. The two-rank validation was run in the installed
FEniCSx 0.10.0 environment.
