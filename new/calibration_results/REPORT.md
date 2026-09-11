# ANSYS mismatch audit and attempted thermal calibration

The present elastic model cannot be calibrated to reproduce both displacement
and stress in `old/ansys-results/250_500nm_3.5_SiO2.csv` by adjusting uniform
Cu and SiO2 thermal strains. No fitted thermal parameters were promoted to
solver defaults. The export has strong signatures consistent with inelastic
Cu deformation or a mismatch between exported result states. This conflicts
with the intended elastic-only setup and needs verification in the ANSYS model.
It does not establish that plasticity was enabled.

## Confirmed findings

Normal stresses fit isotropic Hooke's law applied to the exported **elastic**
strains with relative residuals below 0.003%:

| Material | E | nu |
|---|---:|---:|
| Cu | 129 GPa | 0.34 |
| SiO2 base and annulus | 73 GPa | 0.17 |

This constitutive check confirms elastic constants, but does **not** prove that
total strain was purely elastic plus thermal. An elastoplastic solve also
obeys the elastic constitutive relation when its elastic strains are used.

Cu shear stress divided by exported shear strain gives 48.1343 GPa, matching
G = E/[2(1+nu)], so the exported shear strain is engineering gamma, not tensor
epsilon_rz. The CSV energy column matches exactly:

    W = 0.5 * (sigma_x*epsilon_x_el + sigma_y*epsilon_y_el
               + sigma_z*epsilon_z_el + tau_xy*gamma_xy_el)

The old zz-only formula has 66% relative error across this CSV. The solver's
`W_ansys_like_field` now uses the verified full energy definition. This change
affects the comparison diagnostic, not displacement or `W_field` mechanics.

## Why uniform thermal calibration fails

On the symmetry axis, the small-strain identity is:

    du_z/dz - epsilon_zz_el = thermal strain + other inelastic strain

Finite differences of exported displacement and interval-averaged elastic
strain give these median residuals inside Cu:

| Axis height interval | Residual axial strain |
|---|---:|
| 525–750 nm | 0.006588 |
| 750–900 nm | 0.005741 |
| 900–950 nm | 0.004788 |
| 970–997 nm | 0.003879 |

For constant CTE and uniform temperature in a purely thermoelastic material,
this residual should be constant. The difference is much larger than output
rounding, though nodal averaging and displacement differentiation introduce
uncertainty. A nonuniform temperature or inconsistent result-time export is
another possibility. The top value is compatible with 17e-6/K times about
228 K, suggesting a 22 C reference **if** that CTE and 250 C temperature apply.
This is a hypothesis, not an identified stress-free temperature.

About 50% of interior Cu nodes have von Mises stress between 170 and 180 MPa.
The 90th percentile is 174.27 MPa. The supplied material sheet defines a
172.3 MPa yield strength and 517 MPa bilinear tangent modulus. This clustering,
together with excess axial strain, supports checking whether that hardening
law was active in the exported run despite the intended elastic baseline.
ANSYS documents that this law produces plastic strain beyond yield:
[ANSYS bilinear isotropic hardening](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/ans_arch/archbisodesc.html).

As an additional conditional check, assuming thermal strain 0.003876 and
approximately monotonic axisymmetric J2 flow, the inferred axial plastic strain
between 510 and 850 nm is about 0.0018–0.0029. The material-sheet hardening law
sigma_y = 172.3 MPa + [E*Et/(E-Et)]*equivalent_plastic_strain predicts the
interval stress differences with about 0.25 MPa RMSE. This is supporting
evidence, not a reconstructed loading history or a substitute for plastic-strain
output. The flow assumption can fail during unloading/nonproportional loading.

## Fits and independent checks

Two elastic basis solves apply thermal strain 0.001 separately to Cu and SiO2.
Linear superposition allows fitting their amplitudes without changing stiffness,
geometry, or the confirmed boundary conditions. Alternating surface radii train
the surface fit; the remaining radii are held out. These are spatial holdouts
from the same simulation, not independent experimental validation. Interior uz
and all four stress components provide additional checks.

Results on the refined 80 x 162 cell mesh:

| Candidate | Cu thermal strain | SiO2 thermal strain | Held-out surface RMSE (nm) | Interior stress relative L2 error |
|---|---:|---:|---:|---:|
| Current 20 C reference | 0.003910 | 0.0001265 | 0.519 | 25.0% |
| 22 C reference hypothesis | 0.003876 | 0.0001254 | 0.536 | 24.3% |
| Unconstrained surface fit | 0.006006 | -0.0003339 | 0.065 | 82.6% |
| Nonnegative surface fit | 0.005271 | 0 | 0.082 | 60.9% |
| Joint interior displacement/stress fit | 0.003527 | 0.0003232 | 0.503 | 18.9% |

The joint objective sums squared relative L2 residuals of interior displacement
and stress. All metrics are unweighted over sampled nodes/components; they are
not volume integrals or uncertainty-weighted estimates. Stress samples compare
DG0 cell averages with ANSYS nodal values, so local peak errors are sensitive to
mesh and averaging. The 40-resolution surface fit also failed stress validation
(82.2% error); refinement does not remove this problem.

The unconstrained fit corresponds to Cu alpha = 26.11e-6/K and SiO2 alpha =
-1.452e-6/K if T_ref is fixed to 20 C. A surface match alone would disguise the
interior mismatch. The nonnegative fit hits the zero-SiO2-expansion boundary.
Neither is accepted as physical calibration. One temperature snapshot identifies
only alpha*(T-T_ref); alpha and T_ref cannot be separately inferred by this fit.

## Next validation evidence

1. Check the Cu material **assigned to the solved body**, including whether
   bilinear isotropic hardening was active; inspect equivalent plastic strain.
2. Export total strain, elastic strain, temperature, and displacement from the
   same mesh/result time. Their decomposition can distinguish inelastic strain
   from thermal gradients or an inconsistent combined CSV.
3. Confirm the stress-free temperature and CTE values. A 20-to-22 C correction
   alone is too small to resolve this discrepancy.
4. If the target is elastic-only parity, export a verified elastic-only ANSYS
   run with the confirmed geometry and constraints. If the intended target
   includes plasticity, it needs a validated history-dependent constitutive
   model, not fitted elastic constants or output scaling.

## Reproduce

From the repository root in the FEniCSx environment (serial execution):

```sh
python new/calibrate_baseline.py --resolution 80 --output new/calibration_results
python new/validate_baseline.py --check --output /tmp/bonding-baseline-check
```

`calibration.json` contains fit parameters and metrics; `surface_candidates.csv`
contains profile comparisons; `axis_nonelastic_strain.csv` contains the axial
compatibility diagnostic. Existing Stage 15 equations and source CSVs are unchanged.
