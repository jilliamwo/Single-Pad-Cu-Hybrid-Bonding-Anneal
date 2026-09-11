# Calibration across temperatures and materials

Native FEA was run for all 11 exports using each case’s inferred Cu and dielectric
elastic constants. Three promising 250 C cases were refined to 160 radial cells.
Provisional fits are available for SiCOH and Polyimide A/B; no universal
calibration has been established, and solver defaults were not changed.

## Refined candidate fits at 250 C

All CTEs below are in 1e-6/K and assume a shared **22 C stress-free reference**.
They are fitted values, not independently measured material properties.

| Case | Cu CTE | Dielectric CTE | Surface holdout RMSE (nm) | Stress holdout error | Cu W holdout error |
|---|---:|---:|---:|---:|---:|
| PolyimideA | 17.486 | 54.812 | 0.0087 | 3.54% | 4.02% |
| PolyimideB | 17.166 | 22.935 | 0.0036 | 2.81% | 3.16% |
| SiCOH | 17.031 | 11.970 | 0.0041 | 1.95% | 3.22% |

Stress and energy errors decrease with refinement. For Polyimide A, stress
error drops from 7.15% to 3.54%; for Polyimide B, 5.62% to 2.81%; for SiCOH,
3.02% to 1.95%. These cases are compatible with the elastic model to the tested
accuracy. Nodal stress recovery versus cell-average comparisons still limit
interpretation of remaining error.

The estimates obtained from exported strain compatibility, without fitting the
surface, imply dielectric CTEs near 55, 23, and 12 for Polyimide A, Polyimide B,
and SiCOH respectively. Their refined surface RMSEs are 0.0128, 0.0066, and
0.0019 nm. The close agreement supports these parameter scales; it does not
independently identify the reference temperature.

`provisional_material_profiles.json` records the refined joint fits. To use one,
pass separate Cu/dielectric `MaterialParams` with its E, nu, alpha and assumed
T_ref to `solve(..., material_dielectric=...)`. No global Cu default should be
overwritten: the Polyimide files use Cu E=117 GPa, nu=0.30, while SiCOH uses
129 GPa and 0.34. These profiles have only been tested at 250 C.

## All-case joint fit results at resolution 80

| Export | Surface holdout RMSE (nm) | Stress holdout error | Cu W holdout error |
|---|---:|---:|---:|
| 200_500nm_3.5.csv | 0.077 | 9.4% | 21.2% |
| 250C_500nm_3.5recess.csv | 0.316 | 48.1% | 37.0% |
| 250_500nm_3.5_BCB.csv | 0.164 | 10.2% | 13.2% |
| 250_500nm_3.5_PolyimideA.csv | 0.035 | 7.2% | 8.5% |
| 250_500nm_3.5_PolyimideB.csv | 0.014 | 5.6% | 6.4% |
| 250_500nm_3.5_SiCN.csv | 2.161 | 35.9% | 53.0% |
| 250_500nm_3.5_SiCOH.csv | 0.007 | 3.0% | 5.9% |
| 250_500nm_3.5_SiO2.csv | 0.503 | 18.8% | 36.7% |
| 250_500nm_3.5_diamond.csv | 0.004 | 5.9% | 12.0% |
| 300_500nm_3.5.csv | 1.230 | 25.7% | 44.6% |
| 350_500nm_3.5.csv | 2.129 | 31.2% | 48.4% |

The diamond fit is **diagnostic only**. It reproduces the exported 1050 Pa
dielectric modulus; this is not a validated physical diamond material model.
Its residual void stiffness was reduced to 0.00105 Pa so numerical void stiffness
would not compete with the anomalously soft dielectric.

## The two 250 C SiO2 exports are not equivalent

The older `250C_500nm_3.5recess.csv` has center uz=1.8884 nm; the labeled
`250_500nm_3.5_SiO2.csv` has 2.8213 nm. At the outer dielectric side, the older
file has **positive** hoop elastic strain (median 0.001006), while the labeled
file has negative hoop elastic strain (median -0.00012812).

For an elastic dielectric, hoop total strain is ur/r. With positive thermal
expansion, fixed ur=0 requires negative hoop elastic strain. Thus the older
export is incompatible with that same constraint unless other eigenstrains or
result-state inconsistencies are present. A free outer radial boundary is a
plausible diagnostic hypothesis, not a recovered ANSYS setting.

With the free outer radial variant, the older file’s joint-fit stress error drops from 48.1% to 15.0%, and surface RMSE from 0.316 to 0.093 nm.
The remaining error prevents claiming full parity. This variant did not alter
baseline constraints. The older file was excluded from the main temperature-transfer
test; the explicitly labeled SiO2 export was used at 250 C instead.

## Predicting a temperature excluded from calibration

Two shared, nonnegative CTEs were fitted at fixed T_ref=22 C to three SiO2
exports, then used to predict the fourth. The temperature being tested contributes
no displacement/stress data to that CTE fit. Elastic constants remain inferred
from each file and are effectively the same across this series.

| Held-out temperature | Surface RMSE (nm) | Interior stress error | Cu W error |
|---|---:|---:|---:|
| 200 C | 0.118 | 19.2% | 31.1% |
| 250 C | 0.614 | 18.5% | 32.7% |
| 300 C | 1.284 | 26.1% | 47.9% |
| 350 C | 1.971 | 44.5% | 91.9% |

A shared elastic thermal calibration does not reproduce the SiO2 temperature
series. Error increases substantially at higher temperatures, consistent with
an unmodeled nonlinear contribution, but not proof of plasticity. Actual
plastic-strain, temperature, total-strain and load-step exports are still needed
to distinguish material nonlinearity from result/input differences.

## Method and limits

- Effective thermal strains theta_Cu and theta_dielectric are the fitted variables.
  A single-temperature state cannot separately identify CTE and T_ref.
- Linear basis solves provide displacement and stress. A third combined solve
  provides the cross term for quadratic elastic energy. This avoids incorrectly
  scaling W linearly with temperature or fitted amplitude.
- Additional direct-solve tests passed for displacement, stress and quadratic W
  reconstruction in the SiCOH check case at amplitudes 2 and 0.4.
- The joint objective balances relative displacement and stress norms globally.
  Half the interior nodes in each material region are held out with a fixed seed.
  Region-specific stress metrics are saved because global stress norms can hide
  larger errors in softer materials. Surface-fit candidates use alternate radii.
- Spatial holdouts remain from the same simulation and are not independent
  experimental validation. The temperature test excludes entire files.
- Full DG0 energy/stress fields are sampled at interior export locations, away
  from interfaces. Nodal averaging and corner singularities affect local errors.
- Raw exports are preserved. Exact repeated rows are removed for calibration;
  surface radii are deduplicated. The nominal temperature comes from filenames.
- The assumed geometry and perfectly bonded interface model are unchanged.

## Reproduce

```sh
python new/calibrate_sweep.py --resolution 80
python new/calibrate_sweep.py --resolution 160 --case 250_500nm_3.5_SiCOH.csv --case 250_500nm_3.5_PolyimideA.csv --case 250_500nm_3.5_PolyimideB.csv --output new/calibration_results/sweep_refined
python new/calibrate_sweep.py --resolution 80 --case 250C_500nm_3.5recess.csv --legacy-outer-free --output new/calibration_results/legacy_outer_free
python new/calibrate_sweep.py --resolution 12 --case 250_500nm_3.5_SiCOH.csv --verify-basis --output new/calibration_results/sweep_basis_check
python new/summarize_sweep.py
```

JSON files contain all candidate parameters and region-level metrics; each case
also has a CSV of reference and predicted surface profiles.
