# Optional plasticity: implementation and ANSYS comparison

Optional small-strain J2 plasticity is now integrated into `thermoelastic_solver.py`.
The material-sheet Cu law greatly reduces the SiO2 mismatch without fitting yield
stress or hardening to the exported displacement. It does not establish exact
ANSYS parity or improve every material case.

## Model and parameter provenance

- Associated von Mises flow with bilinear isotropic hardening.
- Cu initial yield stress: **172.3 MPa**; uniaxial tangent modulus: **517 MPa**,
  taken from the supplied Material_Sheet.xml. These are trial material-law inputs;
  the ANSYS model assignment still has not been independently inspected.
- Plastic hardening slope is H = E*Et/(E-Et), not Et itself.
- Cu and dielectric E/nu are inferred separately from each export. Most Cu uses
  129 GPa / 0.34; Polyimide exports use 117 GPa / 0.30.
- Cu CTE: 17e-6/K. Dielectric CTEs: SiO2 0.55, BCB 42, Polyimide A 55,
  Polyimide B 23, SiCN 3, SiCOH 12, and diagnostic diamond 1, all in 1e-6/K.
  Non-SiO2 CTE scales are rounded compatibility estimates from the preceding
  export audits, not new fits to displacement.
- The comparison assumes **22 C stress-free temperature**. Existing solver
  material defaults remain unchanged; the driver explicitly supplies this value.
- Temperature increases uniformly from the reference to the filename temperature.
  Each increment solves quasi-static mechanical equilibrium. No inertia, creep,
  transient heat equation, or measured ANSYS temperature history is included.
- Dielectrics remain elastic in these comparisons. The API permits a separate
  J2 dielectric law if justified parameters are supplied, but no dielectric yield
  parameters were invented. In particular, J2 need not be an adequate polymer law.
- The older `250C_500nm_3.5recess.csv` is evaluated with the previously identified
  free-outer-radial-boundary hypothesis. Other cases use the confirmed constraints.
- Diamond is diagnostic only because its exported dielectric E implies 1050 Pa.

## What changed in the solver

Plastic strain and accumulated equivalent plastic strain are stored at degree-4
integration points. A radial-return update enforces the yield surface. The global
Newton tangent is differentiated from the return-mapped stress. History stays
fixed throughout Newton iterations and is committed only after convergence.
Failed increments restore the previous displacement and history before raising.
Snapshots/restoration are exposed on `PlasticState`.

Outputs use elastic strain = total strain - thermal strain - plastic strain.
`W_field` and `W_ansys_like_field` remain **elastic stored energy** to match the
CSV definition. They do not include plastic dissipation or isotropic-hardening
stored energy. The result additionally exposes integration-point `plastic_state`
and a DG0 `equivalent_plastic_strain_field`.

Every plastic solve currently starts from virgin material and integrates its
complete supplied path. Repeated independent calls do not continue prior plastic
history. Coupling to evolving phi is explicitly rejected for now; external
history continuation and phase-field timestep rollback integration are future work.

## Refined 250 C SiO2 result

160 radial x 322 vertical cells, 20 thermal increments:

| Quantity | Elastic | Plastic | ANSYS |
|---|---:|---:|---:|
| Cu-center uz (nm) | 2.173182 | 2.753426 | 2.821300 |
| Surface RMSE (nm) | 0.534960 | 0.046717 | Reference |
| Interior stress relative L2 error | 24.10% | 5.97% | Reference |
| Cu elastic-energy relative L2 error | 49.82% | 6.44% | Reference |

The center error is about 2.4%. The surface-profile RMSE falls by about 91%.
This supports missing inelastic deformation as a major contributor to the prior
mismatch; it is not proof that this exact law was active in ANSYS.

## Mesh and load-step checks

| Radial cells | Thermal increments | Center uz (nm) | Surface RMSE (nm) | Stress error | Cu W error |
|---|---:|---:|---:|---:|---:|
| 40 | 20 | 2.756724 | 0.095323 | 13.17% | 16.37% |
| 80 | 20 | 2.752019 | 0.056766 | 8.34% | 9.38% |
| 80 | 40 | 2.752063 | 0.056840 | 8.32% | 9.36% |
| 160 | 20 | 2.753426 | 0.046717 | 5.97% | 6.44% |

Doubling load increments from 20 to 40 at resolution 80 changes center
displacement by only about 0.000044 nm. Refinement improves stress and energy
agreement, but it does not eliminate the center difference.

The maximum integration-point equivalent plastic strain increases from about
0.019 on the 40-cell mesh to 0.069 on the 160-cell mesh. Local corner peaks
are not converged and can challenge the small-strain approximation locally.
They should not be interpreted as validated failure/damage predictions.

## All 11 exports: same fixed Cu hardening law

40 radial cells and 20 thermal increments. These are direct comparisons with
fixed declared parameters, not per-case optimization of yield or hardening.

| Export | ANSYS center (nm) | Plastic center (nm) | Elastic surface RMSE (nm) | Plastic surface RMSE (nm) | Elastic stress error | Plastic stress error |
|---|---:|---:|---:|---:|---:|---:|
| 200_500nm_3.5.csv | 1.7745 | 1.7714 | 0.084 | 0.034 | 11.9% | 7.9% |
| 250C_500nm_3.5recess.csv | 1.8884 | 1.8897 | 0.063 | 0.009 | 19.6% | 11.9% |
| 250_500nm_3.5_BCB.csv | 13.6940 | 13.6931 | 0.149 | 0.087 | 14.8% | 23.5% |
| 250_500nm_3.5_PolyimideA.csv | 17.2390 | 17.2579 | 0.015 | 0.028 | 14.3% | 18.4% |
| 250_500nm_3.5_PolyimideB.csv | 8.2775 | 8.2808 | 0.005 | 0.006 | 11.2% | 11.1% |
| 250_500nm_3.5_SiCN.csv | 5.2181 | 4.8972 | 1.730 | 0.255 | 72.6% | 26.5% |
| 250_500nm_3.5_SiCOH.csv | 4.3859 | 4.3864 | 0.001 | 0.001 | 5.4% | 5.1% |
| 250_500nm_3.5_SiO2.csv | 2.8213 | 2.7567 | 0.537 | 0.095 | 25.1% | 13.2% |
| 250_500nm_3.5_diamond.csv | 1.8601 | 1.8630 | 0.003 | 0.003 | 10.6% | 10.6% |
| 300_500nm_3.5.csv | 4.1794 | 3.9794 | 1.208 | 0.180 | 41.1% | 20.5% |
| 350_500nm_3.5.csv | 5.6495 | 5.2952 | 1.952 | 0.272 | 58.6% | 28.1% |

SiO2, SiCN, and the older 250 C case improve substantially. At 200 C the
center is reproduced closely; at 300/350 C a residual difference remains.
BCB and Polyimide A have worse stress agreement with this generic Cu law even
where center displacement improves. A good center match alone is insufficient.
These cases need their actual material laws and thermal histories verified.

Metrics use interior nodes away from material interfaces and compare DG0
cell averages to nodal exports. They are sample-weighted errors, not volume
integrals. No held-out parameter fitting is claimed for this fixed-law experiment.

## Verification

Both serial and two-rank MPI plasticity suites passed:

- Analytic uniaxial bilinear stress-strain response, including Et-to-H conversion.
- Stress-free unloading retains residual plastic strain.
- Very high yield stress recovers the existing elastic displacement and W.
- Return-mapped stress satisfies yield consistency.
- Plastic strain is symmetric and trace-free; equivalent plastic strain never decreases.
- Elastic energy is nonnegative; history snapshots restore exactly.
- A heat/cool cycle retains plastic strain and residual deformation.
- Existing elastic baseline checks still pass with plasticity disabled.

All recorded calibration increments converged to the nonlinear residual tolerance.
Individual increment iterations and yield violations are retained in each JSON.

## Run

```sh
python new/validate_plasticity.py
python new/calibrate_plasticity.py --resolution 40 --steps 20
python new/calibrate_plasticity.py --resolution 160 --steps 20 --case 250_500nm_3.5_SiO2.csv --output new/calibration_results/plasticity/baseline160
```

API opt-in:

```python
from plasticity import BilinearPlasticity
result = solve(..., plasticity_cu=BilinearPlasticity(
    yield_stress_Pa=172.3e6, tangent_modulus_Pa=517e6), load_steps=40)
```

For a heat/cool path, pass explicit `load_factors` (fractions of the target
thermal increment) ending at zero. They specify path order, not physical time.
Positive temperature factors start from the separately supplied material
reference temperatures.

## Formulation references

- [ANSYS rate-independent plasticity and bilinear isotropic hardening](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/ans_mat/amp8sq21dldm.html)
- [Elastoplastic return mapping, Jérémy Bleyer](https://bleyerj.github.io/dolfinx_materials/jax_elastoplasticity.html)

Actual ANSYS temperature histories, total/plastic strains, assigned hardening
law, and large-deflection/contact settings are still needed for exact replication.
