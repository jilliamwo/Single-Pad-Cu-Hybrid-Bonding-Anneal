"""Write a reviewable calibration report and provisional per-case input profiles."""
import json
from pathlib import Path

root=Path(__file__).resolve().parent/'calibration_results'
sweep=json.loads((root/'sweep/sweep.json').read_text())
refined=json.loads((root/'sweep_refined/sweep.json').read_text())
legacy=json.loads((root/'legacy_outer_free/sweep.json').read_text())
checks=json.loads((root/'sweep_basis_check/sweep.json').read_text())
lines=['# Calibration across temperatures and materials', '',
'Native FEA was run for all 11 exports using each case’s inferred Cu and dielectric',
'elastic constants. Three promising 250 C cases were refined to 160 radial cells.',
'Provisional fits are available for SiCOH and Polyimide A/B; no universal',
'calibration has been established, and solver defaults were not changed.', '',
'## Refined candidate fits at 250 C', '',
'All CTEs below are in 1e-6/K and assume a shared **22 C stress-free reference**.',
'They are fitted values, not independently measured material properties.', '',
'| Case | Cu CTE | Dielectric CTE | Surface holdout RMSE (nm) | Stress holdout error | Cu W holdout error |',
'|---|---:|---:|---:|---:|---:|']
profiles={}
for name,d in refined['cases'].items():
 c=d['candidates']['joint_nonnegative']
 label=name.removeprefix('250_500nm_3.5_').removesuffix('.csv')
 lines.append(f"| {label} | {c['alpha_Cu_if_Tref22_per_K']*1e6:.3f} | {c['alpha_dielectric_if_Tref22_per_K']*1e6:.3f} | {c['surface_holdout_nm']['rmse']:.4f} | {c['interior_holdout_stress_MPa']['relative_l2']:.2%} | {c['Cu_holdout_W_J_m3']['relative_l2']:.2%} |")
 profiles[label]={'status':'provisional_single_temperature_joint_fit',
     'source_csv':name,'calibration_temperature_C':250.,'T_ref_K_assumed':295.15,
     'temperature_transfer_validated':False,'mesh_resolution':[160,322],
     'Cu':{'E_Pa':d['Cu_E_Pa'],'nu':d['Cu_nu'],'alpha_1_K':c['alpha_Cu_if_Tref22_per_K']},
     'dielectric':{'E_Pa':d['dielectric_E_Pa'],'nu':d['dielectric_nu'],'alpha_1_K':c['alpha_dielectric_if_Tref22_per_K']},
     'holdout_metrics':{key:c[key] for key in ['surface_holdout_nm','interior_holdout_uz_nm','interior_holdout_stress_MPa','Cu_holdout_W_J_m3']}}
lines+=['','Stress and energy errors decrease with refinement. For Polyimide A, stress',
'error drops from 7.15% to 3.54%; for Polyimide B, 5.62% to 2.81%; for SiCOH,',
'3.02% to 1.95%. These cases are compatible with the elastic model to the tested',
'accuracy. Nodal stress recovery versus cell-average comparisons still limit',
'interpretation of remaining error.', '',
'The estimates obtained from exported strain compatibility, without fitting the',
'surface, imply dielectric CTEs near 55, 23, and 12 for Polyimide A, Polyimide B,',
'and SiCOH respectively. Their refined surface RMSEs are 0.0128, 0.0066, and',
'0.0019 nm. The close agreement supports these parameter scales; it does not',
'independently identify the reference temperature.', '',
'`provisional_material_profiles.json` records the refined joint fits. To use one,',
'pass separate Cu/dielectric `MaterialParams` with its E, nu, alpha and assumed',
'T_ref to `solve(..., material_dielectric=...)`. No global Cu default should be',
'overwritten: the Polyimide files use Cu E=117 GPa, nu=0.30, while SiCOH uses',
'129 GPa and 0.34. These profiles have only been tested at 250 C.', '',
'## All-case joint fit results at resolution 80', '',
'| Export | Surface holdout RMSE (nm) | Stress holdout error | Cu W holdout error |',
'|---|---:|---:|---:|']
for name,d in sweep['cases'].items():
 c=d['candidates']['joint_nonnegative']
 lines.append(f"| {name} | {c['surface_holdout_nm']['rmse']:.3f} | {c['interior_holdout_stress_MPa']['relative_l2']:.1%} | {c['Cu_holdout_W_J_m3']['relative_l2']:.1%} |")
lines+=['', 'The diamond fit is **diagnostic only**. It reproduces the exported 1050 Pa',
'dielectric modulus; this is not a validated physical diamond material model.',
'Its residual void stiffness was reduced to 0.00105 Pa so numerical void stiffness',
'would not compete with the anomalously soft dielectric.', '',
'## The two 250 C SiO2 exports are not equivalent', '',
'The older `250C_500nm_3.5recess.csv` has center uz=1.8884 nm; the labeled',
'`250_500nm_3.5_SiO2.csv` has 2.8213 nm. At the outer dielectric side, the older',
'file has **positive** hoop elastic strain (median 0.001006), while the labeled',
'file has negative hoop elastic strain (median -0.00012812).', '',
'For an elastic dielectric, hoop total strain is ur/r. With positive thermal',
'expansion, fixed ur=0 requires negative hoop elastic strain. Thus the older',
'export is incompatible with that same constraint unless other eigenstrains or',
'result-state inconsistencies are present. A free outer radial boundary is a',
'plausible diagnostic hypothesis, not a recovered ANSYS setting.', '']
c=next(iter(legacy['cases'].values()))['candidates']['joint_nonnegative']
lines += [f"With the free outer radial variant, the older file’s joint-fit stress error drops from 48.1% to {c['interior_holdout_stress_MPa']['relative_l2']:.1%}, and surface RMSE from 0.316 to {c['surface_holdout_nm']['rmse']:.3f} nm.",
'The remaining error prevents claiming full parity. This variant did not alter',
'baseline constraints. The older file was excluded from the main temperature-transfer',
'test; the explicitly labeled SiO2 export was used at 250 C instead.', '',
'## Predicting a temperature excluded from calibration', '',
'Two shared, nonnegative CTEs were fitted at fixed T_ref=22 C to three SiO2',
'exports, then used to predict the fourth. The temperature being tested contributes',
'no displacement/stress data to that CTE fit. Elastic constants remain inferred',
'from each file and are effectively the same across this series.', '',
'| Held-out temperature | Surface RMSE (nm) | Interior stress error | Cu W error |',
'|---|---:|---:|---:|']
for name,c in sweep['temperature_transfer'].items():
 lines.append(f"| {name[:3]} C | {c['surface_nm']['rmse']:.3f} | {c['interior_stress_MPa']['relative_l2']:.1%} | {c['Cu_W_J_m3']['relative_l2']:.1%} |")
lines+=['', 'A shared elastic thermal calibration does not reproduce the SiO2 temperature',
'series. Error increases substantially at higher temperatures, consistent with',
'an unmodeled nonlinear contribution, but not proof of plasticity. Actual',
'plastic-strain, temperature, total-strain and load-step exports are still needed',
'to distinguish material nonlinearity from result/input differences.', '',
'## Method and limits', '',
'- Effective thermal strains theta_Cu and theta_dielectric are the fitted variables.',
'  A single-temperature state cannot separately identify CTE and T_ref.',
'- Linear basis solves provide displacement and stress. A third combined solve',
'  provides the cross term for quadratic elastic energy. This avoids incorrectly',
'  scaling W linearly with temperature or fitted amplitude.',
'- Additional direct-solve tests passed for displacement, stress and quadratic W',
'  reconstruction in the SiCOH check case at amplitudes 2 and 0.4.',
'- The joint objective balances relative displacement and stress norms globally.',
'  Half the interior nodes in each material region are held out with a fixed seed.',
'  Region-specific stress metrics are saved because global stress norms can hide',
'  larger errors in softer materials. Surface-fit candidates use alternate radii.',
'- Spatial holdouts remain from the same simulation and are not independent',
'  experimental validation. The temperature test excludes entire files.',
'- Full DG0 energy/stress fields are sampled at interior export locations, away',
'  from interfaces. Nodal averaging and corner singularities affect local errors.',
'- Raw exports are preserved. Exact repeated rows are removed for calibration;',
'  surface radii are deduplicated. The nominal temperature comes from filenames.',
'- The assumed geometry and perfectly bonded interface model are unchanged.', '',
'## Reproduce', '', '```sh',
'python new/calibrate_sweep.py --resolution 80',
'python new/calibrate_sweep.py --resolution 160 --case 250_500nm_3.5_SiCOH.csv --case 250_500nm_3.5_PolyimideA.csv --case 250_500nm_3.5_PolyimideB.csv --output new/calibration_results/sweep_refined',
'python new/calibrate_sweep.py --resolution 80 --case 250C_500nm_3.5recess.csv --legacy-outer-free --output new/calibration_results/legacy_outer_free',
'python new/calibrate_sweep.py --resolution 12 --case 250_500nm_3.5_SiCOH.csv --verify-basis --output new/calibration_results/sweep_basis_check',
'python new/summarize_sweep.py', '```', '',
'JSON files contain all candidate parameters and region-level metrics; each case',
'also has a CSV of reference and predicted surface profiles.']
assert checks['basis_verification'] and all(c['passed'] for c in checks['basis_verification'])
(root/'SWEEP_REPORT.md').write_text('\n'.join(lines)+'\n')
(root/'provisional_material_profiles.json').write_text(json.dumps(profiles,indent=2,allow_nan=False)+'\n')
print('Wrote SWEEP_REPORT.md and provisional_material_profiles.json')
