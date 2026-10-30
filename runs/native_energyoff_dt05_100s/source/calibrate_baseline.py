"""Audit the ANSYS export and fit effective thermal strains, without changing defaults.

Run serially in the FEniCSx environment. Only alpha*(T-T_ref) is identifiable
from one uniform-temperature elastic state, so fitted CTEs are conditional on T_ref.
"""
import argparse
import json
from pathlib import Path
from dataclasses import replace
import numpy as np
from dolfinx import geometry as dg
from validate_baseline import create_baseline_mesh
from thermoelastic_solver import ReferenceGeometry, cu_at_250C, sio2_dielectric, solve


def load_ansys(path):
    """Resolve column names; older temperature exports have a different order."""
    import csv
    columns = [('Node Number',), ('X Location (m)',), ('Y Location (m)',),
               ('Z Location (m)',), ('Directional Deformation Y (m)', 'Y Deformation (m)'),
               ('Normal Stress X (Pa)', 'Stress X (Pa)'), ('Normal Stress Y (Pa)', 'Stress Y (Pa)'),
               ('Normal Stress Z (Pa)', 'Stress Z (Pa)'), ('Shear Stress XY (Pa)',),
               ('Elastic Strain X', 'Normal Elastic Strain X'),
               ('Elastic Strain Y', 'Normal Elastic Strain Y'),
               ('Elastic Strain Z', 'Normal Elastic Strain Z'),
               ('Shear Elastic Strain XY', 'Shear Strain XY'),
               ('W (J/m^3)', 'W Strain Energy Density (J/m^3)')]
    with path.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        names = [next((n for n in aliases if n in reader.fieldnames), None) for aliases in columns]
        if None in names:
            raise ValueError('Missing ANSYS displacement, stress, strain, or W columns')
        data = np.array([[float(row[n]) for n in names] for row in reader])
    if len(data) == 0 or not np.all(np.isfinite(data)):
        raise ValueError('ANSYS data must be nonempty and finite')
    return data


def evaluate(f, points):
    mesh = f.function_space.mesh
    tree = dg.bb_tree(mesh, mesh.topology.dim)
    candidates = dg.compute_collisions_points(tree, points)
    collisions = dg.compute_colliding_cells(mesh, candidates, points)
    cells = np.array([collisions.links(i)[0] if len(collisions.links(i)) else -1
                      for i in range(len(points))], dtype=np.int32)
    if np.any(cells < 0):
        raise ValueError('Some ANSYS points are outside the mechanics mesh')
    return f.eval(points, cells)


def audit(a):
    r, z = a[:,1]*1e9, a[:,2]*1e9
    masks = {'Cu': (r<240)&(z>510)&(z<986),
             'SiO2_base': (r<490)&(z>10)&(z<490),
             'SiO2_annulus': (r>260)&(r<490)&(z>510)&(z<990)}
    fits = {}
    for name, mask in masks.items():
        stress, strain = a[mask,5:8], a[mask,9:12]
        X = np.column_stack((np.repeat(strain.sum(axis=1),3), 2*strain.ravel()))
        lam, mu = np.linalg.lstsq(X, stress.ravel(), rcond=None)[0]
        fits[name] = {'n_nodes': int(mask.sum()),
                      'E_GPa': float(mu*(3*lam+2*mu)/(lam+mu)/1e9),
                      'nu': float(lam/(2*(lam+mu))),
                      'relative_residual': float(np.linalg.norm(X@[lam,mu]-stress.ravel())/np.linalg.norm(stress))}
    s = a[masks['Cu'],5:8]
    vm = np.sqrt(((s[:,0]-s[:,1])**2+(s[:,1]-s[:,2])**2+(s[:,2]-s[:,0])**2)/2
                 +3*a[masks['Cu'],8]**2)/1e6
    # On the axis, duz/dz - epsilon_zz_el = thermal + other eigenstrain.
    line = a[(np.abs(r)<.001)&(z>510)]
    line = line[np.argsort(line[:,2])]
    heights, unique = np.unique(line[:,2], return_index=True)
    line = line[unique]
    midpoint = (heights[1:]+heights[:-1])/2*1e9
    nonelastic = np.diff(line[:,4])/np.diff(heights) - (line[1:,10]+line[:-1,10])/2
    bins = {}
    for low, high in [(525,750),(750,900),(900,950),(970,997)]:
        values = nonelastic[(midpoint>low)&(midpoint<high)]
        bins[f'{low}-{high}_nm'] = {'median':float(np.median(values)),
                                  'min':float(values.min()), 'max':float(values.max())}
    # Engineering vs tensor shear convention inferred from sigma_xy/strain_xy.
    m = masks['Cu'] & (np.abs(a[:,12])>0.)
    denominator = np.dot(a[m,12],a[m,12])
    shear_slope = float(np.dot(a[m,12],a[m,8])/denominator) if denominator > 0 else None
    # Compare the exported W against its apparent reconstruction, not total energy.
    legacy = .5*(a[:,6]*a[:,10]+2*a[:,8]*a[:,12])
    w_error = float(np.linalg.norm(legacy-a[:,13])/np.linalg.norm(a[:,13]))
    full_energy = .5*((a[:,5:8]*a[:,9:12]).sum(axis=1)+a[:,8]*a[:,12])
    full_error = float(np.linalg.norm(full_energy-a[:,13])/np.linalg.norm(a[:,13]))
    return {'elastic_constants':fits, 'Cu_vm_quantiles_MPa':np.quantile(vm,[0,.1,.5,.9,1]).tolist(),
            'Cu_fraction_vm_170_to_180_MPa':float(np.mean((vm>170)&(vm<180))),
            'axis_nonelastic_axial_strain':bins,
            'Cu_shear_stress_over_exported_shear_strain_GPa':shear_slope/1e9 if shear_slope is not None else None,
            'exported_W_legacy_formula_relative_error':w_error,
            'exported_W_full_engineering_shear_formula_relative_error':full_error}, np.column_stack((midpoint,nonelastic)), masks


def metrics(prediction, target):
    error = prediction-target
    return {'rmse':float(np.sqrt(np.mean(error**2))),
            'relative_l2':float(np.linalg.norm(error)/np.linalg.norm(target)),
            'max_abs_error':float(np.max(np.abs(error)))}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ansys',type=Path,default=Path('old/ansys-results/250_500nm_3.5_SiO2.csv'))
    parser.add_argument('--resolution',type=int,default=40)
    parser.add_argument('--output',type=Path,default=Path('/tmp/bonding-calibration'))
    args=parser.parse_args()
    a=load_ansys(args.ansys)
    report,axis,masks=audit(a)
    g=ReferenceGeometry();mesh=create_baseline_mesh(g,args.resolution,args.resolution)
    if mesh.comm.size != 1:
        raise ValueError('Calibration sampling currently requires serial execution')
    m=(a[:,1]<=250e-9+1e-12)&(np.abs(a[:,2]-996.5e-9)<1e-12)
    radii,inv=np.unique(a[m,1],return_inverse=True)
    target_surface=np.bincount(inv,weights=a[m,4])/np.bincount(inv)*1e9
    # Fit alternate surface radii; test the remaining radii and interior fields.
    train=np.arange(len(radii))%2==0
    interior=masks['Cu']|masks['SiO2_base']|masks['SiO2_annulus']
    points=np.zeros((interior.sum(),3));points[:,:2]=a[interior,1:3]/g.L0_m
    def response(theta_cu,theta_dielectric):
        cu=replace(cu_at_250C(),alpha_1_K=theta_cu,T_ref_K=522.15)
        dielectric=replace(sio2_dielectric(),alpha_1_K=theta_dielectric,T_ref_K=522.15)
        result=solve(None,mesh,cu,523.15,g.L0_m,.5,.993,geometry=g,material_dielectric=dielectric)
        surface=np.interp(radii,result.UY_r_dimless*g.L0_m,result.UY_dimless*g.L0_nm)
        uz=evaluate(result.u_field,points)[:,1]*1e9
        stress=np.column_stack([evaluate(result.sigma_components[n],points).ravel()/1e6
                               for n in ['sigma_rr','sigma_zz','sigma_tt','sigma_rz']])
        return surface,uz,stress
    basis=[response(.001,0.),response(0.,.001)]
    surface=np.column_stack([b[0] for b in basis])
    displacement=np.column_stack([b[1] for b in basis])
    stress=np.column_stack([b[2].ravel() for b in basis])
    target_u=a[interior,4]*1e9
    target_s=a[interior,5:9].ravel()/1e6
    coef=np.linalg.lstsq(surface[train],target_surface[train],rcond=None)[0]
    # Exact two-variable nonnegative least squares by checking interior and faces.
    possible=[np.zeros(2)]
    if np.all(coef>=0):
        possible.append(coef)
    for j in range(2):
        c=np.zeros(2); column=surface[train,j]
        c[j]=max(0.,np.dot(column,target_surface[train])/np.dot(column,column))
        possible.append(c)
    nonnegative=min(possible,key=lambda c:np.linalg.norm(surface[train]@c-target_surface[train]))
    # A second fit asks whether one parameter pair can match BOTH interior fields.
    X=np.vstack((displacement/np.linalg.norm(target_u),stress/np.linalg.norm(target_s)))
    y=np.concatenate((target_u/np.linalg.norm(target_u),target_s/np.linalg.norm(target_s)))
    joint=np.linalg.lstsq(X,y,rcond=None)[0]
    report['source']=str(args.ansys)
    report['resolution']=args.resolution
    report['surface_basis_condition_number']=float(np.linalg.cond(surface[train]))
    report['candidates']={}
    candidates={'current_20C':np.array([17e-6*230,.55e-6*230])/.001,
                'reference_22C_hypothesis':np.array([17e-6*228,.55e-6*228])/.001,
                'surface_fit':coef,'nonnegative_surface_fit':nonnegative,'joint_interior_fit':joint}
    for name,c in candidates.items():
        predicted_stress=(stress@c).reshape(-1,4)
        entry={'effective_thermal_strain_Cu':float(c[0]*.001),
               'effective_thermal_strain_SiO2':float(c[1]*.001),
               'Cu_alpha_per_K_if_Tref_20C':float(c[0]*.001/230),
               'SiO2_alpha_per_K_if_Tref_20C':float(c[1]*.001/230),
               'surface_training_nm':metrics((surface@c)[train],target_surface[train]),
               'surface_holdout_nm':metrics((surface@c)[~train],target_surface[~train]),
               'interior_uz_nm':metrics(displacement@c,target_u),
               'interior_stress_MPa':metrics(stress@c,target_s),
               'surface_center_nm':float((surface@c)[0]), 'surface_edge_nm':float((surface@c)[-1])}
        entry['interior_stress_by_material']={region:metrics(predicted_stress[mask[interior]],
               a[interior,5:9][mask[interior]]/1e6) for region,mask in masks.items()}
        report['candidates'][name]=entry
    report['accepted_as_physical_calibration']=False
    report['reason']='Surface-only fits must also match independent interior displacement/stress; alpha and T_ref are not separately identifiable.'
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'calibration.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savetxt(args.output/'axis_nonelastic_strain.csv',axis,delimiter=',',header='z_midpoint_nm,duz_dz_minus_elastic_zz',comments='')
    np.savetxt(args.output/'surface_candidates.csv',np.column_stack((radii*1e9,target_surface,
               *[surface@c for c in candidates.values()])),delimiter=',',
               header='r_nm,ansys_uz_nm,'+','.join(candidates),comments='')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
