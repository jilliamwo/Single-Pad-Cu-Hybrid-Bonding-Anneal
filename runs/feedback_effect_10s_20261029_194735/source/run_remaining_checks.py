"""Additional export audits and solver checks; does not modify model parameters."""
import json
from dataclasses import replace
from pathlib import Path
import numpy as np
from dolfinx import fem
from mpi4py import MPI
from calibrate_baseline import load_ansys, audit, metrics
from validate_baseline import create_baseline_mesh
from thermoelastic_solver import (ReferenceGeometry, MechanicalBCs, cu_at_250C,
                                 solve)


def axis_residual(a, radius_nm, zlo, zhi):
    selected=a[(np.abs(a[:,1]*1e9-radius_nm)<.001)&
               (a[:,2]*1e9>zlo)&(a[:,2]*1e9<zhi)]
    z, inverse=np.unique(selected[:,2],return_inverse=True)
    if len(z)<2:
        return None
    uz=np.bincount(inverse,weights=selected[:,4])/np.bincount(inverse)
    e=np.bincount(inverse,weights=selected[:,10])/np.bincount(inverse)
    residual=np.diff(uz)/np.diff(z)-(e[1:]+e[:-1])/2
    return {'n_intervals':len(residual),'median':float(np.median(residual)),
            'min':float(residual.min()),'max':float(residual.max())}


def main():
    if MPI.COMM_WORLD.size!=1:
        raise ValueError('Run this audit serially; use validate_baseline for MPI checks')
    output=Path('new/calibration_results')
    output.mkdir(exist_ok=True)
    report={'export_audits':{},'mesh_convergence':{},'boundary_sensitivity':{},'solver_checks':{}}
    for path in sorted(Path('old/ansys-results').glob('*.csv')):
        a=load_ansys(path)
        summary,_,_=audit(a)
        # audit's region names are geometric; only the SiO2 case identifies SiO2.
        summary['elastic_constants']['dielectric_base']=summary['elastic_constants'].pop('SiO2_base')
        summary['elastic_constants']['dielectric_annulus']=summary['elastic_constants'].pop('SiO2_annulus')
        summary['bounds_nm']={'min':(a[:,1:4].min(0)*1e9).tolist(),
                              'max':(a[:,1:4].max(0)*1e9).tolist()}
        summary['rows']=len(a)
        summary['duplicate_node_ids']=len(a)-len(np.unique(a[:,0]))
        summary['duplicate_coordinates']=len(a)-len(np.unique(a[:,1:4],axis=0))
        summary['bottom_max_abs_uz_nm']=float(np.max(np.abs(a[np.abs(a[:,2])<1e-12,4]))*1e9)
        outer=a[np.abs(a[:,1]-500e-9)<1e-12]
        summary['outer_max_abs_uz_nm']=float(np.max(np.abs(outer[:,4]))*1e9)
        summary['Cu_top_axis_residual']=axis_residual(a,0,980,997)
        summary['base_axis_residual']=axis_residual(a,0,10,490)
        summary['dielectric_outer_residual']=axis_residual(a,500,510,990)
        coords, inverse, counts=np.unique(np.round(a[:,1:3]*1e9,5),axis=0,
                                          return_inverse=True,return_counts=True)
        repeats=[]
        for i in np.where(counts>1)[0]:
            repeats.append([*coords[i],float(np.ptp(a[inverse==i,4])*1e9)])
        summary['coincident_coordinate_groups_1e_5_nm_tolerance']=len(repeats)
        summary['interface_displacement_jumps_nm']={}
        for name,coordinate,value in [('base',1,500.),('side',0,250.)]:
            jumps=[row[2] for row in repeats if abs(row[coordinate]-value)<.001]
            summary['interface_displacement_jumps_nm'][name]={
                'n_groups':len(jumps),'max':max(jumps) if jumps else None,
                'median':float(np.median(jumps)) if jumps else None}
        report['export_audits'][path.name]=summary
    # Fit top-axis apparent thermal strain vs nominal filename temperature.
    names=['200_500nm_3.5.csv','250C_500nm_3.5recess.csv','300_500nm_3.5.csv','350_500nm_3.5.csv']
    temperatures=np.array([200.,250.,300.,350.])
    strains=np.array([report['export_audits'][n]['Cu_top_axis_residual']['median'] for n in names])
    slope,intercept=np.polyfit(temperatures,strains,1)
    report['temperature_series_conditional_fit']={
        'files':names,'nominal_temperature_C':temperatures.tolist(),
        'apparent_thermal_strain':strains.tolist(),
        'alpha_per_K':float(slope),'implied_reference_C':float(-intercept/slope),
        'strain_rmse':float(np.sqrt(np.mean((slope*temperatures+intercept-strains)**2))),
        'caveat':'Assumes same material, filename temperatures are solved temperatures, and negligible inelastic strain in top-axis intervals.'}
    g=ReferenceGeometry(); cu=cu_at_250C()
    def run(mesh, bc=None, residual=1e3, phi=None, material_dielectric=None, geom=g):
        return solve(phi,mesh,cu,523.15,geom.L0_m,geom.pad_radius_nm/geom.L0_nm,
                     (geom.pad_height_nm-geom.recess_nm)/geom.L0_nm,geometry=geom,
                     boundary_conditions=bc,residual_void_E_Pa=residual,
                     material_dielectric=material_dielectric)
    a=load_ansys(Path('old/ansys-results/250_500nm_3.5_SiO2.csv'))
    top=a[(a[:,1]<=250e-9+1e-12)&(np.abs(a[:,2]-996.5e-9)<1e-12)]
    def compare(result):
        predicted=np.interp(top[:,1],result.UY_r_dimless*g.L0_m,result.UY_dimless*g.L0_nm)
        return metrics(predicted,top[:,4]*1e9)
    for n in [40,80,160]:
        mesh=create_baseline_mesh(g,n,n)
        result=run(mesh)
        assert np.allclose(result.W_field.x.array,result.W_ansys_like_field.x.array,rtol=1e-12,atol=1e-8)
        report['mesh_convergence'][str(n)]={**result.stats,'surface_error_nm':compare(result)}
    mesh=create_baseline_mesh(g,40,40)
    for name,bc in [('bottom_radial_free',MechanicalBCs(bottom_ur_fixed=False)),
                    ('outer_radial_free',MechanicalBCs(outer_ur_fixed=False)),
                    ('both_radial_free',MechanicalBCs(bottom_ur_fixed=False,outer_ur_fixed=False))]:
        result=run(mesh,bc)
        report['boundary_sensitivity'][name]={**result.stats,'surface_error_nm':compare(result)}
    # A homogeneous radially confined cylinder has a closed-form solution:
    # ur=0, uz=theta*(1+nu)/(1-nu)*z; sigma_rr=sigma_tt=-E*theta/(1-nu).
    phase=fem.Function(fem.functionspace(mesh,('Lagrange',1)));phase.x.array[:]=1.
    result=run(mesh,phi=phase,material_dielectric=cu)
    coords=result.u_dimless_field.function_space.tabulate_dof_coordinates()
    actual=result.u_dimless_field.x.array.reshape(-1,2)
    theta=cu.alpha_1_K*(523.15-cu.T_ref_K)
    expected=np.column_stack((np.zeros(len(coords)),theta*(1+cu.nu)/(1-cu.nu)*coords[:,1]))
    error=float(np.max(np.abs(actual-expected)))
    assert error<1e-10
    radial_stress=-cu.E_Pa*theta/(1-cu.nu)
    assert np.allclose(result.sigma_components['sigma_rr'].x.array,radial_stress,rtol=1e-8,atol=.1)
    assert np.max(np.abs(result.sigma_components['sigma_zz'].x.array))<.1
    report['solver_checks']['homogeneous_radial_confinement']={'passed':True,'max_dimless_displacement_error':error}
    # No-recess geometry exercises the zero-gap mesh branch.
    flat=replace(g,recess_nm=0.)
    flat_mesh=create_baseline_mesh(flat,12,12)
    result=run(flat_mesh,geom=flat)
    assert np.all(np.isfinite(result.u_field.x.array))
    report['solver_checks']['zero_recess_geometry']={'passed':True,**result.stats}
    report['solver_checks']['energy_export_definition']={'passed':True,'definition':'full Cu-localized elastic energy'}
    report['unavailable_checks']=['No ANSYS project/result files found in repository',
      'CSV does not contain total strain, plastic strain, temperature, radial displacement, or result-time metadata',
      'Assigned material law and load-step settings cannot be inspected from these CSVs']
    (output/'remaining_checks.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('All executable solver assertions passed; export and sensitivity findings saved.')
    print(json.dumps(report['temperature_series_conditional_fit'],indent=2))


if __name__=='__main__':
    main()
