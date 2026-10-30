"""Compare elastic and Cu J2 predictions using the material-sheet hardening law.

No yield/CTE fitting to surface displacement: parameters are declared below.
Nominal temperature ramps are quasi-static, uniform and start at 22 C.
"""
import argparse
import json
from dataclasses import replace
from pathlib import Path
import numpy as np
from calibrate_baseline import load_ansys,audit,evaluate,metrics
from thermoelastic_solver import ReferenceGeometry,MechanicalBCs,cu_at_250C,sio2_dielectric,solve
from plasticity import BilinearPlasticity
from validate_baseline import create_baseline_mesh


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resolution',type=int,default=40)
    parser.add_argument('--steps',type=int,default=20)
    parser.add_argument('--case',action='append')
    parser.add_argument('--yield-mpa',type=float,default=172.3)
    parser.add_argument('--tangent-mpa',type=float,default=517.)
    parser.add_argument('--reference-c',type=float,default=22.)
    parser.add_argument('--output',type=Path,default=Path('new/calibration_results/plasticity/sweep'))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    g=ReferenceGeometry();mesh=create_baseline_mesh(g,args.resolution,args.resolution)
    if mesh.comm.size!=1:raise ValueError('CSV field sampling currently requires serial execution')
    law=BilinearPlasticity(args.yield_mpa*1e6,args.tangent_mpa*1e6)
    report={'resolution':args.resolution,'load_steps':args.steps,'reference_C_assumed':args.reference_c,
            'yield_stress_MPa':args.yield_mpa,'tangent_modulus_MPa':args.tangent_mpa,
            'thermal_path':'uniform monotonic quasi-static ramp from reference to nominal temperature',
            'dielectric_plasticity':False,'parameter_fit_to_displacement':False,'cases':{}}
    for path in sorted(Path('old/ansys-results').glob('*.csv')):
        if args.case and path.name not in args.case:continue
        a=np.unique(load_ansys(path),axis=0);audited,_,masks=audit(a)
        ec=audited['elastic_constants']['Cu'];ed=audited['elastic_constants']['SiO2_annulus']
        alpha_d=.55e-6
        for label,value in [('BCB',42.),('PolyimideA',55.),('PolyimideB',23.),('SiCN',3.),('SiCOH',12.),('diamond',1.)]:
            if label in path.name:alpha_d=value*1e-6
        cu=replace(cu_at_250C(),E_Pa=float(f'{ec["E_GPa"]*1e9:.6g}'),nu=round(ec['nu'],5),
                   alpha_1_K=17e-6,T_ref_K=args.reference_c+273.15)
        dielectric=replace(sio2_dielectric(),E_Pa=float(f'{ed["E_GPa"]*1e9:.6g}'),nu=round(ed['nu'],5),
                           alpha_1_K=alpha_d,T_ref_K=args.reference_c+273.15)
        # Older unlabeled export was inconsistent with the fixed radial boundary.
        outer_fixed=path.name!='250C_500nm_3.5recess.csv'
        temp=float(path.name[:3])+273.15
        interior=masks['Cu']|masks['SiO2_base']|masks['SiO2_annulus']
        points=np.zeros((interior.sum(),3));points[:,:2]=a[interior,1:3]/g.L0_m
        surface=(a[:,1]<=250e-9+1e-12)&(np.abs(a[:,2]-996.5e-9)<1e-12)
        radii,inv=np.unique(a[surface,1],return_inverse=True)
        target=np.bincount(inv,weights=a[surface,4])/np.bincount(inv)*1e9
        case={'temperature_C':temp-273.15,'outer_ur_fixed':outer_fixed,'Cu_E_Pa':cu.E_Pa,'Cu_nu':cu.nu,
              'dielectric_E_Pa':dielectric.E_Pa,'dielectric_nu':dielectric.nu,
              'Cu_alpha_1_K':cu.alpha_1_K,'dielectric_alpha_1_K':alpha_d,
              'suspect_diamond_units':'diamond' in path.name,
              'ansys_center_uz_nm':float(target[0]),'predictions':{}}
        profiles=[]
        for mode,plastic in [('elastic',None),('plastic',law)]:
            result=solve(None,mesh,cu,temp,g.L0_m,.5,.993,geometry=g,material_dielectric=dielectric,
                         boundary_conditions=MechanicalBCs(outer_ur_fixed=outer_fixed),
                         residual_void_E_Pa=min(1000.,dielectric.E_Pa*1e-6),plasticity_cu=plastic,load_steps=args.steps)
            uz=evaluate(result.u_field,points)[:,1]*1e9
            stress=np.column_stack([evaluate(result.sigma_components[n],points).ravel()/1e6
                                    for n in ['sigma_rr','sigma_zz','sigma_tt','sigma_rz']])
            W=evaluate(result.W_full_field,points).ravel()
            profile=np.interp(radii,result.UY_r_dimless*g.L0_m,result.UY_dimless*g.L0_nm);profiles.append(profile)
            region_cu=masks['Cu'][interior]
            case['predictions'][mode]={'stats':result.stats,'surface_nm':metrics(profile,target),
                  'interior_uz_nm':metrics(uz,a[interior,4]*1e9),
                  'interior_stress_MPa':metrics(stress,a[interior,5:9]/1e6),
                  'Cu_W_J_m3':metrics(W[region_cu],a[interior,13][region_cu]),
                  'stress_by_region':{name:metrics(stress[mask[interior]],a[interior,5:9][mask[interior]]/1e6) for name,mask in masks.items()}}
        report['cases'][path.name]=case
        np.savetxt(args.output/(path.stem+'_surface.csv'),np.column_stack((radii*1e9,target,*profiles)),delimiter=',',header='r_nm,ansys_uz_nm,elastic_uz_nm,plastic_uz_nm',comments='')
        (args.output/'results.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        print(path.name,'ANSYS',case['ansys_center_uz_nm'],'elastic',case['predictions']['elastic']['stats']['UY_center_nm'],
              'plastic',case['predictions']['plastic']['stats']['UY_center_nm'],flush=True)
    print('Completed elastic/plastic comparison.',flush=True)


if __name__=='__main__':main()
