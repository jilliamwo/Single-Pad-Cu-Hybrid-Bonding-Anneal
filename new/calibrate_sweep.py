"""Case-specific elastic calibration and leave-one-temperature-out validation.

Only effective uniform thermal strains are fitted. No solver defaults are changed.
All inputs use the confirmed single-pad geometry; nominal temperatures come from
filenames. The diamond export is diagnostic only because its inferred E is 1050 Pa.
"""
import argparse
import json
from dataclasses import replace
from pathlib import Path
import numpy as np
from calibrate_baseline import load_ansys, audit, evaluate, metrics
from validate_baseline import create_baseline_mesh
from thermoelastic_solver import ReferenceGeometry, cu_at_250C, sio2_dielectric, solve, MechanicalBCs


def nnls2(X,y):
    c=np.linalg.lstsq(X,y,rcond=None)[0]
    choices=[np.zeros(2)]
    if np.all(c>=0): choices.append(c)
    for j in range(2):
        v=np.zeros(2)
        v[j]=max(0.,np.dot(X[:,j],y)/np.dot(X[:,j],X[:,j]))
        choices.append(v)
    return min(choices,key=lambda v:np.linalg.norm(X@v-y))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resolution',type=int,default=80)
    parser.add_argument('--output',type=Path,default=Path('new/calibration_results/sweep'))
    parser.add_argument('--verify-basis', action='store_true', help='Check stress and quadratic energy against an additional direct solve')
    parser.add_argument('--case', action='append', help='Restrict to named CSVs')
    parser.add_argument('--legacy-outer-free', action='store_true', help='Diagnostic free outer ur for older 250C export only')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    g=ReferenceGeometry();mesh=create_baseline_mesh(g,args.resolution,args.resolution)
    if mesh.comm.size!=1: raise ValueError('Sampling requires serial execution')
    cache={};cases={};outputs={};basis_checks=[]
    for path in sorted(Path('old/ansys-results').glob('*.csv')):
        if args.case and path.name not in args.case: continue
        outer_free=args.legacy_outer_free and path.name=='250C_500nm_3.5recess.csv'
        a=load_ansys(path)
        # Exact repeated rows are redundant; distinct interface rows are retained.
        a=np.unique(a,axis=0)
        audited,axis,masks=audit(a)
        ec=audited['elastic_constants']['Cu'];ed=audited['elastic_constants']['SiO2_annulus']
        # Round export-reconstructed values to reduce sensitivity to CSV rounding.
        Ecu=float(f'{ec["E_GPa"]*1e9:.6g}'); ncu=round(ec['nu'],5)
        Ed=float(f'{ed["E_GPa"]*1e9:.6g}'); nd=round(ed['nu'],5)
        key=(Ecu,ncu,Ed,nd,outer_free)
        temp=float(path.name[:3])
        interior=masks['Cu']|masks['SiO2_base']|masks['SiO2_annulus']
        points=np.zeros((interior.sum(),3));points[:,:2]=a[interior,1:3]/g.L0_m
        m=(a[:,1]<=250e-9+1e-12)&(np.abs(a[:,2]-996.5e-9)<1e-12)
        radii,inv=np.unique(a[m,1],return_inverse=True)
        target_surface=np.bincount(inv,weights=a[m,4])/np.bincount(inv)*1e9
        target_u=a[interior,4]*1e9; target_s=a[interior,5:9]/1e6
        target_W=a[interior,13]
        is_cu=masks['Cu'][interior]
        # Stress and displacement both train on one half of interior nodes.
        # Fixed seed avoids relying on export node-number/row ordering.
        rng=np.random.default_rng(451)
        train=np.zeros(interior.sum(),dtype=bool)
        for region in masks.values():
            idx=np.flatnonzero(region[interior]);idx=rng.permutation(idx)
            train[idx[:len(idx)//2]]=True
        surface_train=np.arange(len(radii))%2==0
        if key not in cache:
            results=[]
            for th_c,th_d in [(.001,0.),(0.,.001),(.001,.001)]:
                cu=replace(cu_at_250C(),E_Pa=Ecu,nu=ncu,alpha_1_K=th_c,T_ref_K=522.15)
                dielectric=replace(sio2_dielectric(),E_Pa=Ed,nu=nd,alpha_1_K=th_d,T_ref_K=522.15)
                # Keep residual void far softer than even the suspect diamond case.
                residual=min(1000.,Ed*1e-6)
                result=solve(None,mesh,cu,523.15,g.L0_m,.5,.993,geometry=g,
                             material_dielectric=dielectric,residual_void_E_Pa=residual,
                             boundary_conditions=MechanicalBCs(outer_ur_fixed=not outer_free))
                results.append(result)
            if args.verify_basis:
                direct=solve(None,mesh,replace(cu,alpha_1_K=.002),523.15,g.L0_m,.5,.993,geometry=g,
                    material_dielectric=replace(dielectric,alpha_1_K=.0004),residual_void_E_Pa=residual,
                    boundary_conditions=MechanicalBCs(outer_ur_fixed=not outer_free))
                assert np.allclose(direct.u_field.x.array,2*results[0].u_field.x.array+.4*results[1].u_field.x.array,rtol=1e-7,atol=1e-16)
                for component in direct.sigma_components:
                    expected=2*results[0].sigma_components[component].x.array+.4*results[1].sigma_components[component].x.array
                    assert np.allclose(direct.sigma_components[component].x.array,expected,rtol=1e-6,atol=.01)
                e0,e1,e2=[b.W_full_field.x.array for b in results]
                expected=4*e0+.16*e1+.8*(e2-e0-e1)
                assert np.allclose(direct.W_full_field.x.array,expected,rtol=1e-6,atol=1e-6)
                basis_checks.append({'case':path.name,'passed':True,'checks':['displacement','stress','quadratic_energy']})
            cache[key]=results
        result_fields=cache[key]
        surface=np.column_stack([np.interp(radii,b.UY_r_dimless*g.L0_m,b.UY_dimless*g.L0_nm)
                                 for b in result_fields[:2]])
        displacement=np.column_stack([evaluate(b.u_field,points)[:,1]*1e9 for b in result_fields[:2]])
        stress=np.stack([np.column_stack([evaluate(b.sigma_components[n],points).ravel()/1e6
                        for n in ['sigma_rr','sigma_zz','sigma_tt','sigma_rz']]) for b in result_fields[:2]],axis=-1)
        energies=np.column_stack([evaluate(b.W_full_field,points).ravel() for b in result_fields])
        def design(mask):
            unorm=np.linalg.norm(target_u[mask]);snorm=np.linalg.norm(target_s[mask])
            return (np.vstack((displacement[mask]/unorm,stress[mask].reshape(-1,2)/snorm)),
                    np.concatenate((target_u[mask]/unorm,target_s[mask].ravel()/snorm)))
        X,y=design(train)
        outer=(a[:,1]>499.99e-9)&(a[:,2]>510e-9)&(a[:,2]<990e-9)
        # On fixed outer r, epsilon_tt_total=0. This inference assumes elastic
        # dielectric and no additional eigenstrain; it is NOT independent data.
        thd=float(np.median(-a[outer,11]))
        if outer_free:
            line=a[outer];line=line[np.argsort(line[:,2])]
            thd=float(np.median(np.diff(line[:,4])/np.diff(line[:,2])-(line[1:,10]+line[:-1,10])/2))
        thc=float(audited['axis_nonelastic_axial_strain']['970-997_nm']['median'])
        candidates={'export_thermal_estimate':np.array([thc,thd])/.001,
                    'surface_fit':np.linalg.lstsq(surface[surface_train],target_surface[surface_train],rcond=None)[0],
                    'surface_nonnegative':nnls2(surface[surface_train],target_surface[surface_train]),
                    'joint_nonnegative':nnls2(X,y)}
        def assess(c):
            pu=displacement@c;ps=stress@c
            W=energies[:,0]*c[0]**2+energies[:,1]*c[1]**2+(energies[:,2]-energies[:,0]-energies[:,1])*c[0]*c[1]
            return {'theta_Cu':float(c[0]*.001),'theta_dielectric':float(c[1]*.001),
                    'alpha_Cu_if_Tref22_per_K':float(c[0]*.001/(temp-22)),
                    'alpha_dielectric_if_Tref22_per_K':float(c[1]*.001/(temp-22)),
                    'center_uz_nm':float((surface@c)[0]),
                    'surface_holdout_nm':metrics((surface@c)[~surface_train],target_surface[~surface_train]),
                    'interior_holdout_uz_nm':metrics(pu[~train],target_u[~train]),
                    'interior_holdout_stress_MPa':metrics(ps[~train],target_s[~train]),
                    'Cu_holdout_W_J_m3':metrics(W[is_cu&~train],target_W[is_cu&~train]),
                    'stress_holdout_by_region':{name:metrics(ps[mask[interior]&~train],target_s[mask[interior]&~train])
                                                for name,mask in masks.items()}}
        entries={name:assess(c) for name,c in candidates.items()}
        # Verify displacement superposition against the direct
        # combined solve; the quadratic W expression is required for fitted loads.
        direct_surface=np.interp(radii,result_fields[2].UY_r_dimless*g.L0_m,result_fields[2].UY_dimless*g.L0_nm)
        assert np.allclose(surface@np.ones(2),direct_surface,rtol=1e-7,atol=1e-8)
        outputs[path.name]={'temperature_C_from_filename':temp,'Cu_E_Pa':Ecu,'Cu_nu':ncu,
                 'dielectric_E_Pa':Ed,'dielectric_nu':nd,
                 'residual_void_E_Pa':min(1000.,Ed*1e-6),
                 'ansys_center_uz_nm':float(target_surface[0]),
                 'suspect_material_units':Ed<1e6,'outer_ur_fixed':not outer_free,
                 'surface_basis_condition_number':float(np.linalg.cond(surface[surface_train])),
                 'joint_basis_condition_number':float(np.linalg.cond(X)),
                 'candidates':entries,'promoted_to_solver_defaults':False}
        # Store arrays for genuine held-out temperature prediction with shared CTEs.
        cases[path.name]={'temp':temp,'X':X,'y':y,'surface':surface,'target_surface':target_surface,
                         'stress':stress,'target_s':target_s,'displacement':displacement,'target_u':target_u,
                         'energies':energies,'target_W':target_W,'is_cu':is_cu}
        np.savetxt(args.output/(path.stem+'_surfaces.csv'),np.column_stack((radii*1e9,target_surface,
                   *[surface@c for c in candidates.values()])),delimiter=',',
                   header='r_nm,ansys_uz_nm,'+','.join(candidates),comments='')
        print(path.name, 'joint stress holdout relative error',entries['joint_nonnegative']['interior_holdout_stress_MPa']['relative_l2'],flush=True)
    names=['200_500nm_3.5.csv','250_500nm_3.5_SiO2.csv','300_500nm_3.5.csv','350_500nm_3.5.csv']
    if not all(n in cases for n in names): names=[]
    transfer={}
    for heldout in names:
        training=[cases[n] for n in names if n!=heldout]
        # Common fixed reference; infer two shared CTEs from other temperatures.
        X=np.vstack([c['X']*((c['temp']-22)/.001) for c in training]);y=np.concatenate([c['y'] for c in training])
        alpha=nnls2(X,y);case=cases[heldout];c=alpha*(case['temp']-22)/.001
        W=case['energies'][:,0]*c[0]**2+case['energies'][:,1]*c[1]**2+(case['energies'][:,2]-case['energies'][:,0]-case['energies'][:,1])*c[0]*c[1]
        transfer[heldout]={'training_files':[n for n in names if n!=heldout],
            'alpha_Cu_per_K':float(alpha[0]),'alpha_dielectric_per_K':float(alpha[1]),'Tref_C_assumed':22.,
            'surface_nm':metrics(case['surface']@c,case['target_surface']),
            'interior_stress_MPa':metrics(case['stress']@c,case['target_s']),
            'interior_uz_nm':metrics(case['displacement']@c,case['target_u']),
            'Cu_W_J_m3':metrics(W[case['is_cu']],case['target_W'][case['is_cu']])}
    report={'resolution':args.resolution,'cases':outputs,'temperature_transfer':transfer,
            'defaults_changed':False,'basis_verification':basis_checks,
            'notes':['Case-specific constants inferred from exports; no invented dielectric CTE tables.',
                     'Only thermal strains are identifiable per case; quoted CTEs assume 22 C reference.',
                     'Spatial holdouts are from the same run; temperature transfer holds out whole files.',
                     'Temperature labels and common geometry assumed; temperature fields unavailable.',
                     'Interior objective balances relative displacement and stress norms globally.',
                     'DG0 cell energy/stress compared to nodal exports; peaks remain mesh sensitive.',
                     'Diamond case is diagnostic only and excluded from physical calibration.']}
    (args.output/'sweep.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('Sweep and held-out temperature calibration complete.',flush=True)


if __name__=='__main__':
    main()
