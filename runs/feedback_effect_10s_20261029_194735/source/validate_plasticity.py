"""Constitutive and finite-element verification of optional J2 plasticity."""
import json
from dataclasses import replace
from pathlib import Path
import numpy as np
import ufl
from basix.ufl import quadrature_element
from dolfinx import fem
from mpi4py import MPI
from plasticity import BilinearPlasticity, radial_return
from thermoelastic_solver import ReferenceGeometry, cu_at_250C, sio2_dielectric, solve
from validate_baseline import create_baseline_mesh


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('new/calibration_results/plasticity/verification.json'))
    args=parser.parse_args()
    g=ReferenceGeometry();mesh=create_baseline_mesh(g,12,12)
    cu=replace(cu_at_250C(),T_ref_K=295.15);dielectric=replace(sio2_dielectric(),T_ref_K=295.15)
    Q=fem.functionspace(mesh,quadrature_element(mesh.basix_cell(),value_shape=(3,3),degree=4))
    P=fem.functionspace(mesh,quadrature_element(mesh.basix_cell(),degree=4))
    strain=fem.Constant(mesh,np.zeros((3,3)));ep=fem.Constant(mesh,np.zeros((3,3)));p=fem.Constant(mesh,0.)
    C=lambda e: 2*cu.lame_mu_Pa*e+cu.lame_lambda_Pa*ufl.tr(e)*ufl.Identity(3)
    law=BilinearPlasticity();H=law.hardening_modulus(cu.E_Pa)
    sigma,dep,dp=radial_return(strain-ep,C,cu.lame_mu_Pa,p,law.yield_stress_Pa,H,1.)
    se=fem.Expression(sigma,Q.element.interpolation_points);ee=fem.Expression(dep,Q.element.interpolation_points)
    pe=fem.Expression(dp,P.element.interpolation_points)
    sf=fem.Function(Q);ef=fem.Function(Q);pf=fem.Function(P)
    history=[]
    for axial in [.0005,.002,.004,.007]:
        expected_stress=min(cu.E_Pa*axial,law.yield_stress_Pa+law.tangent_modulus_Pa*(axial-law.yield_stress_Pa/cu.E_Pa))
        expected_p=axial-expected_stress/cu.E_Pa
        lateral=-cu.nu*expected_stress/cu.E_Pa-expected_p/2
        strain.value=np.diag([axial,lateral,lateral])
        sf.interpolate(se);ef.interpolate(ee);pf.interpolate(pe)
        actual=sf.x.array.reshape(-1,3,3)[0]
        assert np.allclose(actual,np.diag([expected_stress,0.,0.]),rtol=1e-9,atol=.001)
        ep.value=ep.value+ef.x.array.reshape(-1,3,3)[0]
        p.value=float(p.value)+float(pf.x.array[0])
        assert abs(float(p.value)-expected_p)<1e-12
        assert abs(np.trace(ep.value))<1e-12
        history.append({'axial_strain':axial,'stress_Pa':float(actual[0,0]),'equivalent_plastic_strain':float(p.value)})
    # Stress-free unloading retains the plastic deformation.
    strain.value=ep.value.copy();sf.interpolate(se);pf.interpolate(pe)
    assert np.max(np.abs(sf.x.array))<.001
    assert np.max(np.abs(pf.x.array))<1e-14
    def run(law=None,steps=20,factors=None):
        return solve(None,mesh,cu,523.15,g.L0_m,.5,.993,geometry=g,material_dielectric=dielectric,
                     plasticity_cu=law,load_steps=steps,load_factors=factors)
    elastic=run();high=run(BilinearPlasticity(yield_stress_Pa=1e14),steps=4)
    assert np.allclose(elastic.u_field.x.array,high.u_field.x.array,rtol=1e-7,atol=1e-16)
    assert np.allclose(elastic.W_field.x.array,high.W_field.x.array,rtol=1e-6,atol=1e-3)
    assert high.stats['max_equivalent_plastic_strain']==0.
    plastic20=run(law,20);plastic40=run(law,40)
    assert plastic40.stats['max_equivalent_plastic_strain']>0
    assert all(h['max_yield_violation_Pa']<.01 for h in plastic40.stats['plastic_load_history'])
    ep_values=plastic40.plastic_state.strain.x.array.reshape(-1,3,3)
    assert np.max(np.abs(np.trace(ep_values,axis1=1,axis2=2)))<1e-12
    assert np.allclose(ep_values,np.swapaxes(ep_values,1,2),atol=1e-12)
    assert np.min(plastic40.W_field.x.array)>=-1e-6
    snapshot=plastic40.plastic_state.snapshot()
    plastic40.plastic_state.strain.x.array[:]=9.
    plastic40.plastic_state.equivalent_strain.x.array[:]=9.
    plastic40.plastic_state.restore(snapshot)
    assert np.array_equal(plastic40.plastic_state.strain.x.array,snapshot[0])
    assert np.array_equal(plastic40.plastic_state.equivalent_strain.x.array,snapshot[1])
    cooled=run(law,factors=np.r_[np.linspace(.05,1.,20),np.linspace(.95,0.,20)])
    accumulated=[h['max_equivalent_plastic_strain'] for h in cooled.stats['plastic_load_history']]
    assert np.min(np.diff(accumulated))>=-1e-12
    assert cooled.stats['max_equivalent_plastic_strain']>0
    assert np.max(np.abs(cooled.u_field.x.array))>1e-12
    step_change=abs(plastic40.stats['UY_center_nm']-plastic20.stats['UY_center_nm'])
    assert step_change<.01
    report={'passed':True,'mpi_ranks':MPI.COMM_WORLD.size,'material_point_uniaxial_history':history,
            'checks':['uniaxial bilinear law','stress-free unloading with residual plastic strain',
                      'elastic-limit displacement and energy','yield consistency','isochoric symmetric plastic strain',
                      'nonnegative elastic energy','history snapshot/restore','thermal unloading history'],
            'center_uz_20_steps_nm':plastic20.stats['UY_center_nm'],
            'center_uz_40_steps_nm':plastic40.stats['UY_center_nm'],
            'load_step_change_nm':step_change,'cooling_residual_center_uz_nm':cooled.stats['UY_center_nm']}
    if MPI.COMM_WORLD.rank==0:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report,indent=2))


if __name__=='__main__':main()
