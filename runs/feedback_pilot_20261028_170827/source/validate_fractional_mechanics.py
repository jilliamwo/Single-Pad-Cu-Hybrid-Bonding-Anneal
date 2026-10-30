"""Small-mesh verification of experimental constituent-weighted mechanics."""
import json
from pathlib import Path
from dataclasses import replace
import numpy as np
import ufl
from dolfinx import fem
from thermoelastic_solver import ReferenceGeometry,cu_at_250C,sio2_dielectric,solve,axisymmetric_strain
from plasticity import BilinearPlasticity
from validate_baseline import create_baseline_mesh
from evolving_mechanics import EvolvingMechanics,prepare_history


def main():
    g=ReferenceGeometry();mesh=create_baseline_mesh(g,8,8)
    cu=replace(cu_at_250C(),T_ref_K=295.15);d=replace(sio2_dielectric(),T_ref_K=295.15)
    law=BilinearPlasticity()
    initial=solve(None,mesh,cu,523.15,g.L0_m,.5,.993,geometry=g,material_dielectric=d,plasticity_cu=law,load_steps=4)
    engine=EvolvingMechanics(mesh,g,cu,d,523.15,law,initial,occupancy_mode='fractional',new_cu_history='virgin_at_reference')
    radii=initial.UY_r_dimless;surface=.993+initial.UY_dimless
    snapshot=initial.plastic_state.snapshot()
    same=engine.propose(radii,surface)
    np.testing.assert_allclose(same.result.u_dimless_field.x.array,initial.u_dimless_field.x.array,atol=1e-10,rtol=1e-6)
    np.testing.assert_allclose(same.result.W_field.x.array,initial.W_field.x.array,atol=.02,rtol=1e-5)
    eroded=engine.propose(radii,surface-1e-5)
    grown=engine.propose(radii,surface+1e-5)
    for trial in (eroded,grown):
        assert trial.diagnostics['conservative_geometry']['fractional_cells']>0
        assert np.max(np.abs(trial.result.u_dimless_field.x.array-initial.u_dimless_field.x.array))>1e-10
        assert np.max(np.abs(trial.result.W_field.x.array-initial.W_field.x.array))>1.
        assert np.min(trial.result.W_field.x.array)>-1e-6
        assert trial.result.stats['plastic_load_history'][-1]['max_yield_violation_Pa']<1.
        r=ufl.SpatialCoordinate(mesh)[0]
        strain=axisymmetric_strain(trial.result.u_dimless_field,r)-cu.alpha_1_K*(523.15-cu.T_ref_K)*ufl.Identity(3)-trial.result.plastic_state.strain
        stress=2*cu.lame_mu_Pa*strain+cu.lame_lambda_Pa*ufl.tr(strain)*ufl.Identity(3)
        dx=ufl.dx(metadata={'quadrature_degree':4})
        expected=fem.assemble_scalar(fem.form(trial.occupancy*.5*ufl.inner(stress,strain)*r*dx))
        actual=fem.assemble_scalar(fem.form(trial.result.W_field*r*dx))
        np.testing.assert_allclose(actual,expected,rtol=1e-10,atol=1e-8)
    for a,b in zip(snapshot,initial.plastic_state.snapshot()):np.testing.assert_array_equal(a,b)
    # Direct history transition: partial removal retains intrinsic history;
    # doubling fraction dilutes intrinsic history by half, conserving f*p.
    old=fem.Function(engine.occupancy.function_space);new=fem.Function(old.function_space)
    old.x.array[:]=.25;new.x.array[:]=.5
    diluted,_=prepare_history(initial.plastic_state,old,new,'virgin_at_reference')
    np.testing.assert_allclose(diluted.equivalent_strain.x.array,snapshot[1]*.5)
    np.testing.assert_allclose(diluted.strain.x.array,snapshot[0]*.5)
    retained,_=prepare_history(initial.plastic_state,new,old,'reject_growth')
    np.testing.assert_array_equal(retained.equivalent_strain.x.array,snapshot[1])
    new.x.array[:]=0
    cleared,_=prepare_history(initial.plastic_state,old,new,'reject_growth')
    assert np.all(cleared.equivalent_strain.x.array==0)
    report={'passed':True,'checks':['binary limit matches original solve','0.005 nm erosion and growth change equilibrium and energy',
        'intrinsic Cu yield consistency','Cu energy weighted exactly once','nonnegative stored Cu energy','caller history unchanged',
        'partial growth dilutes history; partial loss retains history; total loss clears history'],
        'erosion':eroded.diagnostics,'growth':grown.diagnostics,
        'scope':'8 by 8 per layer; homogeneous mixture approximation, not physical validation'}
    Path('new/verification_results/fractional_mechanics.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
