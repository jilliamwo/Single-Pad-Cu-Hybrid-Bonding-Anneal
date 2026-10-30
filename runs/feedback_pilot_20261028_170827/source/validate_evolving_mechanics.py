"""Small-mesh state-transition and mechanics-feedback verification."""
import json
from pathlib import Path
from dataclasses import replace
from unittest.mock import patch
import numpy as np
from thermoelastic_solver import ReferenceGeometry,cu_at_250C,sio2_dielectric,solve
from plasticity import BilinearPlasticity
from validate_baseline import create_baseline_mesh
from evolving_mechanics import EvolvingMechanics,occupancy_from_surface,prepare_history


def main():
    g=ReferenceGeometry();mesh=create_baseline_mesh(g,8,8)
    cu=replace(cu_at_250C(),T_ref_K=295.15);d=replace(sio2_dielectric(),T_ref_K=295.15);law=BilinearPlasticity()
    initial=solve(None,mesh,cu,523.15,g.L0_m,.5,.993,geometry=g,material_dielectric=d,plasticity_cu=law,load_steps=4)
    engine=EvolvingMechanics(mesh,g,cu,d,523.15,law,initial)
    radii=initial.UY_r_dimless;surface=.993+initial.UY_dimless
    snapshot=initial.plastic_state.snapshot();u0=initial.u_dimless_field.x.array.copy()
    same=engine.propose(radii,surface)
    assert same.diagnostics['changed_cells']==0
    np.testing.assert_allclose(same.result.u_dimless_field.x.array,u0,atol=1e-10,rtol=1e-7)
    for a,b in zip(same.result.plastic_state.snapshot(),snapshot):np.testing.assert_allclose(a,b,atol=1e-10,rtol=1e-6)
    np.testing.assert_allclose(same.result.W_field.x.array,initial.W_field.x.array,atol=.01,rtol=1e-5)
    # A large synthetic erosion resolves material loss even on this tiny test mesh.
    eroded_surface=surface-.08
    eroded=engine.propose(radii,eroded_surface)
    assert eroded.diagnostics['removed_integration_points']>0
    assert not np.allclose(eroded.result.W_field.x.array,initial.W_field.x.array)
    for a,b in zip(initial.plastic_state.snapshot(),snapshot):np.testing.assert_array_equal(a,b)
    np.testing.assert_array_equal(initial.u_dimless_field.x.array,u0)
    assert engine.version==0 and engine.result is initial
    with patch('evolving_mechanics.solve',side_effect=RuntimeError('Injected failure')):
        try:engine.propose(radii,eroded_surface)
        except RuntimeError:pass
        else:raise AssertionError('Failed mechanics proposal must propagate')
    assert engine.result is initial and engine.version==0
    engine.commit(eroded)
    try:engine.commit(same)
    except ValueError:pass
    else:raise AssertionError('Stale trial must not commit')
    try:engine.propose(radii,surface)
    except ValueError as error:assert 'explicit history policy' in str(error)
    else:raise AssertionError('Growth must not silently invent history')
    # Inspect reset rule directly, before solving heated newly occupied material.
    reset,counts=prepare_history(engine.result.plastic_state,engine.occupancy,same.occupancy,'virgin_at_reference')
    assert counts['gained_integration_points']>0
    assert np.all(reset.equivalent_strain.x.array>=0)
    permissive=EvolvingMechanics(mesh,g,cu,d,523.15,law,initial,new_cu_history='virgin_at_reference')
    growing=permissive.propose(radii,surface+.006)
    assert growing.diagnostics['gained_integration_points']>0
    assert np.all(np.isfinite(growing.result.W_field.x.array))
    assert permissive.result is initial
    report={'passed':True,'checks':['unchanged occupancy reproduces calibrated sharp solve',
        'erosion updates equilibrium and energy','trial proposals do not mutate accepted history',
        'failed mechanics proposal leaves accepted state intact','stale commit rejected',
        'growth rejected by default','explicit virgin-reference growth produces finite mechanics'],
        'erosion':eroded.diagnostics,'growth':growing.diagnostics,
        'scope':'synthetic 8-by-8-per-layer mechanics tests; no coupled phase evolution or physical validation'}
    out=Path('new/verification_results/evolving_mechanics.json');out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
