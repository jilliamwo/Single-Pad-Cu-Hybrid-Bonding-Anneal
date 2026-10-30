"""Small-mesh tests of restartable mechanics, a prerequisite for feedback coupling."""
import argparse
import json
from pathlib import Path
from dataclasses import replace
from unittest.mock import patch
import numpy as np
from plasticity import BilinearPlasticity,NonlinearProblem
from thermoelastic_solver import ReferenceGeometry,cu_at_250C,sio2_dielectric,solve
from validate_baseline import create_baseline_mesh


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('/tmp/model-history-continuation.json'))
    args=parser.parse_args()
    g=ReferenceGeometry();mesh=create_baseline_mesh(g,8,8)
    cu=replace(cu_at_250C(),T_ref_K=295.15);dielectric=replace(sio2_dielectric(),T_ref_K=295.15)
    law=BilinearPlasticity()
    def run(factors,previous=None,mesh_override=None):
        return solve(None,mesh if mesh_override is None else mesh_override,cu,523.15,g.L0_m,.5,.993,
                     geometry=g,material_dielectric=dielectric,plasticity_cu=law,load_factors=factors,
                     initial_plastic_state=None if previous is None else previous.plastic_state,
                     initial_displacement=None if previous is None else previous.u_dimless_field)
    full=run([.25,.5,.75,1.])
    first=run([.25,.5]);saved=first.plastic_state.snapshot();saved_u=first.u_dimless_field.x.array.copy()
    second=run([.75,1.],first)
    def same(a,b):
        np.testing.assert_allclose(a.u_dimless_field.x.array,b.u_dimless_field.x.array,rtol=1e-7,atol=1e-10)
        for x,y in zip(a.plastic_state.snapshot(),b.plastic_state.snapshot()):
            np.testing.assert_allclose(x,y,rtol=1e-6,atol=1e-10)
        np.testing.assert_allclose(a.W_field.x.array,b.W_field.x.array,rtol=1e-5,atol=.01)
    same(full,second)
    assert second.stats['history_continued']
    for x,y in zip(first.plastic_state.snapshot(),saved):np.testing.assert_array_equal(x,y)
    np.testing.assert_array_equal(first.u_dimless_field.x.array,saved_u)
    repeated=run([1.],second);same(second,repeated)
    cooled_full=run([.25,.5,.75,1.,.75,.5,.25,0.])
    cooled_split=run([.75,.5,.25,0.],second);same(cooled_full,cooled_split)
    # Fail after one accepted increment inside the continuation call.
    original=NonlinearProblem.solve;calls=0
    def fail_after_one(problem):
        nonlocal calls
        calls+=1
        if calls==2:raise RuntimeError('Injected mechanics failure')
        return original(problem)
    with patch.object(NonlinearProblem,'solve',fail_after_one):
        try:run([.75,1.],first)
        except RuntimeError as error:assert 'Injected' in str(error)
        else:raise AssertionError('Failure must propagate')
    for x,y in zip(first.plastic_state.snapshot(),saved):np.testing.assert_array_equal(x,y)
    np.testing.assert_array_equal(first.u_dimless_field.x.array,saved_u)
    other_mesh=create_baseline_mesh(g,8,8)
    try:run([1.],first,other_mesh)
    except ValueError as error:assert 'same mesh' in str(error)
    else:raise AssertionError('Different mesh must require explicit history remapping')
    report={'passed':True,'mesh':'8 radial subdivisions; 8 per main vertical layer',
            'checks':['split heating equals continuous path','input history and displacement immutable',
                      'same-load re-equilibration preserves state','split cooling equals continuous path',
                      'failed continuation leaves caller state unchanged','different mesh rejected'],
            'split_vs_full_max_displacement_difference_dimless':float(np.max(abs(second.u_dimless_field.x.array-full.u_dimless_field.x.array))),
            'heated_max_equivalent_plastic_strain':second.stats['max_equivalent_plastic_strain'],
            'scope':'fixed geometry only; not an evolving-phase two-way simulation'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
