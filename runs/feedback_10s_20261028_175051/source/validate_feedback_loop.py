"""Tiny end-to-end tests of the experimental phase/mechanics feedback loop."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from dataclasses import replace
from unittest.mock import patch
import numpy as np
from model_config import smoke_config
from static_coupling import build_static_mechanics
from thermoelastic_solver import ReferenceGeometry,cu_at_250C,sio2_dielectric
from plasticity import BilinearPlasticity
from evolving_mechanics import EvolvingMechanics
from feedback_coupling import MechanicalFeedback
from phase_field_solver import run_phase_field
from dolfinx import mesh as dm,fem
from mpi4py import MPI


def main(occupancy_mode="binary",transfer_mode="contour"):
    config=smoke_config();config.write_fields=False
    mechanics=build_static_mechanics(config);initial=mechanics.native_result
    g=ReferenceGeometry();cu=replace(cu_at_250C(),T_ref_K=config.reference_C+273.15)
    dielectric=replace(sio2_dielectric(),T_ref_K=config.reference_C+273.15)
    def adapter(policy=None,volume_limit=.01,mode=None):
        policy=policy or ('virgin_at_reference' if occupancy_mode=='fractional' else 'reject_growth')
        engine=EvolvingMechanics(initial.u_dimless_field.function_space.mesh,g,cu,dielectric,
            config.temperature_C+273.15,BilinearPlasticity(),initial,new_cu_history=policy,occupancy_mode=occupancy_mode)
        return MechanicalFeedback(engine,interval_steps=1,max_volume_change=volume_limit,transfer_mode=mode or transfer_mode)
    # A synthetic resolved growth tests that energy transfer does not retain
    # the static sampler's old Cu-top cutoff. This is not a phase prediction.
    phase=dm.create_rectangle(MPI.COMM_WORLD,[np.array([0.,0.]),np.array([.5,1.1])],[20,44])
    V=fem.functionspace(phase,('Lagrange',1));radii=np.linspace(0,.5,201)
    height=mechanics.surface_height(radii)
    growing=adapter('virgin_at_reference',mode='contour')
    energy=growing.initialize(V,radii,height,1.)
    old_energy=energy.x.array.copy()
    growth=growing.propose(height+.006,.1,1)
    assert growth.mechanical.diagnostics['changed_cells']>0
    assert np.any((old_energy==0)&(growth.energy>0))
    assert growing.engine.version==0
    growing.commit(growth,energy)
    assert growing.engine.version==1
    np.testing.assert_array_equal(energy.x.array,growth.energy)
    limited=adapter('virgin_at_reference',.0001,mode='contour')
    limited.initialize(V,radii,height,1.)
    try:limited.propose(height+.006,.1,1)
    except ValueError as error:assert 'volume drift' in str(error)
    else:raise AssertionError('Excessive volume change must not commit')
    assert limited.engine.version==0
    with TemporaryDirectory(prefix='feedback-loop-') as tmp:
        root=Path(tmp)
        normal=root/'normal';retry=root/'retry';failed=root/'failed'
        for folder in (normal,retry,failed):folder.mkdir()
        finer=replace(config,dt_s=.025,max_dt_s=.025)
        clean=adapter();baseline=run_phase_field(finer,mechanics,normal,feedback=clean)
        assert baseline['accepted_steps']==4 and baseline['mechanical_updates']==4
        if transfer_mode=='column_mass':
            assert max(abs(record['cumulative_reference_volume_change']) for record in clean.records)<1e-10
        checkpoint=initial.plastic_state.snapshot()
        recovered=adapter();original=recovered.propose;calls=0
        def interrupted(surface,time_s,step,**kwargs):
            nonlocal calls
            calls+=1
            trial=original(surface,time_s,step,**kwargs)
            if calls==1:raise RuntimeError('Injected failure after mechanical trial, before commit')
            return trial
        with patch.object(recovered,'propose',interrupted):
            result=run_phase_field(config,mechanics,retry,feedback=recovered)
        assert result['rejected_steps']==1 and result['accepted_steps']==4
        assert recovered.engine.version==4 and len(recovered.records)==4
        np.testing.assert_allclose(np.loadtxt(normal/'final_surface.csv',delimiter=',',skiprows=1),
                                   np.loadtxt(retry/'final_surface.csv',delimiter=',',skiprows=1),rtol=0,atol=1e-9)
        np.testing.assert_allclose(np.loadtxt(normal/'latest_energy_samples.csv',delimiter=',',skiprows=1),
                                   np.loadtxt(retry/'latest_energy_samples.csv',delimiter=',',skiprows=1),rtol=1e-7,atol=.01)
        for actual,expected in zip(initial.plastic_state.snapshot(),checkpoint):np.testing.assert_array_equal(actual,expected)
        rejected=adapter()
        with patch.object(rejected,'propose',side_effect=ValueError('Injected domain violation')):
            try:run_phase_field(config,mechanics,failed,feedback=rejected)
            except RuntimeError as error:assert 'policy rejected' in str(error)
            else:raise AssertionError('Invalid geometry must stop')
        assert rejected.engine.version==0 and not rejected.records
        assert json.loads((failed/'progress.json').read_text())['time_s']==0
        np.testing.assert_allclose(np.loadtxt(failed/'final_surface.csv',delimiter=',',skiprows=1),
                                   np.loadtxt(failed/'initial_surface.csv',delimiter=',',skiprows=1)[:,:2])
        report={'passed':True,'occupancy_mode':occupancy_mode,'transfer_mode':transfer_mode,'checks':['real phase/mechanics updates across four accepted steps',
            'post-mechanical-trial failure rolls back phase/history/energy','retry agrees with uninterrupted smaller-step run',
            'input mechanics unchanged','domain rejection stops with accepted state intact',
            'resolved synthetic growth transfers nonzero energy into newly occupied region',
            'volume-limit rejection leaves accepted mechanics untouched'],
            'normal_summary':baseline,'retry_summary':result,'updates':clean.records,
            'scope':'tiny workflow validation only; not physical validation or converged two-way results'}
    path=Path('new/verification_results/feedback_loop'+('_fractional' if occupancy_mode=='fractional' else '')+('_column_mass' if transfer_mode=='column_mass' else '')+'.json');path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--occupancy-mode',choices=['binary','fractional'],default='binary')
    parser.add_argument('--transfer-mode',choices=['contour','column_mass'],default='contour')
    args=parser.parse_args()
    main(args.occupancy_mode,args.transfer_mode)
