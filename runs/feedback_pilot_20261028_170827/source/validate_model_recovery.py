"""Fault-injection checks for PETSc recovery and accepted-state preservation."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import numpy as np
from petsc4py import PETSc
from model_config import smoke_config
from static_coupling import build_static_mechanics
from phase_field_solver import run_phase_field,NonlinearProblem


def main():
    config=smoke_config();config.mechanics_source='none';config.write_fields=False
    mechanics=build_static_mechanics(config)
    original=NonlinearProblem.solve
    calls=0
    def fail_second_solve(problem):
        nonlocal calls
        calls+=1
        result=original(problem)
        if calls==2:
            problem.u.x.array[:]=999. # Simulate a corrupted failed iterate.
            raise PETSc.Error(73)
        return result
    with TemporaryDirectory(prefix='model-recovery-') as tmp:
        root=Path(tmp)
        normal=root/'normal';recovered=root/'recovered';failed=root/'failed'
        for folder in (normal,recovered,failed):folder.mkdir()
        baseline=run_phase_field(config,mechanics,normal)
        with patch.object(NonlinearProblem,'solve',fail_second_solve):
            result=run_phase_field(config,mechanics,recovered)
        assert result['accepted_steps']==baseline['accepted_steps']==2
        assert result['solver_rebuilds']==1 and result['rejected_steps']==1
        for name in ('initial_surface.csv','final_surface.csv'):
            np.testing.assert_allclose(np.loadtxt(normal/name,delimiter=',',skiprows=1),
                                       np.loadtxt(recovered/name,delimiter=',',skiprows=1),atol=1e-10,rtol=1e-12)
        assert abs(result['relative_mass_drift'])<1e-10
        errors=json.loads((recovered/'solver_errors.json').read_text())
        assert len(errors)==1 and errors[0]['error_code']==73 and errors[0]['time_s']==.05
        with patch.object(NonlinearProblem,'solve',side_effect=PETSc.Error(73)):
            try:run_phase_field(config,mechanics,failed)
            except RuntimeError as error:assert 'after solver rebuild' in str(error)
            else:raise AssertionError('Repeated PETSc error must terminate')
        errors=json.loads((failed/'solver_errors.json').read_text())
        assert len(errors)==2 and errors[-1]['after_rebuild']
        assert json.loads((failed/'progress.json').read_text())['time_s']==0
        with np.load(failed/'last_accepted_state.npz') as state:
            assert state['time_s']==0 and np.all(np.isfinite(state['mixed_dofs']))
        np.testing.assert_allclose(np.loadtxt(failed/'final_surface.csv',delimiter=',',skiprows=1),
                                   np.loadtxt(normal/'initial_surface.csv',delimiter=',',skiprows=1)[:,:2])
    print('PASS: one-shot solver rebuild, exact rollback, unchanged final surface, bounded repeated failure, saved accepted-state diagnostics')


if __name__=='__main__':main()
