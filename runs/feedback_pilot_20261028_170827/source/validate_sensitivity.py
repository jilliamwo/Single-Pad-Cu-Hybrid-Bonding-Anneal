"""Check study isolation and surface-change metrics without running simulations."""
import contextlib
import io
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from run_model import compare
from model_config import ModelConfig


class SensitivityChecks(unittest.TestCase):
    def test_study_isolation(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp)
            def make(name,**changes):
                path=root/name;path.mkdir()
                config=ModelConfig(duration_s=100).to_dict()|changes
                (path/'config.json').write_text(json.dumps(config))
                (path/'status.json').write_text(json.dumps({'status':'completed'}))
                (path/'summary.json').write_text(json.dumps({'final_time_s':100}))
                for stage,height in [('initial',497.),('final',498.)]:
                    np.savetxt(path/f'{stage}_surface.csv',[[0,height],[250,height]],delimiter=',',header='r_nm,z_phase_nm',comments='')
                return path
            base=make('base');half=make('half',dt_s=.5,max_dt_s=.5)
            fine=make('fine',dt_s=.5,max_dt_s=.5,phase_nr=200,phase_nz=440)
            energyoff=make('energyoff',mechanical_energy_multiplier=0.)
            wrong_physics=make('wrong',dt_s=.5,max_dt_s=.5,reference_C=20.)
            for index,(left,right,study,valid) in enumerate([
                (base,half,'time',True),(half,fine,'mesh',True),
                (base,half,'inputs',False),(base,fine,'time',False),
                (base,fine,'mesh',False),(base,wrong_physics,'time',False),
                (base,base,'time',False),(base,energyoff,'energy',True),
                (half,energyoff,'energy',False)]):
                args=SimpleNamespace(left=str(left),right=str(right),study=study,output=str(root/f'comparison{index}'))
                with contextlib.redirect_stdout(io.StringIO()):
                    if valid:
                        self.assertEqual(compare(args),0)
                        result=json.loads((Path(args.output)/'comparison.json').read_text())
                        self.assertEqual(result['study'],study)
                        self.assertEqual(result['evolution_difference']['surface_change_rmse_nm'],0.)
                    else:
                        with self.assertRaises(ValueError):compare(args)
            (half/'summary.json').write_text(json.dumps({'final_time_s':50}))
            with self.assertRaises(ValueError):
                compare(SimpleNamespace(left=str(base),right=str(half),study='time',output=str(root/'incomplete')))


if __name__=='__main__':unittest.main()
