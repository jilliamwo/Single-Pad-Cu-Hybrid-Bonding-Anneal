#!/usr/bin/env python3
"""Short experimental two-way update-interval study; serial execution only."""
import argparse
import contextlib
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback
import zipfile
from run_model import ROOT,Tee,environment,write_json
from model_config import ModelConfig,smoke_config


def archive(out):
    target=out/'analysis.zip'
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED) as z:
        for path in sorted(out.rglob('*')):
            if path.is_file() and path.suffix in ('.json','.csv','.log','.md','.py'):
                z.write(path,path.relative_to(out))
    return target


def run(args):
    import numpy as np
    from static_coupling import build_static_mechanics
    from thermoelastic_solver import ReferenceGeometry,cu_at_250C,sio2_dielectric
    from plasticity import BilinearPlasticity
    from evolving_mechanics import EvolvingMechanics
    from feedback_coupling import MechanicalFeedback
    from phase_field_solver import run_phase_field
    # Fixed time steps isolate the update interval. Failure stops the pilot;
    # it cannot silently change one case's time discretization.
    config=(replace(smoke_config(),dt_s=.025,min_dt_s=.025,max_dt_s=.025) if args.smoke else
        replace(ModelConfig(),duration_s=2.,dt_s=.25,min_dt_s=.25,max_dt_s=.25,
                measure_interval_s=.25,save_interval_s=2.,write_fields=False))
    if args.duration_s is not None:
        duration=args.duration_s
        if not np.isfinite(duration) or duration<=0:
            raise ValueError('Duration must be finite and positive')
        steps=duration/config.dt_s
        if not np.isclose(steps,round(steps),rtol=0,atol=1e-9) or round(steps)%2:
            raise ValueError('Duration must span an even number of phase steps so both cases update at the endpoint')
        config=replace(config,duration_s=duration,save_interval_s=duration)
    config.validate()
    out=Path(args.output).resolve()
    if out.exists() and any(out.iterdir()):raise ValueError(f'Refusing to overwrite nonempty directory: {out}')
    out.mkdir(parents=True,exist_ok=True)
    write_json(out/'config.json',config.to_dict());write_json(out/'environment.json',environment())
    write_json(out/'study.json',{'experimental':True,'smoke':args.smoke,'update_intervals_steps':[1,2],
        'occupancy':'fractional','transfer':'column_mass','new_Cu_history':'virgin_at_reference',
        'initial_reference_C':config.reference_C,'ANSYS_reference_verified':False,
        'interpretation':'Pairwise interval sensitivity; not physical validation or convergence proof'})
    source=out/'source';source.mkdir();hashes={}
    for path in [ROOT/'run_feedback_pilot.py',ROOT/'run_model.py',*sorted((ROOT/'new').glob('*.py'))]:
        content=path.read_bytes();(source/path.name).write_bytes(content)
        hashes[str(path.relative_to(ROOT))]=hashlib.sha256(content).hexdigest()
    write_json(out/'source_hashes.json',hashes)
    code=0
    with (out/'run.log').open('w',buffering=1) as log,contextlib.redirect_stdout(Tee(sys.stdout,log)),contextlib.redirect_stderr(Tee(sys.stderr,log)):
        try:
            write_json(out/'status.json',{'status':'running'})
            print(f'EXPERIMENTAL feedback pilot: {out}',flush=True)
            print(f'Phase mesh {config.phase_nr}x{config.phase_nz}; mechanics {config.mechanics_nr}x{config.mechanics_nz} per layer; duration {config.duration_s} s.',flush=True)
            print('Thermal reference 22 C remains unverified. Fractional history uses homogeneous virgin-Cu mixing.',flush=True)
            start=time.perf_counter();mechanics=build_static_mechanics(config)
            write_json(out/'mechanics_summary.json',mechanics.metadata)
            g=ReferenceGeometry();cu=replace(cu_at_250C(),T_ref_K=config.reference_C+273.15)
            dielectric=replace(sio2_dielectric(),T_ref_K=config.reference_C+273.15)
            summaries=[];energies=[];profiles=[];states=[]
            for interval in (1,2):
                folder=out/f'interval_{interval}';folder.mkdir()
                write_json(folder/'status.json',{'status':'running'})
                initial=mechanics.native_result
                engine=EvolvingMechanics(initial.u_dimless_field.function_space.mesh,g,cu,dielectric,
                    config.temperature_C+273.15,BilinearPlasticity(),initial,
                    new_cu_history='virgin_at_reference',occupancy_mode='fractional')
                adapter=MechanicalFeedback(engine,interval_steps=interval,transfer_mode='column_mass')
                print(f'Running mechanics update every {interval} phase step(s), nominal {interval*config.dt_s:g} s.',flush=True)
                begun=time.perf_counter()
                summary=run_phase_field(config,mechanics,folder,feedback=adapter)
                expected=round(config.duration_s/config.dt_s)
                if summary['accepted_steps']!=expected or summary['rejected_steps'] or abs(summary['final_time_s']-config.duration_s)>1e-10:
                    raise RuntimeError('Pilot did not complete the required fixed-step path')
                if summary['mechanical_updates']!=expected//interval or adapter.records[-1]['time_s']!=summary['final_time_s']:
                    raise RuntimeError('Pilot mechanics did not reach the final phase time')
                summary['wall_seconds']=time.perf_counter()-begun
                summary['max_relative_mechanics_volume_drift']=max(abs(r['cumulative_reference_volume_change']) for r in adapter.records)
                summary['max_mass_vs_contour_height_difference_nm']=500*max(r['max_mass_vs_contour_height_difference'] for r in adapter.records)
                write_json(folder/'summary.json',summary)
                write_json(folder/'final_mechanics_stats.json',engine.result.stats)
                write_json(folder/'status.json',{'status':'completed'})
                summaries.append(summary)
                profiles.append(np.loadtxt(folder/'final_surface.csv',delimiter=',',skiprows=1))
                energies.append(np.loadtxt(folder/'latest_energy_samples.csv',delimiter=',',skiprows=1))
                states.append(engine.result.plastic_state.equivalent_strain.x.array.copy())
            np.testing.assert_allclose(profiles[0][:,0],profiles[1][:,0],rtol=0,atol=0)
            np.testing.assert_allclose(energies[0][:,:2],energies[1][:,:2],rtol=0,atol=0)
            delta=profiles[0][:,1]-profiles[1][:,1]
            de=energies[0][:,-1]-energies[1][:,-1]
            report={'difference_convention':'interval 1 minus interval 2','wall_seconds':time.perf_counter()-start,
                'surface_rms_difference_nm':float(np.sqrt(np.mean(delta**2))),
                'surface_max_difference_nm':float(np.max(np.abs(delta))),
                'center_difference_nm':float(delta[0]),'edge_difference_nm':float(delta[-1]),
                'nodal_energy_rms_difference_J_m3':float(np.sqrt(np.mean(de**2))),
                'nodal_energy_relative_L2_difference':float(np.linalg.norm(de)/max(np.linalg.norm(energies[0][:,-1]),1e-30)),
                'max_equivalent_plastic_strain_difference':float(np.max(np.abs(states[0]-states[1]))),
                'run_summaries':summaries,
                'interpretation':'Two update intervals on the same fixed-step phase path; no automatic accuracy/pass threshold'}
            np.savetxt(out/'surface_comparison.csv',np.column_stack((profiles[0][:,0],profiles[0][:,1],profiles[1][:,1],delta)),delimiter=',',header='r_nm,interval1_z_nm,interval2_z_nm,difference_nm',comments='')
            write_json(out/'comparison.json',report);write_json(out/'status.json',{'status':'completed'})
            print(json.dumps(report,indent=2),flush=True)
        except (Exception,KeyboardInterrupt) as error:
            code=1;traceback.print_exc()
            failure={'status':'failed','error':str(error),'type':type(error).__name__}
            write_json(out/'status.json',failure)
            for folder in out.glob('interval_*'):
                status=folder/'status.json'
                if status.exists() and json.loads(status.read_text())['status']=='running':write_json(status,failure)
    print(f'Send this archive for analysis: {archive(out)}',flush=True)
    return code


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    parser.add_argument('--duration-s',type=float,help='Study duration; routine preset requires positive multiples of 0.5 s (default 2 s)')
    parser.add_argument('--smoke',action='store_true',help='Tiny workflow test instead of routine-mesh pilot')
    args=parser.parse_args()
    try:return run(args)
    except Exception as error:print(f'Error: {error}',file=sys.stderr);return 1

if __name__=='__main__':raise SystemExit(main())
