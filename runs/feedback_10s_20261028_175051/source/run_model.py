#!/usr/bin/env python3
"""Run the model and collect reproducible, shareable diagnostics."""
import argparse
import contextlib
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
import traceback
import zipfile

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'new'))
from model_config import ModelConfig,smoke_config


def write_json(path,data):
    path.write_text(json.dumps(data,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x))+'\n')


def package(folder):
    folder=Path(folder).resolve()
    if not folder.is_dir():raise ValueError(f'No run directory: {folder}')
    target=folder/'analysis.zip'
    files=[p for p in folder.iterdir() if p.is_file() and p.suffix in ('.json','.csv','.log','.md')]
    if not files:raise ValueError('No analysis files to package')
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(files):archive.write(path,path.name)
        snapshot=folder/'source'
        if snapshot.is_dir():
            for path in sorted(snapshot.glob('*.py')):archive.write(path,'source/'+path.name)
    return target


class Tee:
    def __init__(self,*streams):self.streams=streams
    def write(self,text):
        for stream in self.streams:stream.write(text);stream.flush()
        return len(text)
    def flush(self):
        for stream in self.streams:stream.flush()


def environment():
    data={'python':sys.version,'executable':sys.executable,'platform':platform.platform()}
    for name in ('fenics-dolfinx','fenics-ufl','fenics-basix','numpy','scipy','petsc4py','mpi4py'):
        try:data[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:data[name]='not installed'
    result=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,capture_output=True,text=True)
    data['git_commit']=result.stdout.strip() or 'unavailable'
    return data


def run(args):
    data=(smoke_config() if args.preset=='smoke' else ModelConfig()).to_dict()
    if args.config:data.update(json.loads(Path(args.config).read_text()))
    for key in ('duration_s','dt_s','ansys_csv','ansys_interpolation','reference_C'):
        value=getattr(args,key)
        if value is not None:data[key]=value
    if args.source:data['mechanics_source']=args.source
    if args.no_plasticity:data['plasticity']=False
    if args.no_fields:data['write_fields']=False
    config=ModelConfig(**data).validate()
    csv=Path(config.ansys_csv)
    if not csv.is_absolute():config.ansys_csv=str(ROOT/csv)
    out=Path(args.output).resolve()
    if out.exists() and any(out.iterdir()):raise ValueError(f'Refusing to overwrite nonempty directory: {out}')
    out.mkdir(parents=True,exist_ok=True)
    write_json(out/'config.json',config.to_dict());write_json(out/'environment.json',environment())
    snapshot=out/'source';snapshot.mkdir()
    hashes={}
    for path in [ROOT/'run_model.py',*(ROOT/'new'/name for name in ('model_config.py','static_coupling.py','phase_field_solver.py','thermoelastic_solver.py','plasticity.py','validate_baseline.py','calibrate_baseline.py'))]:
        content=path.read_bytes();(snapshot/path.name).write_bytes(content);hashes[str(path.relative_to(ROOT))]=hashlib.sha256(content).hexdigest()
    write_json(out/'source_hashes.json',hashes)
    code=0
    with (out/'run.log').open('w',buffering=1) as log,contextlib.redirect_stdout(Tee(sys.stdout,log)),contextlib.redirect_stderr(Tee(sys.stderr,log)):
        try:
            write_json(out/'status.json',{'status':'running'})
            print(f'Model output: {out}\nMechanics: {config.mechanics_source}; phase duration: {config.duration_s:g} s',flush=True)
            if config.mechanics_source=='native':
                print(f'Thermal reference: {config.reference_C:g} C (assumed input; ANSYS reference has not been verified).',flush=True)
            if args.preset=='smoke':print('SMOKE PRESET: coarse workflow test, not a physical prediction.',flush=True)
            import numpy as np
            from static_coupling import build_static_mechanics
            from phase_field_solver import run_phase_field
            mechanics=build_static_mechanics(config)
            write_json(out/'mechanics_summary.json',mechanics.metadata)
            np.savetxt(out/'mechanical_surface.csv',np.column_stack((mechanics.radii*500,mechanics.uz*500,mechanics.surface_height(mechanics.radii)*500)),delimiter=',',header='reference_r_nm,uz_nm,z_phase_nm',comments='')
            np.savetxt(out/'dielectric_surface.csv',np.column_stack((mechanics.dielectric_radii*500,mechanics.dielectric_uz*500)),delimiter=',',header='reference_r_nm,uz_nm',comments='')
            summary=run_phase_field(config,mechanics,out)
            write_json(out/'summary.json',summary);write_json(out/'status.json',{'status':'completed'})
        except (Exception,KeyboardInterrupt) as error:
            code=1;traceback.print_exc()
            write_json(out/'status.json',{'status':'failed','error':str(error),'type':type(error).__name__})
    print(f'Send this archive for analysis: {package(out)}')
    return code


def compare(args):
    import numpy as np
    left,right=Path(args.left),Path(args.right)
    configs=[json.loads((p/'config.json').read_text()) for p in (left,right)]
    for p in (left,right):
        if json.loads((p/'status.json').read_text())['status']!='completed':raise ValueError(f'Run is not complete: {p}')
    study=args.study
    allowed={'time':{'dt_s','max_dt_s'},'mesh':{'phase_nr','phase_nz'},'energy':{'mechanical_energy_multiplier'}}
    ignored=({'write_fields'}|allowed[study]) if study!='inputs' else {
        'mechanics_source','ansys_csv','ansys_interpolation','plasticity',
        'mechanics_nr','mechanics_nz','thermal_steps','write_fields'}
    differences={k:[configs[0].get(k),configs[1].get(k)] for k in configs[0].keys()|configs[1].keys()
                 if k not in ignored and configs[0].get(k)!=configs[1].get(k)}
    if differences:raise ValueError(f'{study} comparison has incompatible settings: {differences}')
    if study!='inputs' and not any(configs[0][k]!=configs[1][k] for k in allowed[study]):
        raise ValueError(f'{study} study requires different {sorted(allowed[study])}')
    summaries=[json.loads((p/'summary.json').read_text()) for p in (left,right)]
    if any(abs(s['final_time_s']-c['duration_s'])>1e-9 for s,c in zip(summaries,configs)):
        raise ValueError('A run did not reach its configured final time')
    out=Path(args.output).resolve()
    if out.exists() and any(out.iterdir()):raise ValueError(f'Refusing to overwrite: {out}')
    out.mkdir(parents=True,exist_ok=True)
    result={'left':str(left.resolve()),'right':str(right.resolve()),'difference_convention':'left minus right','source_configs':configs,'study':study}
    deltas={}
    for stage in ('initial','final'):
        a=np.loadtxt(left/f'{stage}_surface.csv',delimiter=',',skiprows=1)
        b=np.loadtxt(right/f'{stage}_surface.csv',delimiter=',',skiprows=1)
        if not np.all(np.isfinite(a[:,:2])) or not np.all(np.isfinite(b[:,:2])):raise ValueError('Surface extraction has missing points; analyze individual run diagnostics first')
        if a[0,0]<b[0,0] or a[-1,0]>b[-1,0]:raise ValueError('Surface domains differ')
        other=np.interp(a[:,0],b[:,0],b[:,1]);delta=a[:,1]-other
        deltas[stage]=delta
        result[stage]={'surface_rmse_nm':float(np.sqrt(np.mean(delta**2))),'surface_max_abs_nm':float(np.max(np.abs(delta))),'center_difference_nm':float(delta[0]),'edge_difference_nm':float(delta[-1])}
        np.savetxt(out/f'{stage}_comparison.csv',np.column_stack((a[:,0],a[:,1],other,delta)),delimiter=',',header='r_nm,left_z_phase_nm,right_z_phase_nm,left_minus_right_nm',comments='')
    evolution=deltas['final']-deltas['initial']
    result['evolution_difference']={'surface_change_rmse_nm':float(np.sqrt(np.mean(evolution**2))),'center_change_difference_nm':float(evolution[0]),'edge_change_difference_nm':float(evolution[-1])}
    result['run_summaries']=summaries
    if study!='inputs':
        result['interpretation']='Pairwise sensitivity only; no convergence order or automatic pass/fail claim.'
    write_json(out/'comparison.json',result)
    print(json.dumps({k:result[k] for k in ('initial','final','evolution_difference')},indent=2))
    print(f'Send this archive with both run archives: {package(out)}')
    return 0


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    doctor=sub.add_parser('doctor',help='Check the active Python environment')
    p=sub.add_parser('run',help='Run static mechanics followed by phase-field evolution')
    p.add_argument('--preset',choices=['baseline','smoke'],default='baseline')
    p.add_argument('--config',help='JSON settings; overrides preset')
    p.add_argument('--source',choices=['native','ansys','none'])
    p.add_argument('--duration-s',type=float);p.add_argument('--dt-s',type=float)
    p.add_argument('--reference-c',dest='reference_C',type=float,help='Stress-free reference temperature; default 22 C is an unverified calibration assumption')
    p.add_argument('--ansys-csv');p.add_argument('--ansys-interpolation',choices=['linear','nearest'])
    p.add_argument('--no-plasticity',action='store_true');p.add_argument('--no-fields',action='store_true')
    p.add_argument('--output',required=True)
    p=sub.add_parser('compare',help='Compare matching native and ANSYS-driven phase runs')
    p.add_argument('left');p.add_argument('right');p.add_argument('--output',required=True)
    p.add_argument('--study',choices=['inputs','time','mesh','energy'],default='inputs',
                   help='Allow only the numerical settings varied by this sensitivity study')
    p=sub.add_parser('package',help='Rebuild the small archive, including after an interrupted run')
    p.add_argument('folder')
    args=parser.parse_args()
    try:
        if args.command=='run':return run(args)
        if args.command=='compare':return compare(args)
        if args.command=='package':print(package(args.folder));return 0
        print(json.dumps(environment(),indent=2))
        import dolfinx,numpy,scipy,ufl
        from mpi4py import MPI
        from petsc4py import PETSc
        from dolfinx.fem.petsc import NonlinearProblem
        print(f'Imports OK. DOLFINx {dolfinx.__version__}; MPI ranks {MPI.COMM_WORLD.size}. Run this model with one rank.')
        return 0
    except Exception as error:
        print(f'Error: {error}',file=sys.stderr);return 1


if __name__=='__main__':sys.exit(main())
