"""Audit existing native energy samples against one fixed ANSYS interpolant."""
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import LinearNDInterpolator


def main():
    ansys_run=Path('runs/ansys_100s_retry')
    config=json.loads((ansys_run/'config.json').read_text())
    raw=np.loadtxt(config['ansys_csv'],delimiter=',',skiprows=1,encoding='utf-8-sig')
    surface=np.loadtxt(ansys_run/'mechanical_surface.csv',delimiter=',',skiprows=1)
    points=np.column_stack((raw[:,1]/500e-9,(raw[:,2]-500e-9)/500e-9))
    mask=(points[:,0]>=0)&(points[:,0]<=.5)&(points[:,1]>=0)&(points[:,1]<=.993+np.interp(points[:,0],surface[:,0]/500,surface[:,1]/500))
    points,inverse=np.unique(points[mask],axis=0,return_inverse=True)
    values=np.bincount(inverse,weights=raw[mask,13])/np.bincount(inverse)
    interp=LinearNDInterpolator(points,values,fill_value=0)
    output=Path('runs/corner_sampling_analysis');output.mkdir(exist_ok=True)
    rows=[];reference=None
    for name in ['native_dt025_100s','native_meshfine_dt025_100s','native_meshfiner_dt025_100s']:
        run=Path('runs')/name
        cfg=json.loads((run/'config.json').read_text())
        a=np.loadtxt(run/'coupling_samples.csv',delimiter=',',skiprows=1)
        r,z=a[:,0],a[:,1]
        near=(r>1e-7)&(r<250-1e-7)&(z>486.5+1e-7)&(z<496.5-1e-7)
        edge=near&(r>225)
        expected=interp(a[:,:2]/500);delta=a[:,2]-expected
        # Common coordinates distinguish sample location changes from value changes.
        lookup={tuple(x):w for x,w in zip(np.round(a[:,:2],7),a[:,2])}
        if reference is None:reference={tuple(x):w for x,w in zip(np.round(a[near,:2],7),a[near,2])}
        common=np.array([lookup[k]-v for k,v in reference.items()])
        row=dict(run=name,phase_nr=cfg['phase_nr'],upper_nodes=int(near.sum()),
                 upper_relative_rms=float(np.linalg.norm(delta[near])/np.linalg.norm(expected[near])),
                 upper_absolute_rms_J_m3=float(np.sqrt(np.mean(delta[near]**2))),
                 outer25nm_squared_error_fraction=float(np.sum(delta[edge]**2)/np.sum(delta[near]**2)),
                 common_node_native_rms_change_J_m3=float(np.sqrt(np.mean(common**2))),
                 common_node_native_max_change_J_m3=float(np.max(abs(common))))
        rows.append(row)
        np.savetxt(output/f'{name}_corner.csv',np.column_stack((r[edge],z[edge],a[edge,2],expected[edge],delta[edge])),delimiter=',',header='r_nm,z_phase_nm,native_W_J_m3,ansys_interpolated_W_J_m3,difference_J_m3',comments='')
    (output/'sampling.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
