"""Postprocess saved mechanical energy inputs; never run a mechanics/PF solve."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import LinearNDInterpolator


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native',type=Path,default=Path('runs/native_100s_retry2'))
    parser.add_argument('--ansys',type=Path,default=Path('runs/ansys_100s_retry'))
    parser.add_argument('--output',type=Path,default=Path('runs/energy_transfer_analysis'))
    args=parser.parse_args()
    a,b=[np.loadtxt(p/'coupling_samples.csv',delimiter=',',skiprows=1) for p in (args.native,args.ansys)]
    np.testing.assert_allclose(a[:,:2],b[:,:2],rtol=0,atol=1e-10)
    r,z=a[:,0],a[:,1];native,ansys=a[:,2],b[:,2]
    interior=(r>1e-7)&(r<250-1e-7)&(z>1e-7)&(z<496.5-1e-7)
    near=interior&(z>486.5)
    def metrics(mask):
        x,y=native[mask],ansys[mask];d=x-y
        return dict(n=int(mask.sum()),native_mean_J_m3=float(x.mean()),ansys_mean_J_m3=float(y.mean()),
                    rms_difference_J_m3=float(np.sqrt(np.mean(d*d))),
                    relative_rms_difference=float(np.linalg.norm(d)/np.linalg.norm(y)) if np.linalg.norm(y)>0 else None,
                    native_zero_nodes=int(np.sum(x==0)),ansys_zero_nodes=int(np.sum(y==0)))
    bands={name:metrics(mask) for name,mask in {
        'Cu_interior':interior,'upper_10nm':near,'upper_10nm_axis_r_lt25':near&(r<25),
        'upper_10nm_middle_r25to225':near&(r>=25)&(r<=225),
        'upper_10nm_edge_r_gt225':near&(r>225)}.items()}
    # Reconstruct the archived ANSYS interpolation, then change only its source mask.
    cfg=json.loads((args.ansys/'config.json').read_text())
    raw=np.loadtxt(cfg['ansys_csv'],delimiter=',',skiprows=1,encoding='utf-8-sig')
    points=np.column_stack((raw[:,1]/500e-9,(raw[:,2]-500e-9)/500e-9))
    surf=np.loadtxt(args.ansys/'mechanical_surface.csv',delimiter=',',skiprows=1)
    top=.993+np.interp(points[:,0],surf[:,0]/500,surf[:,1]/500)
    common=(points[:,0]>=0)&(points[:,0]<=.5)&(points[:,1]>=0)
    legacy=common&(points[:,1]<=top);strict=common&(points[:,1]<=.993)
    def interpolate(mask):
        p,inverse=np.unique(points[mask],axis=0,return_inverse=True)
        w=np.bincount(inverse,weights=np.maximum(raw[mask,13],0))/np.bincount(inverse)
        return np.asarray(LinearNDInterpolator(p,w,fill_value=0)(a[:,:2]/500))
    legacyW,strictW=interpolate(legacy),interpolate(strict)
    check=interior
    reproduction=float(np.max(np.abs(legacyW[check]-ansys[check])))
    if reproduction>1e-3:raise ValueError(f'ANSYS reconstruction does not match saved interior samples: {reproduction}')
    difference=strictW[near]-legacyW[near]
    # Inspect saved initial-interface localization, without treating node sums as FE integrals.
    h=np.interp(r,surf[:,0],surf[:,2])
    eps_nm=cfg['interface_width_nm']
    phi=.5*(1-np.tanh((z-h)/(np.sqrt(2)*eps_nm)))
    localization=16*phi**2*(1-phi)**2
    forcing_delta=15/8*(eps_nm*1e-9)/1.7*localization*(native-ansys)
    report=dict(native=str(args.native),ansys=str(args.ansys),bands=bands,
        ansys_reconstruction_max_interior_difference_J_m3=reproduction,
        outer25nm_fraction_of_upper10nm_squared_error=float(np.sum((native[near&(r>225)]-ansys[near&(r>225)])**2)/np.sum((native[near]-ansys[near])**2)),
        source_mask_extra_points=int(np.sum(legacy&~strict)),
        upper10nm_mask_change_rms_J_m3=float(np.sqrt(np.mean(difference*difference))),
        upper10nm_mask_change_max_J_m3=float(np.max(abs(difference))),
        upper10nm_initial_localized_dimensionless_forcing_difference_rms=float(np.sqrt(np.mean(forcing_delta[near]**2))),
        limitations='Unweighted nodal diagnostics, not FE volume integrals. Forcing uses common ANSYS analytic initial phi. No causal phase-evolution attribution.')
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'energy_transfer.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savetxt(args.output/'upper_surface_samples.csv',np.column_stack((r[near],z[near],native[near],ansys[near],strictW[near])),delimiter=',',header='r_nm,z_phase_nm,native_W_J_m3,ansys_W_J_m3,ansys_strict_mask_W_J_m3',comments='')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
