"""Reuse saved phase fields to compare contour and mass-equivalent surface changes.

Run with a Python containing numpy/scipy/h5py/matplotlib. No simulation is run.
"""
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import h5py
from column_mass_transfer import NodalColumnMassTransfer
from conservative_occupancy import graph_volume


def read_saved_phase(folder):
    root=ET.parse(folder/'phi.xdmf').getroot()
    def read(item):
        filename,dataset=item.text.strip().split(':',1)
        with h5py.File(folder/filename,'r') as f:return f[dataset][:]
    grid=root.find("./Domain/Grid[@GridType='Uniform']")
    coords=read(grid.find('./Geometry/DataItem'));cells=read(grid.find('./Topology/DataItem'))
    snapshots=[]
    for grid in root.findall("./Domain/Grid[@GridType='Collection']/Grid"):
        time=float(grid.find('Time').attrib['Value'])
        snapshots.append((time,read(grid.find('./Attribute/DataItem')).ravel()))
    snapshots.sort(key=lambda item:item[0])
    if len(snapshots)<2 or snapshots[0][0]!=0:raise ValueError('Require initial and later saved phase fields')
    return coords,cells,snapshots


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,default=Path('runs/native_dt05_100s'))
    parser.add_argument('--output',type=Path,default=Path('runs/column_mass_audit'))
    args=parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('Output directory must be empty')
    coords,cells,snapshots=read_saved_phase(args.run)
    initial=snapshots[0][1];final=snapshots[-1][1]
    summary=json.loads((args.run/'summary.json').read_text())
    if not np.isclose(snapshots[-1][0],summary['final_time_s']):raise ValueError('Saved field is not the final simulation state')
    a=np.loadtxt(args.run/'initial_surface.csv',delimiter=',',skiprows=1)
    b=np.loadtxt(args.run/'final_surface.csv',delimiter=',',skiprows=1)
    np.testing.assert_allclose(a[:,0],b[:,0],atol=1e-10,rtol=0)
    query=np.linspace(0,.5,2001);contour=np.interp(query,a[:,0]/500,(b[:,1]-a[:,1])/500)
    # Independent exact P1*r triangle integral, not the transfer matrix.
    t=coords[cells];v=t[:,1]-t[:,0];w=t[:,2]-t[:,0]
    area=.5*np.abs(v[:,0]*w[:,1]-v[:,1]*w[:,0])
    def mass(values):
        f=values[cells];r=t[:,:,0]
        return float(np.sum(area/12*(f.sum(axis=1)*r.sum(axis=1)+(f*r).sum(axis=1))))
    phase_change=mass(final-initial);base=.5*.5**2*.993
    reports=[];profiles=[]
    for n in (51,101,201,401):
        radii=np.linspace(0,.5,n)
        transfer=NodalColumnMassTransfer(coords,cells,initial,radii)
        dh,details=transfer.project(final);height=np.interp(query,radii,dh)
        np.testing.assert_allclose(details['phase_mass_change_without_2pi'],phase_change,atol=1e-15,rtol=1e-9)
        error=(height-contour)*500
        reports.append({'radial_nodes':n,**details,
            'mass_balance_error_relative_to_initial_Cu_volume':(graph_volume(radii,dh,0)-phase_change)/base,
            'mass_vs_contour_rms_nm':float(np.sqrt(np.mean(error**2))),
            'mass_vs_contour_max_nm':float(np.max(np.abs(error))),
            'center_mass_height_change_nm':float(dh[0]*500),'edge_mass_height_change_nm':float(dh[-1]*500),
            'reference_height_in_domain':bool(np.all(1.993+dh>1)&np.all(1.993+dh<=2))})
        profiles.append(height)
        print(f'Completed {n} radial nodes',flush=True)
    for i,row in enumerate(reports):
        delta=(profiles[i]-profiles[-1])*500
        row['rms_vs_401_nodes_nm']=float(np.sqrt(np.mean(delta**2)))
        row['max_vs_401_nodes_nm']=float(np.max(np.abs(delta)))
    report={'source_run':str(args.run),'saved_times_s':[x[0] for x in snapshots],
        'phase_mesh_cells':len(cells),'no_phase_or_mechanics_solve':True,
        'independent_phase_mass_change_without_2pi':phase_change,
        'relative_phase_mass_change':phase_change/mass(initial),
        'contour_reference_volume_change_relative':graph_volume(a[:,0]/500,(b[:,1]-a[:,1])/500,0)/base,
        'radial_resolution_results':reports,
        'limitations':['Saved static one-way trajectory, not two-way feedback prediction',
            'Only saved snapshots available; update-interval and plastic-path sensitivity are not tested',
            '401 radial graph nodes is a projection comparison, not a refined phase or mechanics mesh']}
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savetxt(args.output/'surface_changes.csv',np.column_stack((query*500,contour*500,*[h*500 for h in profiles])),delimiter=',',header='r_nm,contour_change_nm,mass_51_nm,mass_101_nm,mass_201_nm,mass_401_nm',comments='')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(7,7),sharex=True)
    axes[0].plot(query*500,contour*500,label='phi=0.5 contour',color='black')
    for row,h in zip(reports,profiles):
        axes[0].plot(query*500,h*500,label=f"Mass graph, {row['radial_nodes']} nodes",alpha=.8)
        axes[1].plot(query*500,(h-contour)*500,label=f"{row['radial_nodes']} nodes")
    axes[0].set_ylabel('Height change (nm)');axes[0].legend(fontsize=8)
    axes[1].set_ylabel('Mass graph minus contour (nm)');axes[1].set_xlabel('Reference radius (nm)')
    fig.suptitle('Saved 100 s one-way trajectory: geometry-transfer audit')
    fig.tight_layout();fig.savefig(args.output/'surface_comparison.png',dpi=180);plt.close(fig)
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
