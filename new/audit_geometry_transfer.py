"""Compare binary and conservative graph transfer using an existing run; no solve."""
import argparse
import json
from pathlib import Path
import numpy as np
import ufl
from dolfinx import fem
from thermoelastic_solver import ReferenceGeometry
from validate_baseline import create_baseline_mesh
from evolving_mechanics import occupancy_from_surface
from conservative_occupancy import conservative_occupancy


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,default=Path('runs/native_dt05_100s'))
    parser.add_argument('--output',type=Path,default=Path('runs/geometry_transfer_audit'))
    args=parser.parse_args()
    cfg=json.loads((args.run/'config.json').read_text())
    g=ReferenceGeometry();mesh=create_baseline_mesh(g,cfg['mechanics_nr'],cfg['mechanics_nz'])
    a=np.loadtxt(args.run/'initial_surface.csv',delimiter=',',skiprows=1)
    b=np.loadtxt(args.run/'final_surface.csv',delimiter=',',skiprows=1)
    np.testing.assert_allclose(a[:,0],b[:,0],rtol=0,atol=1e-10)
    r=a[:,0]/g.L0_nm;flat=np.full(len(r),g.cu_top)
    height=flat+(b[:,1]-a[:,1])/g.L0_nm
    initial,first=conservative_occupancy(mesh,g,r,flat)
    final,last=conservative_occupancy(mesh,g,r,height)
    binary0=occupancy_from_surface(mesh,g,r,flat);binary1=occupancy_from_surface(mesh,g,r,height)
    x=ufl.SpatialCoordinate(mesh)
    volume=lambda f:float(fem.assemble_scalar(fem.form(f*x[0]*ufl.dx)))
    v0,v1=volume(binary0),volume(binary1)
    history=np.genfromtxt(args.run/'time_history.csv',delimiter=',',names=True)
    report={'source_run':str(args.run),'no_mechanical_or_phase_solve':True,
        'initial_geometry':first,'final_geometry':last,
        'binary_changed_cells':int(np.count_nonzero(binary0.x.array!=binary1.x.array)),
        'fraction_changed_cells':int(np.count_nonzero(abs(final.x.array-initial.x.array)>1e-12)),
        'binary_initial_volume_without_2pi':v0,'binary_final_volume_without_2pi':v1,
        'binary_final_error_relative_to_graph':(v1-last['graph_volume_without_2pi'])/last['graph_volume_without_2pi'],
        'reference_graph_volume_change':last['graph_volume_without_2pi']/first['graph_volume_without_2pi']-1,
        'reported_phase_mass_drift':float(np.atleast_1d(history['relative_mass_drift'])[-1]),
        'interpretation':'Fraction transfer conserves reconstructed reference graph. Graph extraction/frozen-column mapping is not a conservative transfer of diffuse phase mass.'}
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'geometry_transfer.json').write_text(json.dumps(report,indent=2)+'\n')
    coordinates=final.function_space.tabulate_dof_coordinates()
    np.savetxt(args.output/'cell_fractions.csv',np.column_stack((coordinates[:,:2],initial.x.array,final.x.array,binary1.x.array)),delimiter=',',header='r_reference_dimless,z_reference_dimless,initial_Cu_fraction,final_Cu_fraction,binary_final_Cu',comments='')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
