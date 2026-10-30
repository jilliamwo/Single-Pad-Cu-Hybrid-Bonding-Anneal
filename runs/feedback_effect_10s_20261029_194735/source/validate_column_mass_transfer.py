"""Analytic and independent FE checks for the mass-equivalent geometry map."""
import json
from pathlib import Path
import numpy as np
import ufl
from dolfinx import mesh as dm,fem
from mpi4py import MPI
from column_mass_transfer import ColumnMassTransfer
from conservative_occupancy import graph_volume,conservative_occupancy
from thermoelastic_solver import ReferenceGeometry
from validate_baseline import create_baseline_mesh


def main():
    mesh=dm.create_rectangle(MPI.COMM_WORLD,[[0.,0.],[.5,1.1]],[8,12],cell_type=dm.CellType.triangle)
    V=fem.functionspace(mesh,('Lagrange',1));initial=fem.Function(V);current=fem.Function(V)
    initial.interpolate(lambda x:.5*(1-np.tanh((x[1]-.993)/.02)))
    radii=np.array([0.,.037,.12,.26,.39,.5]) # Deliberately cross phase cell boundaries.
    transfer=ColumnMassTransfer(initial,radii)
    current.x.array[:]=initial.x.array
    height,_=transfer.project(current);np.testing.assert_array_equal(height,0)
    # Linear radial column change is reproduced exactly in the radial P1 space.
    x=V.tabulate_dof_coordinates();delta=1e-4*(1-3*x[:,0])
    current.x.array[:]=initial.x.array+delta
    height,diag=transfer.project(current)
    np.testing.assert_allclose(height,1.1e-4*(1-3*radii),atol=1e-15,rtol=1e-10)
    assert abs(diag['phase_mass_change_without_2pi'])<1e-16
    # A varying field tests total mass with independent FEniCS integration.
    current.x.array[:]=initial.x.array+1e-4*np.sin(7*x[:,0])*np.cos(3*x[:,1])
    height,diag=transfer.project(current)
    r=ufl.SpatialCoordinate(mesh)[0]
    expected=float(fem.assemble_scalar(fem.form((current-initial)*r*ufl.dx)))
    np.testing.assert_allclose(diag['phase_mass_change_without_2pi'],expected,rtol=1e-10,atol=1e-16)
    np.testing.assert_allclose(graph_volume(radii,height,0),expected,rtol=1e-10,atol=1e-16)
    g=ReferenceGeometry();mechanics=create_baseline_mesh(g,8,8)
    fraction,details=conservative_occupancy(mechanics,g,radii,g.cu_top+height)
    base=graph_volume(radii,np.full(len(radii),g.cu_top),1.)
    np.testing.assert_allclose(details['mapped_volume_without_2pi']-base,expected,atol=1e-15,rtol=1e-9)
    # Projection is linear: out-of-range phi is not silently clipped to [0,1].
    current.x.array[:]=initial.x.array+1e-4
    height,_=transfer.project(current)
    np.testing.assert_allclose(height,1.1e-4,atol=1e-15)
    snapshot=transfer.initial.copy();current.x.array[:]=0
    np.testing.assert_array_equal(transfer.initial,snapshot)
    current.x.array[0]=np.nan
    try:transfer.project(current)
    except ValueError:pass
    else:raise AssertionError('Nonfinite phase must be rejected')
    report={'passed':True,'checks':['zero change preserves initial geometry','analytic radial redistribution reproduced',
        'zero total mass redistribution conserves graph volume','nonaligned radial intervals integrate P1 exactly',
        'independent FE mass integral agrees','fractional mechanics volume matches phase mass change',
        'initial phase snapshot immutable','nonfinite input rejected'],'example':diag,
        'scope':'mass-equivalent column projection; not the phi=0.5 contour or physical validation'}
    Path('new/verification_results/column_mass_transfer.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
