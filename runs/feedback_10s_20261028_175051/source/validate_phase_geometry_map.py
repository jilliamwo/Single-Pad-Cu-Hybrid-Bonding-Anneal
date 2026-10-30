"""Small-mesh verification of phase/reference coordinate transfer."""
from pathlib import Path
import json
import numpy as np
from mpi4py import MPI
from dolfinx import mesh as dm,fem
from phase_geometry_map import SurfaceReferenceMap,pull_back_phase_pair


def main():
    r=np.linspace(0,.5,7);uz=.004+.002*r
    mapping=SurfaceReferenceMap(r,uz)
    initial_height=.993+uz
    np.testing.assert_allclose(mapping.reference_surface(r,initial_height),1.993,atol=1e-14)
    change=.001*np.cos(np.pi*r/.5)
    np.testing.assert_allclose(mapping.reference_surface(r,initial_height+change),1.993+change,atol=1e-14)
    # Retain the initial map even if caller-owned arrays are later modified.
    uz[:]=9.
    np.testing.assert_allclose(mapping.shift(r),.004+.002*r)
    for radii,height in [(np.array([.6]),np.array([1.])),(r,np.full(len(r),1.1))]:
        try:mapping.reference_surface(radii,height)
        except ValueError:pass
        else:raise AssertionError('Out-of-domain surface must not be silently clipped')
    phase=dm.create_rectangle(MPI.COMM_WORLD,[np.array([0.,0.]),np.array([.5,1.1])],[4,8])
    mechanics=dm.create_rectangle(MPI.COMM_WORLD,[np.array([0.,0.]),np.array([1.,2.])],[8,8])
    V=fem.functionspace(phase,('Lagrange',1));Q=fem.functionspace(mechanics,('Lagrange',1))
    initial=fem.Function(V);current=fem.Function(V)
    initial.interpolate(lambda x:.7-.2*x[1]+.1*x[0])
    current.interpolate(lambda x:.73-.2*x[1]+.1*x[0])
    now,start,delta=pull_back_phase_pair(mapping,current,initial,Q)
    xyz=Q.tabulate_dof_coordinates();inside=(xyz[:,0]>=-1e-12)&(xyz[:,0]<=.5+1e-12)&(xyz[:,1]>=1.-1e-12)
    query=mapping.phase_points(xyz[inside])
    np.testing.assert_allclose(start.x.array[inside],.7-.2*query[:,1]+.1*query[:,0],atol=1e-14)
    np.testing.assert_allclose(delta.x.array[inside],.03,atol=1e-14)
    np.testing.assert_array_equal(now.x.array[~inside],0.)
    _,_,unchanged=pull_back_phase_pair(mapping,initial,initial,Q)
    np.testing.assert_array_equal(unchanged.x.array,0.)
    bad=SurfaceReferenceMap(r,np.full(len(r),.3))
    try:pull_back_phase_pair(bad,current,initial,Q)
    except ValueError as error:assert 'leaves phase mesh' in str(error)
    else:raise AssertionError('Unmapped points must be reported')
    report={'passed':True,'checks':['initial thermal shift removed once','phase surface change preserved',
        'map immutable from caller arrays','out-of-domain surfaces rejected','affine FE transfer exact',
        'dielectric/below-pad exclusion','identical fields give exact zero change','unmapped phase queries rejected'],
        'scope':'frozen-column coordinate map only; no plastic transport or two-way evolution'}
    out=Path('new/verification_results/phase_geometry_map.json');out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
