"""Geometry-only tests of axisymmetric conservative cell fractions; no FEA solve."""
from pathlib import Path
import json
import numpy as np
import ufl
from dolfinx import fem
from conservative_occupancy import triangle_fractions,conservative_occupancy
from evolving_mechanics import occupancy_from_surface
from thermoelastic_solver import ReferenceGeometry
from validate_baseline import create_baseline_mesh


def main():
    t=np.array([[[0.,0.],[1.,0.],[0.,1.]]])
    f,w=triangle_fractions(t,[0.,1.],[.5,.5],radius=1.,base=0.,top=1.)
    np.testing.assert_allclose(w,[1/6],rtol=0,atol=1e-15)
    np.testing.assert_allclose(f,[7/8],rtol=0,atol=1e-15)
    reverse,_=triangle_fractions(t[:,::-1],[0.,1.],[.5,.5],radius=1.,base=0.,top=1.)
    np.testing.assert_allclose(reverse,f,rtol=0,atol=1e-15)
    g=ReferenceGeometry();r=np.linspace(0,.5,17);flat=np.full(len(r),g.cu_top)
    reports=[]
    for n in (4,8):
        mesh=create_baseline_mesh(g,n,n)
        original,base_report=conservative_occupancy(mesh,g,r,flat)
        binary=occupancy_from_surface(mesh,g,r,flat)
        np.testing.assert_allclose(original.x.array,binary.x.array,atol=2e-12,rtol=0)
        # Movement is far smaller than any cell height or centroid offset.
        moved,moved_report=conservative_occupancy(mesh,g,r,flat+1e-5)
        moved_binary=occupancy_from_surface(mesh,g,r,flat+1e-5)
        np.testing.assert_array_equal(moved_binary.x.array,binary.x.array)
        assert np.max(abs(moved.x.array-original.x.array))>1e-8
        np.testing.assert_allclose(moved_report['mapped_volume_without_2pi']-base_report['mapped_volume_without_2pi'],.5*.5**2*1e-5,atol=1e-15,rtol=1e-9)
        # Linear tilt has exactly zero axisymmetric volume change.
        shaped=flat+1e-4*(1-1.5*r/.5)
        tilted,tilted_report=conservative_occupancy(mesh,g,r,shaped)
        np.testing.assert_allclose(tilted_report['mapped_volume_without_2pi'],base_report['mapped_volume_without_2pi'],atol=1e-14,rtol=0)
        assert tilted_report['fractional_cells']>0
        # Compare actual FEniCS radial integration with independent graph integral.
        x=ufl.SpatialCoordinate(mesh)
        actual=float(fem.assemble_scalar(fem.form(tilted*x[0]*ufl.dx)))
        np.testing.assert_allclose(actual,tilted_report['graph_volume_without_2pi'],atol=1e-14,rtol=0)
        # Nonlinear sampled graph crosses multiple linear segments per coarse cell.
        curved=flat-.01*(r/.5)**2
        curved_field,curved_report=conservative_occupancy(mesh,g,r,curved)
        assert abs(curved_report['relative_graph_volume_error'])<1e-12
        reports.append({'mesh_nr':n,'flat':base_report,'subcell_motion':moved_report,'constant_volume_tilt':tilted_report,'piecewise_linear_curve':curved_report})
        for bad in (np.full(len(r),2.01),np.full(len(r),np.nan)):
            try:conservative_occupancy(mesh,g,r,bad)
            except ValueError:pass
            else:raise AssertionError('Invalid graph must not be clipped silently')
    report={'passed':True,'checks':['exact radial moment (not planar area fraction)','triangle orientation invariance',
        'sharp baseline reproduction','subcell motion detected where binary mask is unchanged',
        'known added volume recovered','zero-volume redistribution preserved','FE radial integral agrees',
        'multiple surface segments per cell','invalid input rejected'],
        'results':reports,'scope':'conservation of supplied reference graph only; not phase mass or fractional plasticity validation'}
    path=Path('new/verification_results/conservative_occupancy.json');path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
