"""Axisymmetric Cu volume fractions for a piecewise-linear reference surface.

Each DG0 coefficient is integral_Cu(r dA) / integral_cell(r dA), computed by
polygon clipping. This preserves the supplied graph volume, not automatically
phase-field mass after coordinate transfer. Fractional plasticity is not implied.
"""
import numpy as np


def clip_polygon(vertices,a,b,c):
    """Clip a convex polygon to a*r + b*z <= c."""
    if len(vertices)==0:return np.empty((0,2))
    output=[]
    previous=vertices[-1];fp=a*previous[0]+b*previous[1]-c
    for point in vertices:
        f=a*point[0]+b*point[1]-c
        if (f<=0)!=(fp<=0):
            output.append(previous+(point-previous)*(fp/(fp-f)))
        if f<=0:output.append(point)
        previous=point;fp=f
    return np.asarray(output,dtype=float).reshape(-1,2)


def radial_moment(vertices):
    """Integral r dA by fan triangulation; 2*pi is omitted throughout."""
    if len(vertices)<3:return 0.
    origin=vertices[0];total=0.
    for i in range(1,len(vertices)-1):
        p,q=vertices[i],vertices[i+1]
        v,w=p-origin,q-origin
        area=.5*abs(v[0]*w[1]-v[1]*w[0])
        total+=area*(origin[0]+p[0]+q[0])/3
    return float(total)


def validate_surface(radii,height,radius,base,top):
    r=np.array(radii,dtype=float,copy=True);h=np.asarray(height,dtype=float)
    if r.ndim!=1 or len(r)<2 or h.shape!=r.shape:
        raise ValueError('Matching one-dimensional surface arrays required')
    if not np.all(np.isfinite(r)) or not np.all(np.isfinite(h)) or np.any(np.diff(r)<=0):
        raise ValueError('Surface must be finite with increasing radii')
    if abs(r[0])>1e-12 or abs(r[-1]-radius)>1e-12 or np.any(h<=base) or np.any(h>top+1e-12):
        raise ValueError('Surface must span the pad exactly and lie in the mechanics domain')
    r[0]=0.;r[-1]=radius
    return r,np.minimum(h,top)


def graph_volume(radii,height,base):
    """Exact integral r*(piecewise-linear height-base) dr using two-point Gauss."""
    left,right=radii[:-1],radii[1:];mid=(left+right)/2;half=(right-left)/2
    total=0.
    for sign in (-1.,1.):
        q=mid+sign*half/np.sqrt(3)
        h=height[:-1]+(height[1:]-height[:-1])*(q-left)/(right-left)
        total+=np.sum(half*q*(h-base))
    return float(total)


def triangle_fractions(triangles,radii,height,*,radius,base,top):
    """Return axisymmetric fractions and full-cell radial volume weights."""
    r,h=validate_surface(radii,height,radius,base,top)
    triangles=np.asarray(triangles,dtype=float)
    if triangles.ndim!=3 or triangles.shape[1:]!=(3,2) or not np.all(np.isfinite(triangles)):
        raise ValueError('Expected finite N-by-3-by-2 affine triangles')
    if np.any(triangles[:,:,0]<-1e-12):raise ValueError('Negative radius in mechanics mesh')
    weights=np.array([radial_moment(t) for t in triangles])
    if np.any(weights<=0):raise ValueError('Degenerate or zero radial-volume cell')
    fractions=np.zeros(len(triangles));minimum=float(h.min());maximum=float(h.max())
    for index,t in enumerate(triangles):
        lo=t.min(axis=0);hi=t.max(axis=0)
        if lo[0]>=radius or hi[1]<=base or lo[1]>=maximum:continue
        if hi[0]<=radius and lo[1]>=base and hi[1]<=minimum:
            fractions[index]=1.;continue
        polygon=clip_polygon(t,-1.,0.,0.)
        polygon=clip_polygon(polygon,1.,0.,radius)
        polygon=clip_polygon(polygon,0.,-1.,-base)
        if len(polygon)<3:continue
        begin=max(0,np.searchsorted(r,lo[0],side='right')-1)
        end=min(len(r)-1,np.searchsorted(r,hi[0],side='left')+1)
        occupied=0.
        for k in range(begin,end):
            if r[k]>=hi[0] or r[k+1]<=lo[0]:continue
            part=clip_polygon(polygon,-1.,0.,-r[k])
            part=clip_polygon(part,1.,0.,r[k+1])
            slope=(h[k+1]-h[k])/(r[k+1]-r[k])
            part=clip_polygon(part,-slope,1.,h[k]-slope*r[k])
            occupied+=radial_moment(part)
        fractions[index]=occupied/weights[index]
    if np.any(fractions < -1e-10) or np.any(fractions > 1+1e-10):
        raise RuntimeError('Clipped cell fraction outside [0,1]')
    return np.clip(fractions,0.,1.),weights


def conservative_occupancy(mesh,g,radii,reference_height):
    """Build a DG0 volume map on a serial affine triangular mechanics mesh."""
    from dolfinx import fem
    if mesh.comm.size!=1:raise ValueError('Conservative occupancy currently requires one rank')
    if mesh.topology.dim!=2 or mesh.geometry.dofmap.shape[1]!=3:
        raise ValueError('Conservative occupancy requires affine triangles')
    radius=g.pad_radius_nm/g.L0_nm;base=g.base_thickness_nm/g.L0_nm
    r,h=validate_surface(radii,reference_height,radius,base,g.top)
    triangles=mesh.geometry.x[mesh.geometry.dofmap,:2]
    fraction,weights=triangle_fractions(triangles,r,h,radius=radius,base=base,top=g.top)
    V=fem.functionspace(mesh,('DG',0));field=fem.Function(V,name='Cu_axisymmetric_volume_fraction')
    for cell,value in enumerate(fraction):field.x.array[V.dofmap.cell_dofs(cell)[0]]=value
    field.x.scatter_forward()
    target=graph_volume(r,h,base);actual=float(np.dot(fraction,weights))
    error=(actual-target)/target
    if abs(error)>1e-9:
        raise ValueError('Mechanics mesh does not cover the supplied Cu graph volume')
    diagnostics={'graph_volume_without_2pi':target,'mapped_volume_without_2pi':actual,
                 'relative_graph_volume_error':error,
                 'fractional_cells':int(np.sum((fraction>1e-12)&(fraction<1-1e-12))),
                 'min_fraction':float(fraction.min()),'max_fraction':float(fraction.max())}
    return field,diagnostics
