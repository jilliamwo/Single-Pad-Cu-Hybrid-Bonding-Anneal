"""Conservative projection of P1 phase-mass changes to a radial height graph.

The radial weighted L2 projection satisfies integral(r*dh dr) =
integral(r*(phi-phi_initial) dA). It preserves changes, anchoring the initial
sharp mechanics geometry; it does not interpret the diffuse field as a contour.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.linalg import cho_factor, cho_solve
from conservative_occupancy import clip_polygon


class ColumnMassTransfer:
    def __init__(self,initial,radii):
        V=initial.function_space;mesh=V.mesh
        from basix.ufl import element
        if mesh.comm.size!=1 or mesh.topology.dim!=2 or mesh.geometry.dofmap.shape[1]!=3:
            raise ValueError('Column transfer requires serial affine triangles')
        if V.ufl_element()!=element('Lagrange',mesh.basix_cell(),1):
            raise ValueError('Column transfer requires scalar P1 phase')
        self.V=V
        coords=V.tabulate_dof_coordinates()[:,:2]
        cells=np.array([V.dofmap.cell_dofs(c) for c in range(mesh.topology.index_map(2).size_local)])
        self.projector=NodalColumnMassTransfer(coords,cells,initial.x.array,radii)
        self.initial=self.projector.initial

    def project(self,current):
        if current.function_space is not self.V:
            raise ValueError('Phase transfer requires the original function space')
        return self.projector.project(current.x.array)


class NodalColumnMassTransfer:
    """Same projection for saved affine-triangle nodal arrays; no FEniCS required."""
    def __init__(self,coords,cells,initial,radii):
        coords=np.asarray(coords,dtype=float);cells=np.asarray(cells)
        initial=np.asarray(initial,dtype=float)
        if coords.ndim!=2 or coords.shape[1]!=2 or not np.all(np.isfinite(coords)):
            raise ValueError('Expected finite N-by-2 coordinates')
        if cells.ndim!=2 or cells.shape[1]!=3 or not np.issubdtype(cells.dtype,np.integer) or np.any(cells<0) or np.any(cells>=len(coords)):
            raise ValueError('Expected valid triangle connectivity')
        if initial.shape!=(len(coords),) or not np.all(np.isfinite(initial)):
            raise ValueError('Expected finite initial nodal values')
        r=np.array(radii,dtype=float,copy=True)
        if r.ndim!=1 or len(r)<2 or not np.all(np.isfinite(r)) or np.any(np.diff(r)<=0) or r[0]!=0:
            raise ValueError('Provide increasing finite radii starting at zero')
        if not np.isclose(coords[:,0].min(),0,atol=1e-12) or not np.isclose(coords[:,0].max(),r[-1],rtol=0,atol=1e-12):
            raise ValueError('Phase mesh must span the complete radial transfer interval')
        self.initial=np.array(initial,dtype=float,copy=True);self.radii=r
        # Exact degree-three quadrature integrates phi(P1)*radial_hat(P1)*r.
        bary=np.array([[1/3]*3,[.6,.2,.2],[.2,.6,.2],[.2,.2,.6]])
        weights=np.array([-27/48,25/48,25/48,25/48])
        rows=[];cols=[];values=[]
        for dofs in cells:
            triangle=coords[dofs]
            inverse=np.linalg.inv(np.column_stack((np.ones(3),triangle)))
            lo,hi=triangle[:,0].min(),triangle[:,0].max()
            begin=max(0,np.searchsorted(r,lo,side='right')-1)
            end=min(len(r)-1,np.searchsorted(r,hi,side='left'))
            for k in range(begin,end):
                polygon=clip_polygon(clip_polygon(triangle,-1,0,-r[k]),1,0,r[k+1])
                block=np.zeros((2,3))
                for j in range(1,len(polygon)-1):
                    t=polygon[[0,j,j+1]];v=t[1]-t[0];w=t[2]-t[0]
                    area=.5*abs(v[0]*w[1]-v[1]*w[0])
                    q=bary@t;shape=np.column_stack((np.ones(4),q))@inverse
                    right=(q[:,0]-r[k])/(r[k+1]-r[k])
                    hats=np.column_stack((1-right,right))
                    block+=hats.T@((area*weights*q[:,0])[:,None]*shape)
                for a in range(2):
                    for b in range(3):rows.append(k+a);cols.append(dofs[b]);values.append(block[a,b])
        self.B=coo_matrix((values,(rows,cols)),shape=(len(r),len(self.initial))).tocsr()
        mass=np.zeros((len(r),len(r)))
        for k in range(len(r)-1):
            half=(r[k+1]-r[k])/2;mid=(r[k+1]+r[k])/2
            for sign in (-1,1):
                q=mid+sign*half/np.sqrt(3);right=(q-r[k])/(2*half)
                hat=np.array([1-right,right]);mass[k:k+2,k:k+2]+=half*q*np.outer(hat,hat)
        self.factor=cho_factor(mass);self.volume_weights=mass.sum(axis=0)

    def project(self,current):
        current=np.asarray(current,dtype=float)
        if current.shape!=self.initial.shape:raise ValueError("Phase array shape changed")
        if not np.all(np.isfinite(current)):raise ValueError('Nonfinite phase state')
        rhs=self.B@(current-self.initial)
        height=cho_solve(self.factor,rhs)
        phase_change=float(rhs.sum());graph_change=float(self.volume_weights@height)
        return height,{'phase_mass_change_without_2pi':phase_change,
                       'graph_volume_change_without_2pi':graph_change,
                       'mass_transfer_absolute_error':graph_change-phase_change,
                       'max_abs_mass_equivalent_height_change':float(np.max(np.abs(height)))}
