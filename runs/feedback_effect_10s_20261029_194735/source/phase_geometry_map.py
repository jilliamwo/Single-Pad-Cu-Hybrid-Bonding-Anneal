"""Reference-coordinate transfer for the model's graph-initialized phase field.

The current initializer uses reference r and a vertical surface displacement,
not a full deformed-coordinate mesh. This map deliberately follows that convention.
It is a frozen column shift, NOT material advection or a conservative remap.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class SurfaceReferenceMap:
    radii: np.ndarray
    initial_uz: np.ndarray
    base_z: float = 1.
    reference_cu_top: float = 1.993
    mechanical_top: float = 2.

    def __post_init__(self):
        r=np.array(self.radii,dtype=float,copy=True)
        u=np.array(self.initial_uz,dtype=float,copy=True)
        if r.ndim!=1 or len(r)<2 or r.shape!=u.shape:
            raise ValueError('Provide matching one-dimensional radial/displacement profiles')
        if not np.all(np.isfinite(r)) or not np.all(np.isfinite(u)) or np.any(np.diff(r)<=0):
            raise ValueError('Profile must be finite with strictly increasing radii')
        if not np.all(np.isfinite([self.base_z,self.reference_cu_top,self.mechanical_top])) or not self.base_z<self.reference_cu_top<=self.mechanical_top:
            raise ValueError('Invalid reference geometry')
        r.setflags(write=False);u.setflags(write=False)
        object.__setattr__(self,'radii',r);object.__setattr__(self,'initial_uz',u)

    def shift(self,r):
        r=np.asarray(r,dtype=float)
        if not np.all(np.isfinite(r)) or np.any(r<self.radii[0]-1e-12) or np.any(r>self.radii[-1]+1e-12):
            raise ValueError('Radius outside the calibrated initial surface profile')
        return np.interp(r,self.radii,self.initial_uz)

    def phase_points(self,reference_points):
        """Return query points in phase coordinates, including initial vertical shift."""
        points=np.array(reference_points,dtype=float,copy=True)
        if points.ndim!=2 or points.shape[1] not in (2,3) or not np.all(np.isfinite(points)):
            raise ValueError('Expected finite N-by-2 or N-by-3 points')
        shift=self.shift(points[:,0])  # Validate before snapping boundary roundoff.
        points[:,0]=np.clip(points[:,0],self.radii[0],self.radii[-1])
        points[:,1]=points[:,1]-self.base_z+shift
        return points

    def reference_surface(self,r,phase_height,*,check_bounds=True):
        """Remove initial displacement once; retain subsequent phase-shape changes."""
        r=np.asarray(r,dtype=float);height=np.asarray(phase_height,dtype=float)
        if height.shape!=r.shape or not np.all(np.isfinite(height)):
            raise ValueError('Surface height and radii must have matching finite shapes')
        reference=height+self.base_z-self.shift(r)
        if check_bounds and (np.any(reference<=self.base_z) or np.any(reference>self.mechanical_top+1e-12)):
            raise ValueError('Evolved surface exceeds the supported mechanical domain; do not clip it')
        return reference


def pull_back_phase_pair(mapping,current,initial,mechanics_scalar_space):
    """Sample both phase states consistently on a mechanics P1 space.

    Returns current, initial, and difference Functions. Outside the pad column
    and below its bottom they are zero: dielectric identity is handled separately.
    No plastic state is changed. Both inputs must use the same phase mesh/space.
    Failure to locate a query point is fatal, never a silent void assignment.
    """
    from dolfinx import fem,geometry
    if current.function_space.mesh is not initial.function_space.mesh or current.ufl_shape or initial.ufl_shape:
        raise ValueError('Both phase states must be scalar and on the same mesh')
    if current.function_space.ufl_element()!=initial.function_space.ufl_element():
        raise ValueError('Phase states must use the same finite element')
    mesh=current.function_space.mesh
    if mesh.comm.size!=1 or mechanics_scalar_space.mesh.comm.size!=1:
        raise ValueError('Geometry transfer currently requires one MPI rank')
    if np.asarray(mechanics_scalar_space.element.value_shape).size:
        raise ValueError('Mechanics target must be scalar')
    coords=mechanics_scalar_space.tabulate_dof_coordinates()
    inside=(coords[:,0]>=mapping.radii[0]-1e-12)&(coords[:,0]<=mapping.radii[-1]+1e-12)&(coords[:,1]>=mapping.base_z-1e-12)
    points=mapping.phase_points(coords[inside])
    tree=geometry.bb_tree(mesh,mesh.topology.dim)
    candidates=geometry.compute_colliding_cells(mesh,geometry.compute_collisions_points(tree,points),points)
    cells=np.array([candidates.links(i)[0] if len(candidates.links(i)) else -1 for i in range(len(points))],dtype=np.int32)
    if np.any(cells<0):
        raise ValueError('Reference-to-phase query leaves phase mesh; extend the domain or define a boundary rule')
    outputs=[]
    for source,name in ((current,'phi_reference_current'),(initial,'phi_reference_initial')):
        target=fem.Function(mechanics_scalar_space,name=name)
        values=source.eval(points,cells).ravel()
        if not np.all(np.isfinite(values)):raise ValueError('Nonfinite transferred phase values')
        target.x.array[inside]=values;target.x.scatter_forward();outputs.append(target)
    delta=fem.Function(mechanics_scalar_space,name='phi_reference_change')
    delta.x.array[:]=outputs[0].x.array-outputs[1].x.array;delta.x.scatter_forward()
    return *outputs,delta
