"""Experimental transactional mechanics on a fixed Cu/void cell map.

Retained Cu keeps local integration-point history; removed material loses it.
Growth requires an explicit new-material policy. This is not plastic advection,
not a conservative diffusion-to-material transfer, and not enabled by run_model.
"""
from dataclasses import dataclass
import numpy as np
import ufl
from dolfinx import fem
from plasticity import PlasticState
from thermoelastic_solver import solve
from phase_geometry_map import SurfaceReferenceMap
from conservative_occupancy import conservative_occupancy


def occupancy_from_surface(mesh,g,radii,reference_height):
    radii=np.asarray(radii,dtype=float);height=np.asarray(reference_height,dtype=float)
    if radii.ndim!=1 or len(radii)<2 or radii.shape!=height.shape or np.any(np.diff(radii)<=0):
        raise ValueError('Surface must have strictly increasing radii and matching heights')
    if not np.all(np.isfinite(radii)) or not np.all(np.isfinite(height)):
        raise ValueError('Surface must be finite')
    base=g.base_thickness_nm/g.L0_nm;radius=g.pad_radius_nm/g.L0_nm
    if radii[0]>1e-12 or radii[-1]<radius-1e-12 or np.any(height<=base) or np.any(height>g.top+1e-12):
        raise ValueError('Surface exceeds the supported mechanics domain or radial coverage')
    occupancy=fem.Function(fem.functionspace(mesh,('DG',0)),name='Cu_occupancy')
    # Binary cell-centroid sampling is intentionally explicit and resolution limited.
    occupancy.interpolate(lambda x:((x[0]<=radius)&(x[1]>=base)&(x[1]<=np.interp(x[0],radii,height))).astype(float))
    occupancy.x.scatter_forward()
    return occupancy


def prepare_history(state,old_occupancy,new_occupancy,policy):
    if policy not in ('reject_growth','virgin_at_reference'):
        raise ValueError('Unknown new-Cu history policy')
    Q=state.equivalent_strain.function_space
    if old_occupancy.function_space.mesh is not Q.mesh or new_occupancy.function_space.mesh is not Q.mesh:
        raise ValueError('Occupancy and history must share a mesh')
    sampled=[]
    for source in (old_occupancy,new_occupancy):
        f=fem.Function(Q)
        f.interpolate(fem.Expression(source,Q.element.interpolation_points))
        if not np.all(np.isfinite(f.x.array)) or np.any(f.x.array<0) or np.any(f.x.array>1):
            raise ValueError('History occupancy must be finite in [0,1]')
        sampled.append(f.x.array.copy())
    old,new=sampled;gained=new>old+1e-12;lost=new<old-1e-12
    if np.any(gained) and policy=='reject_growth':
        raise ValueError('New Cu requires an explicit history policy; no material history is inferred')
    ep=fem.Function(state.strain.function_space,name='plastic_strain_trial')
    p=fem.Function(Q,name='equivalent_plastic_strain_trial')
    ep.x.array[:]=state.strain.x.array;p.x.array[:]=state.equivalent_strain.x.array
    # Homogeneous dilution by newly virgin Cu; loss retains intrinsic history.
    scale=np.divide(np.minimum(old,new),new,out=np.zeros_like(new),where=new>0)
    ep.x.array.reshape(-1,3,3)[:]*=scale[:,None,None];p.x.array[:]*=scale
    ep.x.scatter_forward();p.x.scatter_forward()
    return PlasticState(ep,p),{'gained_integration_points':int(gained.sum()),
        'removed_integration_points':int(lost.sum()),'retained_integration_points':int(((old>0)&(new>0)).sum()),
        'new_cu_history_policy':policy}


@dataclass
class MechanicalTrial:
    owner: object
    version: int
    occupancy: object
    result: object
    diagnostics: dict
    reference_fraction: object = None


class EvolvingMechanics:
    """Propose independently; commit only after all coupled substeps succeed.

    A failed/discarded proposal leaves accepted state unchanged. The outer phase
    integrator must also roll back its phase and energy fields before retrying.
    """
    def __init__(self,mesh,geometry,cu,dielectric,temperature_K,law,initial_result,*,new_cu_history='reject_growth',occupancy_mode='binary'):
        if mesh.comm.size!=1:raise ValueError('Experimental evolving mechanics requires one rank')
        if initial_result.plastic_state is None or initial_result.u_dimless_field.function_space.mesh is not mesh:
            raise ValueError('Provide a plastic mechanics result on the same mesh')
        if new_cu_history not in ('reject_growth','virgin_at_reference'):raise ValueError('Unknown new-Cu history policy')
        if initial_result.stats.get('explicit_cu_occupancy',False):
            raise ValueError('Initialize from the original sharp-geometry result, not an evolved occupancy result')
        if occupancy_mode not in ('binary','fractional'):raise ValueError('Unknown occupancy mode')
        self.occupancy_mode=occupancy_mode
        self.mesh=mesh;self.g=geometry;self.cu=cu;self.dielectric=dielectric
        self.temperature=temperature_K;self.law=law;self.policy=new_cu_history
        self.result=initial_result;self.version=0
        radius=geometry.pad_radius_nm/geometry.L0_nm
        self.occupancy=occupancy_from_surface(mesh,geometry,[0.,radius],[geometry.cu_top]*2)
        self.mapping=SurfaceReferenceMap(initial_result.UY_r_dimless,initial_result.UY_dimless,
            base_z=geometry.base_thickness_nm/geometry.L0_nm,reference_cu_top=geometry.cu_top,mechanical_top=geometry.top)

    def propose(self,radii,phase_surface):
        reference=self.mapping.reference_surface(radii,phase_surface)
        occupancy=occupancy_from_surface(self.mesh,self.g,radii,reference)
        fraction,geometry_diagnostics=conservative_occupancy(self.mesh,self.g,radii,reference)
        if self.occupancy_mode=='fractional':occupancy=fraction
        state,diagnostics=prepare_history(self.result.plastic_state,self.occupancy,occupancy,self.policy)
        diagnostics['occupancy_mode']=self.occupancy_mode
        diagnostics['conservative_geometry']=geometry_diagnostics
        diagnostics['changed_cells']=int(np.count_nonzero(self.occupancy.x.array!=occupancy.x.array))
        r=ufl.SpatialCoordinate(self.mesh)[0]
        previous_volume=float(fem.assemble_scalar(fem.form(self.occupancy*r*ufl.dx)))
        trial_volume=float(fem.assemble_scalar(fem.form(occupancy*r*ufl.dx)))
        diagnostics['previous_reference_Cu_volume_without_2pi']=previous_volume
        diagnostics['trial_reference_Cu_volume_without_2pi']=trial_volume
        diagnostics['relative_reference_Cu_volume_change']=(trial_volume-previous_volume)/previous_volume
        graph_volume=geometry_diagnostics['graph_volume_without_2pi']
        diagnostics[self.occupancy_mode+'_volume_error_relative_to_graph']=(trial_volume-graph_volume)/graph_volume

        result=solve(None,self.mesh,self.cu,self.temperature,self.g.L0_m,
            self.g.pad_radius_nm/self.g.L0_nm,(self.g.pad_height_nm-self.g.recess_nm)/self.g.L0_nm,
            geometry=self.g,material_dielectric=self.dielectric,plasticity_cu=self.law,
            load_factors=[1.],initial_plastic_state=state,
            initial_displacement=self.result.u_dimless_field,cu_occupancy=occupancy,
            fractional_occupancy=self.occupancy_mode=='fractional')
        return MechanicalTrial(self,self.version,occupancy,result,diagnostics,reference_fraction=fraction)

    def commit(self,trial):
        if trial.owner is not self or trial.version!=self.version:
            raise ValueError('Cannot commit a foreign or stale mechanical trial')
        self.result=trial.result;self.occupancy=trial.occupancy;self.version+=1
