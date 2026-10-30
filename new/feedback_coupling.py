"""Experimental staggered feedback adapter; no production CLI default changes."""
from dataclasses import dataclass
import numpy as np
from dolfinx import fem,geometry


@dataclass
class FeedbackTrial:
    mechanical: object
    energy: np.ndarray
    time_s: float
    step: int


class MechanicalFeedback:
    """Phase evolves with old W, then mechanics supplies W for the next step.

    No inner fixed-point iterations or advection of phi by displacement are applied.
    The initial column shift is frozen. Geometry may use contour changes or an opt-in mass-equivalent column graph.
    Constitutive mixing and frozen-coordinate limitations still apply.
    """
    def __init__(self,engine,*,interval_steps=1,max_volume_change=.01,transfer_mode="contour"):
        if not isinstance(interval_steps,int) or interval_steps<1:raise ValueError('Positive update interval required')
        if not np.isfinite(max_volume_change) or not 0<max_volume_change<1:raise ValueError('Provide a fractional volume-change limit between zero and one')
        if transfer_mode not in ('contour','column_mass'):raise ValueError('Unknown geometry transfer mode')
        if transfer_mode=='column_mass' and engine.occupancy_mode!='fractional':
            raise ValueError('Column mass transfer requires fractional mechanics')
        self.transfer_mode=transfer_mode;self.mass_transfer=None
        self.engine=engine;self.interval=interval_steps;self.volume_limit=max_volume_change
        self.records=[];self.initial_surface=None;self.radii=None
        self.initial_volume=None;self.multiplier=None

    def sample_energy(self,result,phase_points):
        mapping=self.engine.mapping
        xyz=np.zeros((len(phase_points),3));xyz[:,:2]=phase_points[:,:2]
        xyz[:,0]=np.clip(xyz[:,0],mapping.radii[0]+1e-10,mapping.radii[-1]-1e-10)
        xyz[:,1]=xyz[:,1]+mapping.base_z-mapping.shift(xyz[:,0])
        # Above the mechanical domain is void under the enforced surface bounds.
        inside=(xyz[:,1]>=0)&(xyz[:,1]<=mapping.mechanical_top)
        output=np.zeros(len(xyz));query=xyz[inside]
        mesh=result.W_field.function_space.mesh
        tree=geometry.bb_tree(mesh,mesh.topology.dim)
        links=geometry.compute_colliding_cells(mesh,geometry.compute_collisions_points(tree,query),query)
        cells=np.array([links.links(i)[0] if len(links.links(i)) else -1 for i in range(len(query))],dtype=np.int32)
        if np.any(cells<0):raise ValueError('Updated energy has unmapped mechanical query points')
        output[inside]=result.W_field.eval(query,cells).ravel()
        if not np.all(np.isfinite(output)) or np.min(output)<-1e-6:raise ValueError('Invalid trial mechanical energy')
        return np.maximum(output,0.)

    def initialize(self,V,radii,sampled_initial_surface,multiplier,*,phase=None):
        if self.initial_surface is not None:raise ValueError('Feedback adapter cannot be reused for a second phase run')
        if self.transfer_mode=='column_mass':
            if phase is None:raise ValueError('Column mass transfer requires initial phase')
            from column_mass_transfer import ColumnMassTransfer
            self.mass_transfer=ColumnMassTransfer(phase,radii)
        self.radii=np.array(radii,copy=True);self.initial_surface=np.array(sampled_initial_surface,copy=True)
        if not np.all(np.isfinite(self.initial_surface)):raise ValueError('Initial surface extraction incomplete')
        self.multiplier=multiplier;self.points=V.tabulate_dof_coordinates()[:,:2]
        W=fem.Function(V,name='mechanical_energy_J_m3')
        W.x.array[:]=multiplier*self.sample_energy(self.engine.result,self.points);W.x.scatter_forward()
        return W

    def propose(self,surface,time_s,step,*,phase=None):
        if step%self.interval:return None
        if not np.all(np.isfinite(surface)):raise ValueError('Feedback requires a complete surface profile')
        # Subtract the sampled initial profile to avoid mistaking contour
        # discretization error for physical movement on the first update.
        analytic=(self.engine.g.cu_top-self.engine.mapping.base_z+self.engine.mapping.shift(self.radii))
        corrected=analytic+surface-self.initial_surface
        transfer_diagnostics={'geometry_transfer':self.transfer_mode}
        if self.mass_transfer is not None:
            if phase is None:raise ValueError('Column mass transfer requires current phase')
            height_change,details=self.mass_transfer.project(phase)
            corrected=analytic+height_change
            transfer_diagnostics.update(details)
            transfer_diagnostics['max_mass_vs_contour_height_difference']=float(np.max(np.abs(height_change-(surface-self.initial_surface))))
        trial=self.engine.propose(self.radii,corrected)
        trial.diagnostics.update(transfer_diagnostics)
        volume=trial.diagnostics['trial_reference_Cu_volume_without_2pi']
        baseline=(trial.diagnostics['previous_reference_Cu_volume_without_2pi']
                  if self.initial_volume is None else self.initial_volume)
        if self.mass_transfer is not None:
            balance=volume-baseline-transfer_diagnostics['phase_mass_change_without_2pi']
            trial.diagnostics['mechanics_phase_mass_balance_error_without_2pi']=balance
            if abs(balance)>1e-10*baseline:
                raise ValueError('Mass-equivalent graph and mechanics occupancy lost phase mass balance')
        drift=(volume-baseline)/baseline
        trial.diagnostics['cumulative_reference_volume_change']=drift
        if abs(drift)>self.volume_limit:
            raise ValueError(f'Mechanics volume drift {drift:.3%} exceeds limit {self.volume_limit:.3%}; transfer needs improvement')
        energy=self.multiplier*self.sample_energy(trial.result,self.points)
        return FeedbackTrial(trial,energy,time_s,step)

    def commit(self,trial,W):
        if trial is None:return
        if trial.energy.shape!=W.x.array.shape or not np.all(np.isfinite(trial.energy)):
            raise ValueError('Invalid energy array at commit')
        if self.initial_volume is None:
            baseline=trial.mechanical.diagnostics['previous_reference_Cu_volume_without_2pi']
        else:baseline=self.initial_volume
        self.engine.commit(trial.mechanical)
        W.x.array[:]=trial.energy;W.x.scatter_forward()
        self.initial_volume=baseline
        self.records.append({'time_s':trial.time_s,'step':trial.step,**trial.mechanical.diagnostics})
