"""Axisymmetric surface diffusion with static or experimental staggered energy."""
import csv
import json
import math
from contextlib import ExitStack
import numpy as np
import ufl
from basix.ufl import element, mixed_element
from dolfinx import fem, mesh as dm, geometry, io
from dolfinx.fem.petsc import NonlinearProblem
from mpi4py import MPI
from petsc4py import PETSc
from static_coupling import transfer_energy


def kinetic_rate(config):
    T=config.temperature_C+273.15
    eps=config.interface_width_nm*1e-9
    D=.26e-4*math.exp(-(20.7*4184)/(8.314462618*T))
    mobility=config.diffusivity_factor*D*.30e-9/(1.18e-29*1.380649e-23*T)
    return 3/(2*math.sqrt(2))*1.18e-29*mobility*(1.7*1.18e-29/eps)/(eps*(500e-9)**2)


def residual(config,u,u_old,W,dt,ds):
    phi,mu=ufl.split(u);old,_=ufl.split(u_old)
    q,v=ufl.TestFunctions(u.function_space)
    r=ufl.SpatialCoordinate(u.function_space.mesh)[0]
    dx=ufl.Measure('dx',domain=u.function_space.mesh)
    eps=config.interface_width_nm/500.
    delta_old=16*old**2*(1-old)**2
    delta=16*phi**2*(1-phi)**2
    df=2*phi*(1-phi)*(1-2*phi)
    dwall=-.5*eps*(config.wall_energy_J_m2/1.7)*6*phi*(1-phi)
    return ((phi-old)*q*r*dx+dt*kinetic_rate(config)*delta_old*ufl.inner(ufl.grad(mu),ufl.grad(q))*r*dx
            +mu*v*r*dx-3/math.sqrt(2)*eps**2*ufl.inner(ufl.grad(phi),ufl.grad(v))*r*dx
            -3/math.sqrt(2)*df*v*r*dx+15/8*delta*(config.interface_width_nm*1e-9/1.7)*W*v*r*dx
            -dwall*v*r*ds(1))


class SurfaceSampler:
    def __init__(self,mesh):
        self.r=np.linspace(0,.5,201);self.z=np.linspace(.85,1.1-1e-10,401)
        rr,zz=np.meshgrid(np.clip(self.r,1e-10,.5-1e-10),self.z,indexing='ij')
        self.points=np.column_stack((rr.ravel(),zz.ravel(),np.zeros(rr.size)))
        tree=geometry.bb_tree(mesh,2)
        links=geometry.compute_colliding_cells(mesh,geometry.compute_collisions_points(tree,self.points),self.points)
        self.cells=np.array([links.links(i)[0] if len(links.links(i)) else -1 for i in range(len(self.points))],dtype=np.int32)
        if np.any(self.cells<0):raise RuntimeError('Unmapped surface diagnostic points')

    def sample(self,phi):
        values=phi.eval(self.points,self.cells).reshape(len(self.r),len(self.z))
        height=np.full(len(self.r),np.nan)
        for i,row in enumerate(values):
            crossings=np.flatnonzero((row[:-1]>=.5)&(row[1:]<.5))
            if len(crossings):
                j=crossings[-1]
                height[i]=self.z[j]+(.5-row[j])/(row[j+1]-row[j])*(self.z[j+1]-self.z[j])
        return height


def run_phase_field(config,mechanics,out,*,feedback=None):
    mesh=dm.create_rectangle(MPI.COMM_WORLD,[np.array([0.,0.]),np.array([.5,1.1])],
                             [config.phase_nr,config.phase_nz],cell_type=dm.CellType.triangle)
    P1=element('Lagrange',mesh.basix_cell(),1)
    V=fem.functionspace(mesh,P1);mixed=fem.functionspace(mesh,mixed_element([P1,P1]))
    u=fem.Function(mixed);old=fem.Function(mixed)
    u.sub(0).interpolate(lambda x:.5*(1-np.tanh((x[1]-mechanics.surface_height(x[0]))/(np.sqrt(2)*config.interface_width_nm/500))))
    old.x.array[:]=u.x.array;old.x.scatter_forward()
    phi=fem.Function(V,name='phi');phi.interpolate(u.sub(0))
    sampler=SurfaceSampler(mesh)
    initial_surface=sampler.sample(phi)
    W=(transfer_energy(mechanics,V,config.mechanical_energy_multiplier) if feedback is None else
       feedback.initialize(V,sampler.r,initial_surface,config.mechanical_energy_multiplier,phase=phi))
    frozen=W.x.array.copy()
    coords=V.tabulate_dof_coordinates()
    np.savetxt(out/'coupling_samples.csv',np.column_stack((coords[:,0]*500,coords[:,1]*500,W.x.array)),delimiter=',',header='r_nm,z_phase_nm,W_J_m3',comments='')
    facets=np.sort(dm.locate_entities_boundary(mesh,1,lambda x:np.isclose(x[0],.5)))
    tags=dm.meshtags(mesh,1,facets,np.ones(len(facets),dtype=np.int32))
    ds=ufl.Measure('ds',domain=mesh,subdomain_data=tags)
    dt=fem.Constant(mesh,PETSc.ScalarType(config.dt_s))
    F=residual(config,u,old,W,dt,ds)
    def make_problem():
        return NonlinearProblem(F,u,petsc_options_prefix='model_',petsc_options={
            'snes_type':'newtonls','snes_linesearch_type':'bt','snes_rtol':1e-8,'snes_atol':1e-10,
            'snes_stol':1e-10,'snes_max_it':120,'ksp_type':'preonly','pc_type':'lu'})
    problem=make_problem()
    r=ufl.SpatialCoordinate(mesh)[0]
    massform=fem.form(phi*r*ufl.dx)
    initial_mass=float(fem.assemble_scalar(massform))
    initial_surface=sampler.sample(phi)
    t=0.;step=0;rejected=0;easy=0;next_measure=config.measure_interval_s;next_save=config.save_interval_s
    adaptive_dt=config.dt_s
    solver_errors=[];recoveries=0;rebuild_pending=False
    with ExitStack() as stack:
        historyfile=stack.enter_context((out/'time_history.csv').open('w',newline=''))
        history=csv.writer(historyfile)
        history.writerow(['time_s','mass_axisymmetric_normalized','relative_mass_drift','phi_min','phi_max','center_z_phase_nm','edge_z_phase_nm','center_relative_to_dielectric_edge_nm','missing_surface_fraction'])
        stepsfile=stack.enter_context((out/'step_history.csv').open('w',newline=''));steps=csv.writer(stepsfile)
        steps.writerow(['step','time_before_s','dt_s','accepted','snes_reason','iterations','petsc_error_code'])
        surfacefile=stack.enter_context((out/'surface_history.csv').open('w',newline=''))
        surfaces=csv.writer(surfacefile)
        surfaces.writerow(['time_s','r_nm','z_phase_nm'])
        xdmf=None
        if config.write_fields:
            xdmf=stack.enter_context(io.XDMFFile(mesh.comm,str(out/'phi.xdmf'),'w'));xdmf.write_mesh(mesh)
            with io.XDMFFile(mesh.comm,str(out/'static_energy.xdmf'),'w') as energyfile:
                energyfile.write_mesh(mesh);energyfile.write_function(W)
        def measure():
            phi.interpolate(u.sub(0));height=sampler.sample(phi)
            mass=float(fem.assemble_scalar(massform))
            row=[t,mass,(mass-initial_mass)/initial_mass,float(phi.x.array.min()),float(phi.x.array.max()),height[0]*500,height[-1]*500,(height[0]-mechanics.dielectric_edge_height)*500,float(np.mean(~np.isfinite(height)))]
            history.writerow(row);historyfile.flush()
            surfaces.writerows(zip(np.full(len(height),t),sampler.r*500,height*500))
            surfacefile.flush()
            np.savetxt(out/'final_surface.csv',np.column_stack((sampler.r*500,height*500)),delimiter=',',header='r_nm,z_phase_nm',comments='')
            (out/'progress.json').write_text(json.dumps({'time_s':t,'accepted_steps':step,'rejected_steps':rejected,'solver_rebuilds':recoveries,'relative_mass_drift':row[2],'missing_surface_fraction':row[-1]},indent=2))
            print(f't = {t:g} s | mass drift = {row[2]:.3e} | center relative to dielectric edge = {row[7]:.6g} nm',flush=True)
            return row
        np.savetxt(out/'initial_surface.csv',np.column_stack((sampler.r*500,initial_surface*500,mechanics.surface_height(sampler.r)*500)),delimiter=',',header='r_nm,z_phase_nm,analytic_initial_z_phase_nm',comments='')
        last=measure()
        if xdmf:xdmf.write_function(phi,t)
        while t<config.duration_s-1e-12:
            attempt_dt=min(adaptive_dt,config.duration_s-t);dt.value=attempt_dt
            u.x.array[:]=old.x.array;u.x.scatter_forward()
            failure=None;error_code=None;coupling_error=None;coupling_trial=None
            try:
                problem.solve()
                reason=int(problem.solver.getConvergedReason());iterations=int(problem.solver.getIterationNumber())
                accepted=reason>0 and np.all(np.isfinite(u.x.array))
            except PETSc.Error as error:
                reason=int(problem.solver.getConvergedReason())
                iterations=int(problem.solver.getIterationNumber())
                accepted=False;failure=str(error);error_code=int(error.ierr)
                solver_errors.append({'time_s':t,'dt_s':attempt_dt,'error_code':error_code,
                                      'message':failure,'after_rebuild':rebuild_pending})
                (out/'solver_errors.json').write_text(json.dumps(solver_errors,indent=2))
            if accepted and feedback is not None:
                try:
                    phi.interpolate(u.sub(0))
                    coupling_trial=feedback.propose(sampler.sample(phi),t+attempt_dt,step+1,phase=phi)
                except Exception as error:
                    coupling_error=error;accepted=False
                    with (out/'coupling_failures.log').open('a') as log:
                        log.write(f't={t:g}, dt={attempt_dt:g}, {type(error).__name__}: {error}\n')
            steps.writerow([step+1,t,attempt_dt,int(accepted),reason,iterations,error_code]);stepsfile.flush()
            if not accepted:
                u.x.array[:]=old.x.array;u.x.scatter_forward();rejected+=1;easy=0
                # Capture accepted-state diagnostics before retrying. A failed Newton
                # iterate must never be reported as a physical state.
                last=measure()
                if coupling_error is not None and isinstance(coupling_error,ValueError):
                    raise RuntimeError(f'Feedback geometry/history policy rejected at t={t:g}s: {coupling_error}') from coupling_error
                if failure is not None:
                    if rebuild_pending:
                        np.savez_compressed(out/'last_accepted_state.npz',time_s=t,
                                            mixed_dofs=old.x.array,coordinates=mesh.geometry.x)
                        raise RuntimeError(f'PETSc internal error repeated at t={t:g}s after solver rebuild; '
                                           f'see solver_errors.json. {failure}')
                    # Do not repeatedly halve dt for a corrupted solver object.
                    # A fresh NonlinearProblem owns fresh SNES, matrix and vectors.
                    problem=None
                    problem=make_problem();recoveries+=1;rebuild_pending=True
                    print(f'PETSc error {error_code} at t={t:g}s; rebuilding solver and retrying dt={attempt_dt:g}s',flush=True)
                    continue
                rebuild_pending=False
                if attempt_dt/2<config.min_dt_s:
                    np.savez_compressed(out/'last_accepted_state.npz',time_s=t,
                                        mixed_dofs=old.x.array,coordinates=mesh.geometry.x)
                    raise RuntimeError(f'Step failed at t={t:g}s, dt={attempt_dt:g}s; minimum dt reached. '
                                       f'{coupling_error if coupling_error is not None else f"Phase SNES={reason}"}')
                adaptive_dt=attempt_dt/2;print(f'Retrying at dt={adaptive_dt:g}s',flush=True);continue
            rebuild_pending=False
            if feedback is not None:
                feedback.commit(coupling_trial,W)
                (out/'coupling_history.json').write_text(json.dumps(feedback.records,indent=2))
                np.savetxt(out/'latest_energy_samples.csv',np.column_stack((coords[:,0]*500,coords[:,1]*500,W.x.array)),delimiter=',',header='r_nm,z_phase_nm,W_J_m3',comments='')
            t+=attempt_dt;step+=1;old.x.array[:]=u.x.array;old.x.scatter_forward()
            (out/'progress.json').write_text(json.dumps({'time_s':t,'accepted_steps':step,'rejected_steps':rejected,'solver_rebuilds':recoveries,'last_dt_s':attempt_dt},indent=2))
            easy=easy+1 if iterations<=8 else 0
            if easy>=20:adaptive_dt=min(config.max_dt_s,adaptive_dt*1.25);easy=0
            final=t>=config.duration_s-1e-12
            if t>=next_measure-1e-12 or final:
                last=measure()
                while next_measure<=t+1e-12:next_measure+=config.measure_interval_s
            if t>=next_save-1e-12 or final:
                phi.interpolate(u.sub(0))
                if xdmf:xdmf.write_function(phi,t)
                while next_save<=t+1e-12:next_save+=config.save_interval_s
        if feedback is None and not np.array_equal(frozen,W.x.array):raise RuntimeError('Frozen mechanical energy changed')
    return {'final_time_s':t,'accepted_steps':step,'rejected_steps':rejected,'solver_rebuilds':recoveries,'relative_mass_drift':last[2],
            'phi_min':last[3],'phi_max':last[4],'center_relative_to_dielectric_edge_nm':last[7],
            'missing_surface_fraction':last[-1],'kinetic_rate_per_s':kinetic_rate(config),
            'coupling':('static one-way; no mechanics feedback, contact, or transient heat equation' if feedback is None else
                        'experimental staggered feedback; frozen column coordinates'),
            'mechanical_updates':0 if feedback is None else len(feedback.records),
            'height_reference':'initial deformed dielectric edge' }
