"""Small-strain associative J2 plasticity with linear isotropic hardening.

History is stored at integration points, never updated during Newton iteration.
Temperature paths may start from virgin material or a copied integration-point
history on the same mesh. It is rate independent: load factors are not physical time or creep.
"""
from dataclasses import dataclass
import numpy as np
import ufl
from basix.ufl import quadrature_element
from dolfinx import fem
from dolfinx.fem.petsc import NonlinearProblem
from mpi4py import MPI


@dataclass(frozen=True)
class BilinearPlasticity:
    yield_stress_Pa: float = 172.3e6
    tangent_modulus_Pa: float = 517e6

    def hardening_modulus(self, E):
        if not np.isfinite(self.yield_stress_Pa) or self.yield_stress_Pa <= 0:
            raise ValueError('Yield stress must be finite and positive')
        if not np.isfinite(self.tangent_modulus_Pa) or not 0 <= self.tangent_modulus_Pa < E:
            raise ValueError('Require 0 <= uniaxial tangent modulus < Young modulus')
        # ANSYS input Et is d(sigma)/d(total strain), not d(sigma_y)/dp.
        return E*self.tangent_modulus_Pa/(E-self.tangent_modulus_Pa)


@dataclass
class PlasticState:
    strain: fem.Function
    equivalent_strain: fem.Function

    def snapshot(self):
        return self.strain.x.array.copy(), self.equivalent_strain.x.array.copy()

    def restore(self, snapshot):
        for f, values in zip((self.strain, self.equivalent_strain), snapshot):
            f.x.array[:] = values
            f.x.scatter_forward()


def radial_return(elastic_trial, C, mu, old_p, yield_stress, H, active):
    trial_stress=C(elastic_trial)
    s=ufl.dev(trial_stress)
    # Safeguard the inactive zero-stress branch, including its AD derivative.
    q=ufl.sqrt(ufl.max_value(1., 1.5*ufl.inner(s,s)))
    dp=active*ufl.max_value(q-yield_stress-H*old_p,0.)/(3*mu+H)
    dep=1.5*dp*s/q
    return trial_stress-2*mu*dep, dep, dp


def solve_plastic_path(V, eps, C, mu, theta, cu_mask, dielectric_mask,
                       material_cu, material_dielectric, law_cu, law_dielectric,
                       bcs, load_factors, measure, *, initial_state=None, initial_displacement=None, constituents=None):
    """Continue or initialize a path; commit only converged load increments.

    Uses degree-4 quadrature for state and weak forms. Temperature path starts
    at each material's reference temperature and scales its target increment.
    Factors are absolute, not increments relative to the restart. Input history
    and displacement are copied, so failures never mutate the caller's state.
    """
    mesh=V.mesh
    factors=np.asarray(load_factors,dtype=float)
    if factors.ndim!=1 or not len(factors) or not np.all(np.isfinite(factors)):
        raise ValueError('Supply a nonempty finite sequence of temperature load factors')
    active=0.;Y=0.;H=0.
    for mask,material,law in [(cu_mask,material_cu,law_cu),
                               (dielectric_mask,material_dielectric,law_dielectric)]:
        if law is not None:
            active=active+mask
            Y=Y+mask*law.yield_stress_Pa
            H=H+mask*law.hardening_modulus(material.E_Pa)
    Qp=fem.functionspace(mesh,quadrature_element(mesh.basix_cell(),value_shape=(),degree=4))
    Qe=fem.functionspace(mesh,quadrature_element(mesh.basix_cell(),value_shape=(3,3),degree=4))
    state=PlasticState(fem.Function(Qe,name='plastic_strain'),fem.Function(Qp,name='equivalent_plastic_strain'))
    u=fem.Function(V,name='displacement_dimless')
    if (initial_state is None) != (initial_displacement is None):
        raise ValueError('Continuation requires both plastic history and dimensionless displacement')
    if initial_state is not None:
        for source,target in ((initial_state.strain,state.strain),
                              (initial_state.equivalent_strain,state.equivalent_strain),
                              (initial_displacement,u)):
            if (source.function_space.mesh is not mesh
                    or source.ufl_shape!=target.ufl_shape
                    or source.function_space.ufl_element()!=target.function_space.ufl_element()
                    or source.x.array.shape!=target.x.array.shape):
                raise ValueError('Continuation requires the same mesh and compatible field spaces; no remapping is performed')
            if not np.all(np.isfinite(source.x.array)):
                raise ValueError('Continuation fields must be finite')
            target.x.array[:]=source.x.array
            target.x.scatter_forward()
        if np.any(state.equivalent_strain.x.array<0):
            raise ValueError('Accumulated plastic strain must be nonnegative')
        tensors=state.strain.x.array.reshape(-1,3,3)
        if not np.allclose(tensors,tensors.swapaxes(1,2),atol=1e-12,rtol=0) or not np.allclose(np.trace(tensors,axis1=1,axis2=2),0,atol=1e-10):
            raise ValueError('J2 plastic strain history must be symmetric and traceless')
    factor=fem.Constant(mesh,0.)
    thermal=factor*theta*ufl.Identity(3)
    trial_elastic=eps(u)-thermal-state.strain
    stress,dep,dp=radial_return(trial_elastic,C,mu,state.equivalent_strain,Y,H,active)
    if constituents is not None:
        # Intrinsic Cu return mapping; occupancy weights force, never yield stress.
        cc=constituents
        Ccu=lambda e: 2*material_cu.lame_mu_Pa*e+material_cu.lame_lambda_Pa*ufl.tr(e)*ufl.Identity(3)
        Cd=lambda e: 2*material_dielectric.lame_mu_Pa*e+material_dielectric.lame_lambda_Pa*ufl.tr(e)*ufl.Identity(3)
        thermal=factor*cc['theta_cu']*ufl.Identity(3)
        active=ufl.conditional(ufl.gt(cu_mask,0.),1.,0.)
        Y=law_cu.yield_stress_Pa;H=law_cu.hardening_modulus(material_cu.E_Pa)
        stress_cu,dep,dp=radial_return(eps(u)-thermal-state.strain,Ccu,
            material_cu.lame_mu_Pa,state.equivalent_strain,Y,H,active)
        elastic_d=eps(u)-factor*cc['theta_dielectric']*ufl.Identity(3)
        background=dielectric_mask*Cd(elastic_d)+cc['void']*cc['void_E']*eps(u)
        stress=cu_mask*stress_cu+background
    r=ufl.SpatialCoordinate(mesh)[0]
    residual=(ufl.inner(stress,eps(ufl.TestFunction(V)))/material_cu.E_Pa)*r*measure
    problem=NonlinearProblem(residual,u,bcs=bcs,petsc_options_prefix='j2_path_',
        petsc_options={'snes_type':'newtonls','snes_linesearch_type':'bt',
                       'snes_rtol':1e-9,'snes_atol':1e-11,'snes_stol':1e-12,
                       'snes_max_it':60,'snes_error_if_not_converged':True,
                       'ksp_type':'preonly','pc_type':'lu','pc_factor_mat_solver_type':'mumps',
                       'ksp_error_if_not_converged':True})
    ep_new=fem.Function(Qe);p_new=fem.Function(Qp)
    ep_expr=fem.Expression(state.strain+dep,Qe.element.interpolation_points)
    p_expr=fem.Expression(state.equivalent_strain+dp,Qp.element.interpolation_points)
    final_elastic=eps(u)-thermal-state.strain
    final_stress=C(final_elastic)
    yield_stress=final_stress
    if constituents is not None:
        yield_stress=Ccu(final_elastic)
        final_stress=cu_mask*yield_stress+background
        cc['cu_energy']=cu_mask*.5*ufl.inner(yield_stress,final_elastic)
        cc['energy']=cc['cu_energy']+dielectric_mask*.5*ufl.inner(Cd(elastic_d),elastic_d)+cc['void']*.5*cc['void_E']*ufl.inner(eps(u),eps(u))
    violation_expr=fem.Expression(active*(ufl.sqrt(1.5*ufl.inner(ufl.dev(yield_stress),ufl.dev(yield_stress)))-Y-H*state.equivalent_strain),Qp.element.interpolation_points)
    violation=fem.Function(Qp)
    owned=Qp.dofmap.index_map.size_local*Qp.dofmap.index_map_bs
    history=[]
    for value in factors:
        old_u=u.x.array.copy();old_state=state.snapshot();old_factor=float(factor.value)
        factor.value=float(value)
        try:
            problem.solve()
            if problem.solver.getConvergedReason()<=0:
                raise RuntimeError('Nonlinear mechanics failed to converge')
            u.x.scatter_forward()
            # Evaluate BOTH new fields before changing either old history field.
            ep_new.interpolate(ep_expr);p_new.interpolate(p_expr)
            valid = (np.all(np.isfinite(p_new.x.array)) and np.all(np.isfinite(ep_new.x.array))
                     and not np.any(p_new.x.array < old_state[1]-1e-13))
            if not mesh.comm.allreduce(bool(valid), op=MPI.LAND):
                raise RuntimeError('Plastic state is nonfinite or accumulated strain decreased')
            state.strain.x.array[:]=ep_new.x.array
            state.equivalent_strain.x.array[:]=p_new.x.array
            state.strain.x.scatter_forward();state.equivalent_strain.x.scatter_forward()
            violation.interpolate(violation_expr)
            vmax=mesh.comm.allreduce(float(np.max(violation.x.array[:owned],initial=0.)),op=MPI.MAX)
            pmax=mesh.comm.allreduce(float(np.max(p_new.x.array[:owned],initial=0.)),op=MPI.MAX)
            history.append({'load_factor':float(value),'newton_iterations':int(problem.solver.getIterationNumber()),
                            'snes_reason':int(problem.solver.getConvergedReason()),
                            'max_equivalent_plastic_strain':pmax,'max_yield_violation_Pa':vmax})
        except Exception:
            state.restore(old_state);u.x.array[:]=old_u;u.x.scatter_forward();factor.value=old_factor
            raise
    return u,final_elastic,final_stress,state,history
