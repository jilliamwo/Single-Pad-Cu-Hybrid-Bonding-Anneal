"""Fast coupling checks and weak-form regression against the archived model."""
from pathlib import Path
import math
import numpy as np
import ufl
from mpi4py import MPI
from dolfinx import fem,mesh as dm
from dolfinx.fem import petsc
from basix.ufl import element,mixed_element
from model_config import smoke_config
from phase_field_solver import residual,kinetic_rate
from static_coupling import evaluate_native_energy


def main():
    config=smoke_config()
    mesh=dm.create_rectangle(MPI.COMM_WORLD,[np.array([0.,0.]),np.array([1.,2.])],[10,20])
    DG=fem.functionspace(mesh,('DG',0));field=fem.Function(DG)
    field.interpolate(lambda x:np.where((x[0]<.5)&(x[1]>1.)&(x[1]<2.),7.,99.))
    points=np.array([[0.,0.],[.5,0.],[.5,.993],[.2,.5],[.2,1.01],[.51,.5],[-.1,.2],[.2,-.1]])
    np.testing.assert_allclose(evaluate_native_energy(field,points),[7,7,7,7,0,0,0,0])
    P1=element('Lagrange',mesh.basix_cell(),1);V=fem.functionspace(mesh,mixed_element([P1,P1]))
    u=fem.Function(V);old=fem.Function(V)
    u.sub(0).interpolate(lambda x:.2+.1*x[0]+.05*x[1]);u.sub(1).interpolate(lambda x:.3*x[0]-.2*x[1])
    old.sub(0).interpolate(lambda x:.3+.02*x[0]);old.sub(1).interpolate(lambda x:.1*x[1])
    W=fem.Function(fem.functionspace(mesh,P1));W.interpolate(lambda x:1e6*(1+x[0]))
    facets=np.sort(dm.locate_entities_boundary(mesh,1,lambda x:np.isclose(x[0],1.)))
    tags=dm.meshtags(mesh,1,facets,np.ones(len(facets),dtype=np.int32))
    ds=ufl.Measure('ds',domain=mesh,subdomain_data=tags);dt=fem.Constant(mesh,.03)
    fresh=residual(config,u,old,W,dt,ds)
    phi,mu=ufl.split(u);phi_old,mu_old=ufl.split(old);q,v=ufl.TestFunctions(V)
    variable=ufl.variable(phi)
    ns=dict(ufl=ufl,phi=phi,mu=mu,phi_old=phi_old,q=q,v=v,radial_coordinate=ufl.SpatialCoordinate(mesh)[0],
            dx=ufl.Measure('dx',domain=mesh),ds=ds,wall_marker=1,dt_constant=dt,
            kinetic_rate_per_s=kinetic_rate(config),contact_line_mobility_factor_ufl=1.,
            delta_phi_old=16*phi_old**2*(1-phi_old)**2,delta_phi=16*phi**2*(1-phi)**2,mu_mid=mu,
            interface_width=config.interface_width_nm/500,active_wall_cosine=config.wall_energy_J_m2/1.7,
            surface_energy_prefactor=3/math.sqrt(2),epsilon_squared_over_L0_squared=(config.interface_width_nm/500)**2,
            df_dphi=ufl.diff(variable**2*(1-variable)**2,variable),mechanical_energy_prefactor=15/8,
            W_dimensionless=config.interface_width_nm*1e-9/1.7*W)
    archive=Path(__file__).resolve().parents[1]/'old'/'ch_pad_stage15_with_curvature.py'
    original=archive.read_text().split('# 12. WEAK FORM',1)[1].split('# 13. NONLINEAR SOLVER',1)[0]
    exec(original,ns)
    vectors=[petsc.assemble_vector(fem.form(F)) for F in (fresh,ns['F'])]
    np.testing.assert_allclose(vectors[0].array,vectors[1].array,rtol=1e-12,atol=1e-13)
    matrices=[petsc.assemble_matrix(fem.form(ufl.derivative(F,u))) for F in (fresh,ns['F'])]
    for matrix in matrices:matrix.assemble()
    for a,b in zip(matrices[0].getValuesCSR(),matrices[1].getValuesCSR()):np.testing.assert_allclose(a,b,rtol=1e-12,atol=1e-13)
    for obj in vectors+matrices:obj.destroy()
    print('PASS: reference-to-phase coordinate translation, one-sided Cu sampling, outside-Cu masking')
    print('PASS: residual and Jacobian agree with archived phase-field weak form')


if __name__=='__main__':main()
