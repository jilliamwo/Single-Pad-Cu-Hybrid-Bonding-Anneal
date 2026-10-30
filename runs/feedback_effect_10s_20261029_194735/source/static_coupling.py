"""Static mechanics inputs for the phase-field model; serial transfer only.

Reference mechanics z is 500 nm above the phase-field origin. No deformation
push-forward of W is applied: this preserves the existing frozen-field model.
"""
from dataclasses import dataclass
from dataclasses import replace
from pathlib import Path
import numpy as np
from dolfinx import fem,geometry as dg
from mpi4py import MPI
from thermoelastic_solver import ReferenceGeometry,cu_at_250C,sio2_dielectric,solve
from validate_baseline import create_baseline_mesh
from calibrate_baseline import load_ansys


@dataclass
class StaticMechanics:
    radii: np.ndarray
    uz: np.ndarray
    dielectric_radii: np.ndarray
    dielectric_uz: np.ndarray
    sample_energy: object
    metadata: dict
    native_result: object = None

    def surface_height(self,r):
        return .993+np.interp(r,self.radii,self.uz)

    @property
    def dielectric_edge_height(self):
        # Phase-field coordinates: nominal dielectric surface at z=1.
        return 1.+float(self.dielectric_uz[0])


def deduplicate_profile(r,values):
    unique,inverse=np.unique(r,return_inverse=True)
    return unique,np.bincount(inverse,weights=values)/np.bincount(inverse)


def evaluate_native_energy(field,phase_points):
    """Evaluate DG0 on Cu interior; never sample a dielectric-side interface cell."""
    result=np.zeros(len(phase_points))
    inside=(phase_points[:,0]>=0)&(phase_points[:,0]<=.5)&(phase_points[:,1]>=0)&(phase_points[:,1]<=.993)
    if not np.any(inside):return result
    points=np.zeros((int(inside.sum()),3));points[:,:2]=phase_points[inside,:2]
    # Translation plus one-sided boundary evaluation; no extrapolation into void.
    points[:,0]=np.clip(points[:,0],1e-10,.5-1e-10)
    points[:,1]=np.clip(points[:,1]+1.,1.+1e-10,1.993-1e-10)
    mesh=field.function_space.mesh
    tree=dg.bb_tree(mesh,mesh.topology.dim)
    collisions=dg.compute_colliding_cells(mesh,dg.compute_collisions_points(tree,points),points)
    cells=np.array([collisions.links(i)[0] if len(collisions.links(i)) else -1 for i in range(len(points))],dtype=np.int32)
    if np.any(cells<0):raise RuntimeError('Native energy transfer has unmapped Cu points')
    result[inside]=field.eval(points,cells).ravel()
    return result


def build_static_mechanics(config):
    if MPI.COMM_WORLD.size!=1:
        raise ValueError('The model coupling and surface diagnostics currently require one MPI rank')
    if config.mechanics_source=='none':
        return StaticMechanics(np.array([0.,.5]),np.zeros(2),np.array([.5,1.]),np.zeros(2),
            lambda x:np.zeros(len(x)),{'source':'none','dielectric_reference':'undeformed'})
    if config.mechanics_source=='native':
        from plasticity import BilinearPlasticity
        g=ReferenceGeometry();mesh=create_baseline_mesh(g,config.mechanics_nr,config.mechanics_nz)
        cu=replace(cu_at_250C(),T_ref_K=config.reference_C+273.15)
        dielectric=replace(sio2_dielectric(),T_ref_K=config.reference_C+273.15)
        result=solve(None,mesh,cu,config.temperature_C+273.15,g.L0_m,.5,.993,geometry=g,
                     material_dielectric=dielectric,plasticity_cu=BilinearPlasticity() if config.plasticity else None,
                     load_steps=config.thermal_steps)
        # Recover reference radial coordinates and vertical displacements on the
        # dielectric surface directly from the normalized displacement Function.
        V=result.u_dimless_field.function_space;xyz=V.tabulate_dof_coordinates()
        values=result.u_dimless_field.x.array.reshape(-1,2)
        mask=np.isclose(xyz[:,1],2.,atol=1e-10,rtol=0)&(xyz[:,0]>=.5-1e-10)
        dr,du=deduplicate_profile(xyz[mask,0],values[mask,1])
        return StaticMechanics(result.UY_r_dimless,result.UY_dimless,dr,du,
            lambda x:evaluate_native_energy(result.W_field,x),
            {'source':'native','statistics':result.stats,'L0_nm':500.,'z_origin_offset_nm':500.,
             'transfer':'one-sided DG0 point evaluation into P1; reference coordinates; zero outside reference Cu',
             'Cu_E_Pa':cu.E_Pa,'Cu_nu':cu.nu,'Cu_alpha_1_K':cu.alpha_1_K,
             'SiO2_E_Pa':dielectric.E_Pa,'SiO2_nu':dielectric.nu,'SiO2_alpha_1_K':dielectric.alpha_1_K,
             'reference_C':config.reference_C,'temperature_C':config.temperature_C,
             'yield_stress_Pa':172.3e6 if config.plasticity else None,
             'tangent_modulus_Pa':517e6 if config.plasticity else None},native_result=result)
    a=load_ansys(Path(config.ansys_csv))
    # Match the existing 3 nm Cu-surface sampling tolerance and radial averaging.
    surface=(a[:,1]>=-1e-15)&(a[:,1]<=250e-9+1e-15)&(np.abs(a[:,2]-996.5e-9)<=3e-9)
    if not np.any(surface):raise ValueError('ANSYS CSV has no Cu top-surface samples')
    radii,uz=deduplicate_profile(a[surface,1]/500e-9,a[surface,4]/500e-9)
    wall=(a[:,1]>=250e-9-1e-15)&(a[:,1]<=500e-9+1e-15)&(np.abs(a[:,2]-1000e-9)<1e-12)
    if not np.any(wall):raise ValueError('ANSYS CSV lacks the dielectric top surface for a deformed reference')
    dr,du=deduplicate_profile(a[wall,1]/500e-9,a[wall,4]/500e-9)
    points=np.column_stack((a[:,1]/500e-9,(a[:,2]-500e-9)/500e-9))
    # Preserve the former source/target mask following the displaced initial top.
    source_mask=(points[:,0]>=0)&(points[:,0]<=.5)&(points[:,1]>=0)&(points[:,1]<=.993+np.interp(points[:,0],radii,uz))
    points=points[source_mask];values=np.maximum(a[source_mask,13],0.)
    points,inverse=np.unique(points,axis=0,return_inverse=True)
    values=np.bincount(inverse,weights=values)/np.bincount(inverse)
    if len(points)<3:raise ValueError('Not enough Cu energy samples in ANSYS CSV')
    if config.ansys_interpolation=='linear':
        try:from scipy.interpolate import LinearNDInterpolator
        except ImportError as error:
            raise ImportError('ANSYS linear interpolation requires scipy. Run: conda install -n fenicsx-env -c conda-forge scipy') from error
        interpolator=LinearNDInterpolator(points,values,fill_value=0.)
        sample=lambda q:np.asarray(interpolator(q[:,:2]),dtype=float)
    else:
        # Explicit legacy fallback; blocked queries avoid an enormous N-by-M array.
        def sample(q):
            out=np.zeros(len(q))
            for start in range(0,len(q),128):
                block=q[start:start+128,:2]
                indices=np.argmin(np.sum((block[:,None,:]-points[None,:,:])**2,axis=2),axis=1)
                out[start:start+len(block)]=values[indices]
            return out
    def energy(q):
        out=sample(q)
        inside=(q[:,0]>=0)&(q[:,0]<=.5)&(q[:,1]>=0)&(q[:,1]<=.993+np.interp(q[:,0],radii,uz))
        out[~inside]=0.
        return out
    import hashlib
    return StaticMechanics(radii,uz,dr,du,energy,
        {'source':'ansys','source_file':str(Path(config.ansys_csv).resolve()),
         'source_sha256':hashlib.sha256(Path(config.ansys_csv).read_bytes()).hexdigest(),
         'z_origin_offset_nm':500.,'transfer':config.ansys_interpolation,
         'surface_sampling_tolerance_nm':3.,'reference_C':'not encoded in CSV',
         'plasticity':'not encoded in CSV'})


def transfer_energy(mechanics,V,multiplier=1.):
    W=fem.Function(V,name='static_Cu_elastic_energy_J_m3')
    values=mechanics.sample_energy(V.tabulate_dof_coordinates()[:,:2])
    if len(values)!=len(W.x.array) or not np.all(np.isfinite(values)) or np.min(values)<-1e-6:
        raise ValueError('Invalid energy transfer')
    W.x.array[:]=multiplier*np.maximum(values,0.)
    W.x.scatter_forward()
    return W
