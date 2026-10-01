"""Single-pad axisymmetric thermoelasticity on a normalized reference mesh.

Coordinates and the solved displacement are divided by L0; u_field is
converted to metres for output. The reference material map identifies the
SiO2 base and annulus independently of phi. Elastic or optional J2 elastoplastic pre-contact mechanics.
Stage 15 uses a different origin: z_stage15 = z_reference - base_thickness/L0.
ANSYS parity has not yet been established.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Optional

import numpy as np
from mpi4py import MPI

import ufl
from basix.ufl import element
from dolfinx import default_real_type, fem
from dolfinx.fem import Function, functionspace, Constant
from dolfinx.fem.petsc import LinearProblem


# ---------------------------------------------------------------------------
# 1.  MATERIAL PARAMETER DATACLASSES
# ---------------------------------------------------------------------------

@dataclass
class MaterialParams:
    """
    All dimensional SI material parameters for one material zone.

    Attributes
    ----------
    name : str
        Human-readable label (e.g. 'Cu_250C', 'SiO2', 'Al_250C').
    E_Pa : float
        Young's modulus [Pa].
    nu : float
        Poisson's ratio [-].
    alpha_1_K : float
        Linear CTE [1/K].
    k_W_mK : float
        Thermal conductivity [W/(m·K)].
    rho_cp_J_m3K : float
        Volumetric heat capacity [J/(m³·K)].
    T_ref_K : float
        Stress-free (reference) temperature [K].
        ΔT = T_anneal - T_ref drives the thermoelastic strain.
    gamma_surface_J_m2 : float
        Metal/void surface energy [J/m²].  Used to compute Δγ_w.
    gamma_liner_void_J_m2 : float
        Liner–void interface energy [J/m²].  None = not yet measured.
    gamma_metal_liner_J_m2 : float
        Metal–liner interface energy [J/m²].  None = not yet measured.
    """
    name: str
    E_Pa: float
    nu: float
    alpha_1_K: float
    k_W_mK: float
    rho_cp_J_m3K: float
    T_ref_K: float
    gamma_surface_J_m2: float
    gamma_liner_void_J_m2: Optional[float] = None
    gamma_metal_liner_J_m2: Optional[float] = None

    # ---- derived elastic constants ----------------------------------------
    @property
    def lame_lambda_Pa(self) -> float:
        """First Lamé constant λ [Pa]."""
        return (self.E_Pa * self.nu) / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))

    @property
    def lame_mu_Pa(self) -> float:
        """Shear modulus / second Lamé constant μ [Pa]."""
        return self.E_Pa / (2.0 * (1.0 + self.nu))

    # ---- wall energy -------------------------------------------------------
    def wall_energy_difference_J_m2(self) -> float:
        """
        Δγ_w = γ_liner-void − γ_metal-liner  [J/m²]

        Derived from three independent interface energies.
        Consistent with Young's equation:
            Δγ_w = γ_surface · cos(θ_eq)

        Raises ValueError if either interface energy is not yet set.
        Prints the implied contact angle as a diagnostic.
        """
        if self.gamma_liner_void_J_m2 is None or self.gamma_metal_liner_J_m2 is None:
            raise ValueError(
                f"[{self.name}] gamma_liner_void_J_m2 and gamma_metal_liner_J_m2 "
                "must both be set to compute Δγ_w.  "
                "Derive them from DFT W_sep or solid-state wetting experiments."
            )
        delta = self.gamma_liner_void_J_m2 - self.gamma_metal_liner_J_m2
        cos_theta = delta / self.gamma_surface_J_m2
        if abs(cos_theta) > 1.0:
            raise ValueError(
                f"[{self.name}] Unphysical contact angle: cos(θ) = {cos_theta:.4f}.  "
                "Check interface energy inputs — Young's equation requires "
                "|γ_liner-void − γ_metal-liner| ≤ γ_surface."
            )
        theta_deg = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
        print(
            f"[{self.name}] Δγ_w = {delta:+.4f} J/m²  →  "
            f"θ_eq = {theta_deg:.2f}°  (cos θ = {cos_theta:+.4f})"
        )
        return delta

    def wall_energy_difference_from_angle(self, theta_eq_deg: float) -> float:
        """
        Fallback: compute Δγ_w from a known contact angle θ_eq [deg].

        Use this when γ_liner-void and γ_metal-liner are not individually
        measured but the equilibrium contact angle is known from TEM.
        Sets gamma_liner_void and gamma_metal_liner to None-safe values.
        """
        theta_rad = np.radians(theta_eq_deg)
        delta = self.gamma_surface_J_m2 * np.cos(theta_rad)
        print(
            f"[{self.name}] Contact-angle fallback: θ = {theta_eq_deg:.2f}°  →  "
            f"Δγ_w = {delta:+.4f} J/m²"
        )
        return delta


# ---- Pre-built material tables -------------------------------------------

def cu_at_250C() -> MaterialParams:
    """Cu elastic constants inferred from SiO2 baseline CSV; thermal inputs provisional."""
    return MaterialParams(
        name="Cu_250C",
        E_Pa=129.0e9,
        nu=0.340,
        alpha_1_K=17.0e-6,
        k_W_mK=380.0,
        rho_cp_J_m3K=3.45e6,
        T_ref_K=293.15,        # 20 °C stress-free reference
        gamma_surface_J_m2=1.70,
        # Cu/TaN values: constrained by Hooge TEM (θ=112.5°)
        # γ_liner-void and γ_metal-liner individually underdetermined —
        # set via wall_energy_difference_from_angle(112.5) until DFT values available.
        gamma_liner_void_J_m2=None,
        gamma_metal_liner_J_m2=None,
    )

def sio2_dielectric() -> MaterialParams:
    """SiO₂ dielectric — surrounding material."""
    return MaterialParams(
        name="SiO2",
        E_Pa=73.0e9,
        nu=0.170,
        alpha_1_K=0.55e-6,
        k_W_mK=1.4,
        rho_cp_J_m3K=1.64e6,
        T_ref_K=293.15,
        gamma_surface_J_m2=1.0,   # not used for wall BC
    )

def al_at_200C() -> MaterialParams:
    """Al pad at 200 °C — placeholder for Al-Al bonding extension."""
    return MaterialParams(
        name="Al_200C",
        E_Pa=65.0e9,
        nu=0.330,
        alpha_1_K=23.4e-6,
        k_W_mK=205.0,
        rho_cp_J_m3K=2.43e6,
        T_ref_K=293.15,
        gamma_surface_J_m2=1.05,
        gamma_liner_void_J_m2=None,   # TODO: DFT / Kaloyeros-style experiment
        gamma_metal_liner_J_m2=None,
    )

def indium_at_150C() -> MaterialParams:
    """In pad at 150 °C — quantum application bonding."""
    return MaterialParams(
        name="In_150C",
        E_Pa=10.6e9,
        nu=0.450,
        alpha_1_K=25.0e-6,
        k_W_mK=81.6,
        rho_cp_J_m3K=1.27e6,
        T_ref_K=293.15,
        gamma_surface_J_m2=0.56,
        gamma_liner_void_J_m2=None,   # TODO: AIMD required
        gamma_metal_liner_J_m2=None,
    )

def niobium_at_250C() -> MaterialParams:
    """Nb pad at 250 °C — quantum application bonding."""
    return MaterialParams(
        name="Nb_250C",
        E_Pa=105.0e9,
        nu=0.400,
        alpha_1_K=7.3e-6,
        k_W_mK=53.7,
        rho_cp_J_m3K=2.13e6,
        T_ref_K=293.15,
        gamma_surface_J_m2=2.655,
        gamma_liner_void_J_m2=None,   # TODO: DFT required
        gamma_metal_liner_J_m2=None,
    )


# ---------------------------------------------------------------------------
# 2.  MESH NODE TRACKER
# ---------------------------------------------------------------------------

class MeshNodeTracker:
    """
    Holds the (r, z) positions of every DOF in a scalar function space
    and exposes fast spatial queries.

    This is the structured mesh information store that phase_field_solver.py
    and thermoelastic_solver.py both use to keep track of where data lives,
    enabling diagnostics, masking, and targeted output without repeated
    FEniCSx coordinate queries.

    Parameters
    ----------
    mesh : dolfinx.mesh.Mesh
        The shared FEniCSx mesh.
    L0_m : float
        Physical length scale [m].  All stored coordinates are dimensionless
        (divided by L0_m).  Physical values = coord * L0_m.

    Attributes
    ----------
    r_dimless : ndarray, shape (N,)
        Radial coordinates of all DOFs in dimensionless units.
    z_dimless : ndarray, shape (N,)
        Vertical coordinates of all DOFs in dimensionless units.
    r_m : ndarray
        Radial coordinates in metres.
    z_m : ndarray
        Vertical coordinates in metres.
    n_dofs : int
        Total number of scalar DOFs.
    stats : dict
        Summary statistics written at construction time.
    """

    def __init__(self, mesh, L0_m: float):
        self.mesh = mesh
        self.L0_m = L0_m

        P1 = element("Lagrange", mesh.basix_cell(), 1, dtype=default_real_type)
        V_scalar = functionspace(mesh, P1)
        coords = V_scalar.tabulate_dof_coordinates()   # shape (N, 3)

        self.r_dimless = coords[:, 0].copy()
        self.z_dimless = coords[:, 1].copy()
        self.r_m = self.r_dimless * L0_m
        self.z_m = self.z_dimless * L0_m
        self.n_dofs = len(self.r_dimless)

        self.stats = {
            "n_dofs": self.n_dofs,
            "r_min_nm": float(self.r_m.min()) * 1e9,
            "r_max_nm": float(self.r_m.max()) * 1e9,
            "z_min_nm": float(self.z_m.min()) * 1e9,
            "z_max_nm": float(self.z_m.max()) * 1e9,
        }

    # ---- convenience masks ------------------------------------------------

    def cu_pad_mask(
        self,
        pad_radius_dimless: float,
        pad_top_dimless: float,
        phi: Optional[Function] = None,
        phi_threshold: float = 0.5,
    ) -> np.ndarray:
        """
        Boolean mask for DOFs inside the Cu pad.

        If phi is provided, uses the phase-field value as the primary
        indicator (phi > phi_threshold), falling back to the geometric
        bounding box.  This is the correct mask for the coupled solver
        because the Cu domain evolves during the anneal.

        Parameters
        ----------
        pad_radius_dimless, pad_top_dimless : float
            Geometric bounds (dimensionless).
        phi : Function, optional
            Current phase-field.  When supplied, the geometric mask is
            ANDed with the phi-based mask so void regions are excluded.
        phi_threshold : float
            Phase-field cutoff for Cu vs void.
        """
        geom_mask = (
            (self.r_dimless >= 0.0)
            & (self.r_dimless <= pad_radius_dimless)
            & (self.z_dimless >= 0.0)
            & (self.z_dimless <= pad_top_dimless)
        )
        if phi is None:
            return geom_mask

        # Collapse phi to scalar DOFs and threshold
        V_phi = phi.function_space
        phi_vals = phi.x.array
        phi_mask = phi_vals > phi_threshold

        # phi lives on its own DOF ordering — build r/z for it
        phi_coords = V_phi.tabulate_dof_coordinates()
        phi_r = phi_coords[:, 0]
        phi_z = phi_coords[:, 1]

        # Nearest-neighbour lookup from tracker DOFs onto phi DOFs
        # (both live on P1 of the same mesh, so this is a 1-1 match)
        phi_combined = np.zeros(self.n_dofs, dtype=bool)
        for i in range(self.n_dofs):
            dr = phi_r - self.r_dimless[i]
            dz = phi_z - self.z_dimless[i]
            j = int(np.argmin(dr**2 + dz**2))
            phi_combined[i] = phi_mask[j]

        return geom_mask & phi_combined

    def surface_mask(
        self,
        surface_z_dimless: float,
        tol_dimless: float = 0.006,   # ≈ 3 nm at L0=500 nm
    ) -> np.ndarray:
        """Boolean mask for DOFs within tol of a given z height."""
        return np.abs(self.z_dimless - surface_z_dimless) <= tol_dimless

    def radial_profile_at_z(
        self,
        values: np.ndarray,
        surface_z_dimless: float,
        tol_dimless: float = 0.006,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Extract a radial profile of an array at a given z height.

        Returns sorted (r_nm, values) arrays restricted to DOFs near
        the specified vertical position.
        """
        mask = self.surface_mask(surface_z_dimless, tol_dimless)
        r_sorted_idx = np.argsort(self.r_dimless[mask])
        r_nm = self.r_m[mask][r_sorted_idx] * 1e9
        vals = values[mask][r_sorted_idx]
        return r_nm, vals

    def report(self) -> str:
        s = self.stats
        return (
            f"MeshNodeTracker: {s['n_dofs']} DOFs  "
            f"r=[{s['r_min_nm']:.1f}, {s['r_max_nm']:.1f}] nm  "
            f"z=[{s['z_min_nm']:.1f}, {s['z_max_nm']:.1f}] nm"
        )


# ---------------------------------------------------------------------------
# 3.  THERMOELASTIC SOLVE RESULT
# ---------------------------------------------------------------------------

@dataclass
class ThermoelasticResult:
    """
    All outputs from one thermoelastic solve, ready for consumption
    by phase_field_solver.py.

    Attributes
    ----------
    W_field : Function
        Cu-localized strain-energy density W [J/m³], represented in DG0.
        Transfer to the phase-field mesh requires an explicit coordinate map.
    UY_dimless : ndarray
        Vertical displacement along the Cu top surface in dimensionless
        FEniCSx coordinates.  Shape: (n_surface_points,) paired with
        UY_r_dimless, on the reference mechanics mesh.
    UY_r_dimless : ndarray
        Radial positions (dimensionless) corresponding to UY_dimless.
    u_field : Function
        Full (r, z) displacement vector [m] for diagnostics/ParaView.
    sigma_components : dict of str → Function
        Individual stress components [Pa]:
          'sigma_rr', 'sigma_zz', 'sigma_tt' (hoop), 'sigma_rz'
    stats : dict
        Solver timing, peak values, convergence info.
    """
    W_field: Function
    UY_dimless: np.ndarray
    UY_r_dimless: np.ndarray
    u_field: Function
    sigma_components: dict
    stats: dict
    W_full_field: Function
    W_ansys_like_field: Function
    u_dimless_field: Function
    cu_surface_m: np.ndarray
    dielectric_surface_m: np.ndarray
    plastic_state: object = None
    equivalent_plastic_strain_field: Optional[Function] = None


# ---------------------------------------------------------------------------
# 4.  THERMAL SOLVER
# ---------------------------------------------------------------------------

def solve_thermal(
    mesh,
    tracker: MeshNodeTracker,
    material_cu: MaterialParams,
    material_dielectric: MaterialParams,
    T_anneal_K: float,
    T_initial_K: float = 293.15,
    pad_radius_dimless: float = 0.5,
    pad_height_dimless: float = 1.0,
    ramp_time_s: float = 300.0,
    steady_state: bool = True,
) -> Function:
    """Return uniform anneal temperature; transient heating is not implemented.

    Legacy thermal arguments are retained for callers, but only the uniform
    temperature path is supported during elastic baseline validation.
    """
    if not steady_state:
        raise NotImplementedError("Transient thermal scaling/material map is not validated")
    P1 = element("Lagrange", mesh.basix_cell(), 1, dtype=default_real_type)
    V_T = functionspace(mesh, P1)
    T_field = Function(V_T, name="temperature_K")

    if steady_state:
        T_field.x.array[:] = T_anneal_K
        T_field.x.scatter_forward()
        print(
            f"[Thermal] Steady-state: T = {T_anneal_K:.2f} K "
            f"({T_anneal_K - 273.15:.1f} °C) uniform"
        )
        return T_field

# ---------------------------------------------------------------------------
# 5.  STRUCTURAL (THERMOELASTIC) SOLVER  — AXISYMMETRIC 4-COMPONENT FORM
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ReferenceGeometry:
    L0_nm: float = 500.0
    outer_radius_nm: float = 500.0
    pad_radius_nm: float = 250.0
    base_thickness_nm: float = 500.0
    pad_height_nm: float = 500.0
    recess_nm: float = 3.5

    def __post_init__(self):
        values = (self.L0_nm, self.outer_radius_nm, self.pad_radius_nm,
                  self.base_thickness_nm, self.pad_height_nm, self.recess_nm)
        if not np.all(np.isfinite(values)) or min(values[:5]) <= 0:
            raise ValueError("Geometry lengths must be finite and positive")
        if not 0 <= self.recess_nm < self.pad_height_nm:
            raise ValueError("Recess must lie within the pad height")
        if self.pad_radius_nm >= self.outer_radius_nm:
            raise ValueError("Dielectric outer radius must exceed Cu radius")

    @property
    def L0_m(self):
        return self.L0_nm * 1e-9

    @property
    def cu_top(self):
        return (self.base_thickness_nm + self.pad_height_nm - self.recess_nm) / self.L0_nm

    @property
    def top(self):
        return (self.base_thickness_nm + self.pad_height_nm) / self.L0_nm


@dataclass(frozen=True)
class MechanicalBCs:
    axis_ur_fixed: bool = True
    bottom_uz_fixed: bool = True
    bottom_ur_fixed: bool = True
    outer_ur_fixed: bool = True
    outer_uz_fixed: bool = False


def axisymmetric_strain(u, r):
    # Regular axis limit; r is dimensionless and ur is constrained on the axis.
    hoop = ufl.conditional(ufl.gt(r, 0.0), u[0] / r, u[0].dx(0))
    shear = 0.5 * (u[0].dx(1) + u[1].dx(0))
    return ufl.as_tensor([[u[0].dx(0), shear, 0],
                          [shear, u[1].dx(1), 0], [0, 0, hoop]])


def _linear_solve(a, L, bcs=(), prefix="mechanics_"):
    problem = LinearProblem(
        a, L, bcs=list(bcs), petsc_options_prefix=prefix,
        petsc_options={"ksp_type": "preonly", "pc_type": "lu",
                       "pc_factor_mat_solver_type": "mumps",
                       "ksp_error_if_not_converged": True})
    solution = problem.solve()
    solution.x.scatter_forward()
    return solution


def solve_structural(mesh, tracker, material_cu, T_field,
                     pad_radius_dimless, pad_height_dimless, L0_m,
                     phi=None, contact_state="pre", p_contact=None,
                     surface_z_dimless=None, surface_tol_dimless=0.006,
                     *, geometry=None, material_dielectric=None,
                     boundary_conditions=None, residual_void_E_Pa=1e3,
                     plasticity_cu=None, plasticity_dielectric=None, load_steps=20, load_factors=None):
    """Solve on the full base/annulus reference domain.

    Legacy pad_height_dimless denotes Cu thickness after recess. Explicit
    geometry is required to prevent silently treating a Cu-only mesh as the
    full mechanical domain. phi must be a scalar Function on this mesh.
    Void uses residual stiffness; this is not a contact model.
    Optional plastic laws integrate a virgin-state thermal load path; phi-driven
    evolving plasticity and external history continuation are not supported yet.
    W_field is Cu-localized; W_full_field retains energy in every material.
    """
    if geometry is None:
        raise ValueError("Pass ReferenceGeometry and a full base/annulus mesh")
    g = geometry
    if contact_state != "pre" or p_contact is not None:
        raise NotImplementedError("Contact traction and bonding are not implemented")
    if not np.isfinite(residual_void_E_Pa) or residual_void_E_Pa <= 0:
        raise ValueError("Residual void stiffness must be positive and finite")
    if not np.isclose(L0_m, g.L0_m, rtol=1e-12, atol=0):
        raise ValueError("L0_m must match geometry")
    if not np.isclose(pad_radius_dimless, g.pad_radius_nm/g.L0_nm) or not np.isclose(
            pad_height_dimless, (g.pad_height_nm-g.recess_nm)/g.L0_nm):
        raise ValueError("Legacy pad dimensions must match geometry")
    if phi is not None and (phi.function_space.mesh is not mesh or phi.ufl_shape):
        raise ValueError("phi must be scalar and live on the mechanics mesh")
    for axis, expected in ((0, g.outer_radius_nm/g.L0_nm), (1, g.top)):
        local = mesh.geometry.x[:,axis]
        lo = mesh.comm.allreduce(float(np.min(local, initial=np.inf)), op=MPI.MIN)
        hi = mesh.comm.allreduce(float(np.max(local, initial=-np.inf)), op=MPI.MAX)
        if not np.isclose(lo, 0., atol=1e-12) or not np.isclose(hi, expected):
            raise ValueError("Mesh bounds must match the full reference geometry")
    dielectric = material_dielectric or sio2_dielectric()
    bc_config = boundary_conditions or MechanicalBCs()
    for material in (material_cu, dielectric):
        if material.E_Pa <= 0 or not -1 < material.nu < 0.5:
            raise ValueError("Require E > 0 and -1 < nu < 0.5")
    t_start = time.perf_counter()
    x = ufl.SpatialCoordinate(mesh)
    r, z = x[0], x[1]
    base = g.base_thickness_nm/g.L0_nm
    radius = g.pad_radius_nm/g.L0_nm
    outer = g.outer_radius_nm/g.L0_nm
    # Fixed reference dielectric identity; both dielectric regions deform.
    dielectric_mask = ufl.conditional(ufl.Or(ufl.lt(z, base), ufl.gt(r, radius)), 1., 0.)
    if phi is None:
        h = ufl.conditional(ufl.le(z, g.cu_top), 1., 0.)
    else:
        clipped = ufl.min_value(1., ufl.max_value(0., phi))
        h = clipped**2 * (3. - 2.*clipped)
    cu = (1. - dielectric_mask) * h
    void = (1. - dielectric_mask) * (1. - h)
    lam = cu*material_cu.lame_lambda_Pa + dielectric_mask*dielectric.lame_lambda_Pa
    mu = (cu*material_cu.lame_mu_Pa + dielectric_mask*dielectric.lame_mu_Pa
          + void*residual_void_E_Pa/2.)  # residual void nu = 0
    theta = (cu*material_cu.alpha_1_K*(T_field-material_cu.T_ref_K)
             + dielectric_mask*dielectric.alpha_1_K*(T_field-dielectric.T_ref_K))
    V = functionspace(mesh, element("Lagrange", mesh.basix_cell(), 1, shape=(2,)))
    Q = functionspace(mesh, element("DG", mesh.basix_cell(), 0))
    trial, test = ufl.TrialFunction(V), ufl.TestFunction(V)
    eps = lambda u: axisymmetric_strain(u, r)
    C = lambda e: 2*mu*e + lam*ufl.tr(e)*ufl.Identity(3)
    measure = ufl.dx(metadata={"quadrature_degree": 4})
    a = ufl.inner(C(eps(trial)), eps(test))*r*measure
    L = ufl.inner(C(theta*ufl.Identity(3)), eps(test))*r*measure
    bcs = []
    for enabled, component, axis, value in (
        (bc_config.axis_ur_fixed, 0, 0, 0.),
        (bc_config.bottom_ur_fixed, 0, 1, 0.),
        (bc_config.bottom_uz_fixed, 1, 1, 0.),
        (bc_config.outer_ur_fixed, 0, 0, outer),
        (bc_config.outer_uz_fixed, 1, 0, outer)):
        if enabled:
            collapsed, _ = V.sub(component).collapse()
            zero = Function(collapsed)
            dofs = fem.locate_dofs_geometrical((V.sub(component), collapsed),
                lambda points, axis=axis, value=value: np.isclose(points[axis], value))
            bcs.append(fem.dirichletbc(zero, dofs, V.sub(component)))
    plastic_state = None
    plastic_history = []
    if plasticity_cu is not None or plasticity_dielectric is not None:
        if phi is not None:
            raise NotImplementedError("Plasticity currently requires a sharp fixed reference material map")
        if load_factors is None:
            if not isinstance(load_steps, int) or load_steps < 1:
                raise ValueError("load_steps must be a positive integer")
            load_factors = np.linspace(0., 1., load_steps+1)[1:]
        from plasticity import solve_plastic_path
        u, elastic, sigma, plastic_state, plastic_history = solve_plastic_path(
            V, eps, C, mu, theta, cu, dielectric_mask,
            material_cu, dielectric, plasticity_cu, plasticity_dielectric,
            bcs, load_factors, measure)
    else:
        if load_factors is not None:
            raise ValueError("load_factors is supported only for plastic load paths")
        u = _linear_solve(a, L, bcs)
        u.name = "displacement_dimless"
        elastic = eps(u) - theta*ufl.Identity(3)
        sigma = C(elastic)
    energy = 0.5*ufl.inner(sigma, elastic)

    def project(expr, name):
        q, v = ufl.TrialFunction(Q), ufl.TestFunction(Q)
        result = _linear_solve(q*v*r*measure, expr*v*r*measure, prefix=name+"_")
        result.name = name
        return result

    W_full = project(energy, "W_full_J_m3")
    W = project(cu*energy, "W_Cu_J_m3")
    # The actual SiO2 CSV W is FULL elastic energy, including all normal
    # components and sigma_rz*gamma_rz inside the half factor. The exported
    # shear strain is engineering gamma_rz = 2*epsilon_rz, as verified from G.
    W_legacy = project(cu*energy, "W_ansys_like_J_m3")
    stresses = {name: project(sigma[i,j], name) for name, i, j in
                (("sigma_rr",0,0), ("sigma_zz",1,1), ("sigma_tt",2,2), ("sigma_rz",0,1))}
    physical = Function(V, name="displacement_m")
    physical.x.array[:] = u.x.array*g.L0_m
    physical.x.scatter_forward()

    def surface(z_value, r_min, r_max):
        # Exact mesh surface rows, not an average through a finite z band.
        coords = V.tabulate_dof_coordinates()
        count = V.dofmap.index_map.size_local
        mask = (np.isclose(coords[:count,1], z_value, atol=1e-10, rtol=0)
                & (coords[:count,0] >= r_min-1e-10)
                & (coords[:count,0] <= r_max+1e-10))
        local = np.column_stack((coords[:count,:2][mask], u.x.array.reshape(-1,2)[:count][mask]))
        rows = np.concatenate(mesh.comm.allgather(local))
        if not len(rows):
            raise ValueError("Mesh must contain the requested reference surface row")
        rows = np.unique(rows, axis=0)
        rows = rows[np.argsort(rows[:,0])]
        return rows

    cu_surface = surface(g.cu_top if surface_z_dimless is None else surface_z_dimless, 0., radius)
    dielectric_surface = surface(g.top, radius, outer)
    integrate = lambda expr: mesh.comm.allreduce(fem.assemble_scalar(fem.form(expr*r*measure)), op=MPI.SUM)
    owned = Q.dofmap.index_map.size_local
    maximum = mesh.comm.allreduce(float(np.max(W.x.array[:owned], initial=0)), op=MPI.MAX)
    cu_volume = integrate(cu)
    stats = {"solve_time_s": time.perf_counter()-t_start,
             "W_max_J_m3": maximum,
             "W_mean_J_m3": float(integrate(cu*energy)/cu_volume) if cu_volume > 0 else 0.0,
             "UY_center_nm": float(cu_surface[0,3]*g.L0_nm),
             "UY_edge_nm": float(cu_surface[-1,3]*g.L0_nm),
             "dielectric_edge_UY_nm": float(dielectric_surface[0,3]*g.L0_nm),
             "relative_edge_height_nm": float((cu_surface[-1,1]+cu_surface[-1,3]
                -dielectric_surface[0,1]-dielectric_surface[0,3])*g.L0_nm),
             "residual_void_E_Pa": residual_void_E_Pa}
    stats["plasticity_enabled"] = plastic_state is not None
    if plastic_state is not None:
        stats["plastic_load_history"] = plastic_history
        stats["max_equivalent_plastic_strain"] = plastic_history[-1]["max_equivalent_plastic_strain"]
    p_field = project(plastic_state.equivalent_strain, "equivalent_plastic_strain") if plastic_state is not None else None
    return ThermoelasticResult(W, cu_surface[:,3], cu_surface[:,0], physical,
                              stresses, stats, W_full, W_legacy, u,
                              (cu_surface[:,:2]+cu_surface[:,2:])*g.L0_m,
                              (dielectric_surface[:,:2]+dielectric_surface[:,2:])*g.L0_m, plastic_state, p_field)


def solve(phi, mesh, material_cu, T_anneal_K, L0_m, pad_radius_dimless,
          pad_height_dimless, tracker=None, contact_state="pre", p_contact=None,
          surface_z_dimless=None, steady_state_thermal=True, *, geometry=None,
          material_dielectric=None, boundary_conditions=None, residual_void_E_Pa=1e3,
                     plasticity_cu=None, plasticity_dielectric=None, load_steps=20, load_factors=None):
    """Uniform-temperature mechanics; optional incremental J2 thermal path.

    Plastic paths start from virgin material at the reference temperatures on
    every call. load_factors scales the target thermal increment and can include
    unloading; it does not specify elapsed time or solve transient heat flow.
    """
    if not steady_state_thermal:
        raise NotImplementedError("Transient heating requires a validated thermal material map")
    tracker = tracker or MeshNodeTracker(mesh, L0_m)
    dielectric = material_dielectric or sio2_dielectric()
    temperature = solve_thermal(mesh, tracker, material_cu, dielectric, T_anneal_K)
    return solve_structural(mesh, tracker, material_cu, temperature,
        pad_radius_dimless, pad_height_dimless, L0_m, phi, contact_state,
        p_contact, surface_z_dimless, geometry=geometry, material_dielectric=dielectric,
        boundary_conditions=boundary_conditions, residual_void_E_Pa=residual_void_E_Pa,
        plasticity_cu=plasticity_cu, plasticity_dielectric=plasticity_dielectric,
        load_steps=load_steps, load_factors=load_factors)


if __name__ == "__main__":
    from validate_baseline import main
    main()
