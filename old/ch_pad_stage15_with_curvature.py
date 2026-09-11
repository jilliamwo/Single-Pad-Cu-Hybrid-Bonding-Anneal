"""
Stage 15: Cu surface diffusion initialized from ANSYS thermal
displacement, with optional ANSYS thermoelastic strain-energy W coupling.

This version keeps the physically calibrated Cu/TaN wall-energy term, but
uses the thermally deformed Cu top surface from ANSYS as the initial
condition for phi. This is the Stage 15 test of CTE/thermoelastic
geometry: the ANSYS vertical directional deformation at full anneal
temperature is sampled along the original Cu top surface and added to
the initial recessed surface before surface diffusion begins.

The goal is to separate instantaneous thermoelastic deformation from
slower surface-diffusion evolution.

Implements the modified Cahn-Hilliard structure of Hooge et al.,
"Modeling of Copper Hybrid Bonding Anneal", ECTC 2024, Eqs. 5-6,
using a numerically scaled geometry and chemical potential while
retaining time in seconds.

- uses a static wall free-energy difference
      Delta gamma_w = gamma_liner-void - gamma_Cu-liner = -0.651 J/m^2
- uses a tunable initial CMP profile with full Cu/TaN wall wetting
- adds a reduced mobility zone near the outer contact line/edge
- retains the latest C_D = 0.10 diffusivity-calibration trial as an explicit, independently adjustable kinetic parameter
- ANSYS strain-energy-density hook W(x,t)

Current scope:
- axisymmetric single recessed Cu pad
- curvature-driven surface diffusion
- static liner wall free energy
- imported nodal ANSYS W(r,z) from combined_W.csv
- W is interpolated from ANSYS coordinates onto the FEniCSx mesh

The apparent contact angle is treated as an output measured from the
phi = 0.5 contour. C_D remains a calibration parameter and should be
validated independently against measured bulge-out-versus-time data.
"""

from pathlib import Path
import csv


from mpi4py import MPI
import numpy as np
import pandas as pd
import ufl

from basix.ufl import element, mixed_element
from dolfinx import default_real_type, fem, geometry
from dolfinx.fem import Function, functionspace
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.io import XDMFFile
from dolfinx import mesh as dmesh
from dolfinx.mesh import CellType, create_rectangle

# 1. DIMENSIONAL MATERIAL AND KINETIC DATA
#
# Paper-aligned implementation of the modified Cahn-Hilliard
# equations (paper Eqs. 5-6). Geometry and chemical potential
# are numerically scaled, but time remains in physical seconds.
#
# IMPORTANT: h and the low-temperature extrapolation of D_s are
# model assumptions and should ultimately be calibrated against
# measured bulge-out-versus-time data.

K_BOLTZMANN = 1.380649e-23  # J/K
R_GAS = 8.314462618  # J/(mol K), retained for reporting Q

# Cu surface self-diffusion Arrhenius fit retained from Stage 7.
D0_SURFACE_M2_S = 0.26e-4  # 0.26 cm^2/s
Q_SURFACE_J_MOL = 20.7 * 4184.0  # 20.7 kcal/mol

# Physical Cu parameters.
GAMMA_CU_J_M2 = 1.70
OMEGA_CU_M3_ATOM = 1.18e-29
SURFACE_LAYER_THICKNESS_M = 0.30e-9

# Effective bulk surface-diffusion calibration. C_D = 1 reproduces the
# retained literature Arrhenius fit. Values below one slow the amount of
# redistribution reached during a finite anneal without changing the
# infinite-time equilibrium shape.
DIFFUSIVITY_CALIBRATION_FACTOR = 0.1

# Adaptive physical time stepping. A failed nonlinear step is rolled back
# and retried with a smaller dt. After several easy successful steps, dt is
# increased gradually up to MAX_DT_S.
INITIAL_DT_S = 1.0
MIN_DT_S = 0.03125
MAX_DT_S = 1.0
DT_REDUCTION_FACTOR = 0.5
DT_GROWTH_FACTOR = 1.25
EASY_NEWTON_ITERATIONS = 8
EASY_STEPS_BEFORE_GROWTH = 20
MAX_RETRIES_PER_STEP = 8


def surface_diffusivity(temperature_K: float) -> float:
    """Return Cu surface diffusivity D_s in m^2/s."""

    return D0_SURFACE_M2_S * np.exp(
        -Q_SURFACE_J_MOL / (R_GAS * temperature_K)
    )


def surface_atomic_mobility(temperature_K: float) -> float:
    """
    Return the explicit atomic surface mobility used in Eq. 5.

    The paper writes M_Gamma = D_Gamma h/(R T). Here D_s is a
    conventional diffusivity in m^2/s, so the surface atomic
    concentration is made explicit through 1/Omega_Cu and the
    atomic convention k_B T is used:

        M_Gamma = D_s h / (Omega_Cu k_B T)

    Units are 1/(J s). Multiplication by Omega_Cu in Eq. 5 gives
    the required interfacial transport coefficient.
    """

    D_s = (
        DIFFUSIVITY_CALIBRATION_FACTOR
        * surface_diffusivity(temperature_K)
    )
    return (
        D_s
        * SURFACE_LAYER_THICKNESS_M
        / (OMEGA_CU_M3_ATOM * K_BOLTZMANN * temperature_K)
    )


anneal_temperature_C = 250.0
anneal_temperature_K = anneal_temperature_C + 273.15

D_surface_raw = surface_diffusivity(anneal_temperature_K)
D_surface_effective = (
    DIFFUSIVITY_CALIBRATION_FACTOR * D_surface_raw
)
M_gamma = surface_atomic_mobility(anneal_temperature_K)

# 2. SCALED GEOMETRY

L0_nm = 500.0

L0_m = L0_nm * 1.0e-9

pad_radius = 250.0 / L0_nm
pad_height = 500.0 / L0_nm
recess = 3.5 / L0_nm
void_height = 50.0 / L0_nm

# Tunable post-CMP initial Cu surface profile.
# Key Stage 13C change: the contact angle is not imposed by
# wall energy. Instead, the initial Cu topography supplies curvature that
# can drive surface diffusion.
#
# Allowed modes:
#   "flat"                  -> flat recessed Cu surface
#   "center_dish"           -> center is lower than edge
#   "edge_rolloff"          -> edge/corner is lower than interior
#   "combined_dish_rolloff" -> both effects together
#
# Use positive amplitudes to lower the surface locally.
INITIAL_PROFILE_MODE = "flat"
CENTER_DISH_NM = 0.0
CENTER_PROTRUSION_NM = 0.0
EDGE_ROLLOFF_NM = 0.0
EDGE_ROLLOFF_START_RADIUS_NM = 225.0

# If True, subtract the axisymmetric area-average of the profile correction
# so the total initial Cu volume is comparable to the flat 3.5 nm-recess
# baseline. This isolates shape/curvature effects from simply changing the
# amount of copper.
INITIAL_PROFILE_VOLUME_NORMALIZE = False

# Optional uniform upward shift for later thermal-expansion tests. Keep at
# zero for the pure CMP-profile test.
INITIAL_THERMAL_EXPANSION_SHIFT_NM = 0.0

# Contact-line / edge mobility limiter.
# The surface diffusion mobility is multiplied by a factor smaller than one
# in the outer annular region near the Cu/liner edge. Physically this
# represents contact-line pinning, roughness, oxidation, liner/corner
# barriers, or slower transport near the Ta/TaN/dielectric boundary.
#
# A value of 1.0 recovers the unlimited edge-motion behavior. Values like
# 0.3, 0.1, 0.03 should preserve more memory of the initial CMP profile
# during finite anneal.
#
# Not using as of now and for future modeling.
EDGE_MOBILITY_LIMITER_ENABLED = False
EDGE_CONTACT_LINE_MOBILITY_FACTOR = 1.00
EDGE_MOBILITY_START_RADIUS_NM = 225.0


center_dish = CENTER_DISH_NM / L0_nm
center_protrusion = CENTER_PROTRUSION_NM / L0_nm
edge_rolloff = EDGE_ROLLOFF_NM / L0_nm
edge_rolloff_start_radius = EDGE_ROLLOFF_START_RADIUS_NM / L0_nm
initial_thermal_expansion_shift = INITIAL_THERMAL_EXPANSION_SHIFT_NM / L0_nm
edge_mobility_start_radius = EDGE_MOBILITY_START_RADIUS_NM / L0_nm

# ANSYS thermal displacement initial-condition import.
# Export ANSYS Mechanical > Directional Deformation in the vertical axis
# as a nodal text file. For the current 2D axisymmetric geometry, X Location
# is radial position and Y Location is vertical position.
#
# IMPORTANT: The files uploaded with stress/strain do not contain displacement.
# You need a separate ANSYS directional-deformation export, e.g.
#   Node Number, X Location (m), Y Location (m), Z Location (m), Directional Deformation (m)
# The code samples nodes near the original flat Cu top surface at
# y = ANSYS_CU_BOTTOM_Y_NM + (pad_height - recess).
ANSYS_THERMAL_DISPLACEMENT_INITIAL_PROFILE_ENABLED = True
# Rename this file for each ANSYS material/geometry case, e.g.
# "250C_500nm_3p5recess.csv" or "300C_500nm_3p5recess_SiO2.csv".
DIELECTRIC_MATERIAL_LABEL = "PolyimideA"
ANSYS_RESULTS_CSV = Path("ansys-results/250_500nm_3.5_PolyimideA.csv")
ANSYS_THERMAL_DISPLACEMENT_PATH = ANSYS_RESULTS_CSV
ANSYS_THERMAL_DISPLACEMENT_SURFACE_TOLERANCE_NM = 3.0
ANSYS_THERMAL_DISPLACEMENT_VALUE_SCALE = 1.0
ANSYS_CU_BOTTOM_Y_NM = 500.0
ANSYS_THERMAL_DISPLACEMENT_ZERO_REFERENCE = "none"
# allowed zero references:
#   "none"   -> use the absolute ANSYS vertical displacement
#   "center" -> subtract the r=0 displacement, preserving only shape variation
#   "edge"   -> subtract the r=R displacement, preserving only shape variation

domain_height = pad_height + void_height

interface_width_nm = 5.0
interface_width = interface_width_nm / L0_nm

# Static wall free-energy difference calibrated from the TEM-derived
# apparent angle. It is entered directly rather than prescribing a slope.
# Young relation for the adopted convention:
#   Delta gamma_w = gamma_liner-void - gamma_Cu-liner
#                 = gamma_Cu-surface * cos(theta_eq)
# Reference wall free-energy difference inferred from the paper's
# TEM-measured apparent Cu/Ta-TaN angle. Stage 15 keeps this wetting preference as a material/interface property and
# multiplies the active wall term by WALL_ENERGY_STRENGTH.
WALL_ENERGY_DIFFERENCE_J_M2 = -0.651

# Stage 15: keep the calibrated Cu/TaN wall-energy preference fixed.
# Initial shape now comes from ANSYS thermal displacement rather than
# manually prescribed dishing/edge-rolloff by default.
WALL_ENERGY_STRENGTH = 1.0

reference_wall_cosine = np.clip(
    WALL_ENERGY_DIFFERENCE_J_M2 / GAMMA_CU_J_M2,
    -1.0,
    1.0,
)
reference_wall_equilibrium_angle_rad = np.arccos(reference_wall_cosine)
reference_wall_equilibrium_angle_deg = np.degrees(reference_wall_equilibrium_angle_rad)

active_wall_energy_difference_J_m2 = (
    WALL_ENERGY_STRENGTH * WALL_ENERGY_DIFFERENCE_J_M2
)
active_wall_cosine = np.clip(
    active_wall_energy_difference_J_m2 / GAMMA_CU_J_M2,
    -1.0,
    1.0,
)
active_wall_equilibrium_angle_rad = np.arccos(active_wall_cosine)
active_wall_equilibrium_angle_deg = np.degrees(active_wall_equilibrium_angle_rad)

liner_void_minus_cu_liner_J_m2 = WALL_ENERGY_DIFFERENCE_J_M2

# Numerical noise below this value is not counted as bulge-out.
positive_tolerance_nm = 0.01
positive_tolerance = positive_tolerance_nm / L0_nm


# 3. AXISYMMETRIC MESH

MESH_CELLS_RADIAL = 100
MESH_CELLS_VERTICAL = 220

msh = create_rectangle(
    MPI.COMM_WORLD,
    [
        np.array([0.0, 0.0], dtype=default_real_type),
        np.array(
            [pad_radius, domain_height],
            dtype=default_real_type,
        ),
    ],
    [MESH_CELLS_RADIAL, MESH_CELLS_VERTICAL],
    CellType.triangle,
)


# 4. MARK THE Cu-Ta/TaN SIDEWALL

facet_dim = msh.topology.dim - 1
wall_marker = 1

right_wall_facets = dmesh.locate_entities_boundary(
    msh,
    facet_dim,
    lambda x: np.isclose(x[0], pad_radius),
)

if len(right_wall_facets) == 0:
    raise RuntimeError("No right-wall facets were found.")

facet_tags = dmesh.meshtags(
    msh,
    facet_dim,
    right_wall_facets,
    np.full(
        len(right_wall_facets),
        wall_marker,
        dtype=np.int32,
    ),
)


# 5. MIXED PHI-MU FUNCTION SPACE

P1 = element(
    "Lagrange",
    msh.basix_cell(),
    1,
    dtype=default_real_type,
)

ME = functionspace(
    msh,
    mixed_element([P1, P1]),
)

q, v = ufl.TestFunctions(ME)

u = Function(ME, name="current_solution")
u_old = Function(ME, name="previous_solution")

phi, mu = ufl.split(u)
phi_old, mu_old = ufl.split(u_old)


# 6. INITIAL COPPER SURFACE: TUNABLE POST-CMP PROFILE

def quintic_smoothstep(value):
    """C2-smooth transition from zero to one for scalar/array input."""

    clipped = np.clip(value, 0.0, 1.0)
    return (
        6.0 * clipped**5
        - 15.0 * clipped**4
        + 10.0 * clipped**3
    )


def _as_array(radial_position):
    return np.asarray(radial_position, dtype=float)


def raw_initial_profile_correction(radial_position):
    """
    Return the unnormalized initial surface correction in scaled units.

    Negative correction means the Cu surface is lower than the nominal
    flat recessed surface. The returned correction is added to
    pad_height - recess.
    """

    r = _as_array(radial_position)
    correction = np.zeros_like(r, dtype=float)
    normalized_radius = np.clip(r / pad_radius, 0.0, 1.0)

    if INITIAL_PROFILE_MODE in ("flat", "none"):
        pass

    elif INITIAL_PROFILE_MODE == "center_dish":
        # Smooth parabolic-like dish with zero slope at r=0 and r=R.
        dish_shape = (1.0 - normalized_radius**2) ** 2
        correction -= center_dish * dish_shape
    
    elif INITIAL_PROFILE_MODE == "center_protrusion":
        # Same shape as center_dish, but raised instead of recessed.
        protrusion_shape = (1.0 - normalized_radius**2) ** 2
        correction += center_protrusion * protrusion_shape
    
    elif INITIAL_PROFILE_MODE == "edge_rolloff":
        width = pad_radius - edge_rolloff_start_radius
        if width <= 0.0:
            raise ValueError("EDGE_ROLLOFF_START_RADIUS_NM must be smaller than pad radius.")
        edge_shape = quintic_smoothstep((r - edge_rolloff_start_radius) / width)
        correction -= edge_rolloff * edge_shape

    elif INITIAL_PROFILE_MODE == "combined_dish_rolloff":
        dish_shape = (1.0 - normalized_radius**2) ** 2
        correction -= center_dish * dish_shape

        width = pad_radius - edge_rolloff_start_radius
        if width <= 0.0:
            raise ValueError("EDGE_ROLLOFF_START_RADIUS_NM must be smaller than pad radius.")
        edge_shape = quintic_smoothstep((r - edge_rolloff_start_radius) / width)
        correction -= edge_rolloff * edge_shape
    
    elif INITIAL_PROFILE_MODE == "combined_protrusion_rolloff":
        protrusion_shape = (1.0 - normalized_radius**2) ** 2
        correction += center_protrusion * protrusion_shape

        width = pad_radius - edge_rolloff_start_radius
        if width <= 0.0:
            raise ValueError(
                "EDGE_ROLLOFF_START_RADIUS_NM must be smaller than pad radius."
            )
        
        edge_shape = quintic_smoothstep(
            (r - edge_rolloff_start_radius) / width
        )
        correction -= edge_rolloff * edge_shape

    else:
        raise ValueError(
            f"Unknown INITIAL_PROFILE_MODE={INITIAL_PROFILE_MODE!r}. Use "
            "'flat', 'center_dish', 'edge_rolloff', or 'combined_dish_rolloff'."
        )

    return correction


def axisymmetric_area_average(values, radii):
    """Return 2/R^2 * integral(values * r dr) over the pad radius."""

    return (
        2.0
        / pad_radius**2
        * np.trapezoid(values * radii, radii)
    )


_profile_average_r = np.linspace(0.0, pad_radius, 2001)
_profile_average_raw = raw_initial_profile_correction(_profile_average_r)
INITIAL_PROFILE_AREA_AVERAGE_CORRECTION = axisymmetric_area_average(
    _profile_average_raw,
    _profile_average_r,
)


def _sanitize_column_name(name: str) -> str:
    """Match np.genfromtxt-style sanitized column names robustly."""

    return (
        name.strip()
        .replace(" ", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("/", "_")
        .replace("-", "_")
    )


def _find_column(names, required_fragments, forbidden_fragments=()):
    lowered = {name: name.lower() for name in names}
    for name, low in lowered.items():
        if all(fragment.lower() in low for fragment in required_fragments) and not any(
            fragment.lower() in low for fragment in forbidden_fragments
        ):
            return name
    return None


def _read_ansys_table(path: Path) -> pd.DataFrame:
    """
    Read an ANSYS tabular export robustly.

    This accepts:
    - true comma-separated CSV files,
    - tab-separated ANSYS text exports, and
    - Excel workbooks even if they were accidentally saved with a .csv suffix.

    The latest combined result file is expected to include both
    Y Deformation (m) and W Strain Energy Density (J/m^3) columns.
    """

    if not path.exists():
        raise FileNotFoundError(
            f"ANSYS table was not found: {path.resolve()}\n"
            "Place the combined ANSYS result file at this path or update the "
            "ANSYS_THERMAL_DISPLACEMENT_PATH / ANSYS_W_CSV_PATH settings."
        )

    # Some ANSYS/Excel exports are actually .xlsx workbooks renamed as .csv.
    # The uploaded results1.csv has this structure, so try Excel first if it
    # is a ZIP-style workbook.
    try:
        import zipfile
        if zipfile.is_zipfile(path):
            return pd.read_excel(path, sheet_name=0)
    except Exception:
        pass

    # Try common CSV/text encodings and separators.
    errors = []
    for encoding in ("utf-8-sig", "utf-16", "cp1252", "latin1"):
        for sep in (None, ",", "\t"):
            try:
                return pd.read_csv(path, encoding=encoding, sep=sep, engine="python")
            except Exception as exc:
                errors.append(f"encoding={encoding}, sep={sep}: {type(exc).__name__}: {exc}")

    raise ValueError(
        f"Could not read ANSYS table from {path}. Tried Excel, CSV, and TSV. "
        f"Last errors: {errors[-3:]}"
    )


def _numeric_column(df: pd.DataFrame, col: str) -> np.ndarray:
    """Return a dataframe column as float, tolerating ANSYS text formatting."""
    return pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)


def _read_ansys_directional_displacement_surface(path: Path):
    """
    Read ANSYS vertical/y directional deformation from the combined result
    table and return sorted surface radii and vertical displacements in
    FEniCSx dimensionless coordinates.

    The source coordinates are dimensional meters. The selected source nodes
    are those near the original flat Cu top surface in the ANSYS model.
    """

    df = _read_ansys_table(path)
    names = list(df.columns)

    x_col = _find_column(names, ["X", "Location"])
    y_col = _find_column(names, ["Y", "Location"])

    # Prefer the explicitly signed vertical deformation column from the
    # combined file. Avoid total deformation, because it loses the sign.
    disp_col = _find_column(names, ["Y", "Deformation"], ["Location", "Total"])
    if disp_col is None:
        disp_col = _find_column(names, ["Directional", "Deformation"], ["Location", "Total"])
    if disp_col is None:
        disp_col = _find_column(names, ["Deformation"], ["Location", "Total"])
    if disp_col is None:
        disp_col = _find_column(names, ["UY"], ["Location", "Total"])
    if disp_col is None:
        disp_col = _find_column(names, ["Displacement", "Y"], ["Location", "Total"])

    missing = []
    if x_col is None:
        missing.append("X Location")
    if y_col is None:
        missing.append("Y Location")
    if disp_col is None:
        missing.append("Y/vertical directional deformation")
    if missing:
        raise ValueError(
            f"Missing required columns {missing} in {path}. Available columns: {names}"
        )

    x_m = _numeric_column(df, x_col)
    y_m = _numeric_column(df, y_col)
    uy_m = _numeric_column(df, disp_col) * ANSYS_THERMAL_DISPLACEMENT_VALUE_SCALE

    finite = np.isfinite(x_m) & np.isfinite(y_m) & np.isfinite(uy_m)
    x_m = x_m[finite]
    y_m = y_m[finite]
    uy_m = uy_m[finite]

    ansys_original_top_y_m = (
        ANSYS_CU_BOTTOM_Y_NM * 1.0e-9
        + (pad_height - recess) * L0_m
    )
    surface_tol_m = ANSYS_THERMAL_DISPLACEMENT_SURFACE_TOLERANCE_NM * 1.0e-9

    surface_mask = (
        (x_m >= -1.0e-15)
        & (x_m <= pad_radius * L0_m + 1.0e-15)
        & (np.abs(y_m - ansys_original_top_y_m) <= surface_tol_m)
    )

    if np.count_nonzero(surface_mask) < 3:
        raise RuntimeError(
            "Too few nodes were found near the original Cu top surface for "
            "thermal-displacement initialization. Check "
            "ANSYS_CU_BOTTOM_Y_NM, pad_height/recess, and "
            "ANSYS_THERMAL_DISPLACEMENT_SURFACE_TOLERANCE_NM. "
            f"Expected top y ≈ {ansys_original_top_y_m * 1e9:.3f} nm; "
            f"found {np.count_nonzero(surface_mask)} nodes within "
            f"±{ANSYS_THERMAL_DISPLACEMENT_SURFACE_TOLERANCE_NM:.3f} nm. "
            f"Detected displacement column: {disp_col}."
        )

    r_dimless = x_m[surface_mask] / L0_m
    uy_dimless = uy_m[surface_mask] / L0_m

    order = np.argsort(r_dimless)
    r_dimless = r_dimless[order]
    uy_dimless = uy_dimless[order]

    # Average duplicate radii if present.
    rounded = np.round(r_dimless, decimals=14)
    unique_r, inverse = np.unique(rounded, return_inverse=True)
    if len(unique_r) != len(r_dimless):
        sums = np.zeros(len(unique_r), dtype=float)
        counts = np.zeros(len(unique_r), dtype=float)
        np.add.at(sums, inverse, uy_dimless)
        np.add.at(counts, inverse, 1.0)
        r_dimless = unique_r
        uy_dimless = sums / counts

    if ANSYS_THERMAL_DISPLACEMENT_ZERO_REFERENCE == "center":
        uy_dimless = uy_dimless - np.interp(0.0, r_dimless, uy_dimless)
    elif ANSYS_THERMAL_DISPLACEMENT_ZERO_REFERENCE == "edge":
        uy_dimless = uy_dimless - np.interp(pad_radius, r_dimless, uy_dimless)
    elif ANSYS_THERMAL_DISPLACEMENT_ZERO_REFERENCE == "none":
        pass
    else:
        raise ValueError(
            "ANSYS_THERMAL_DISPLACEMENT_ZERO_REFERENCE must be 'none', "
            "'center', or 'edge'."
        )

    print(f"Read ANSYS thermal IC from {path}: {len(r_dimless)} top-surface radii using column '{disp_col}'")
    return r_dimless, uy_dimless


if ANSYS_THERMAL_DISPLACEMENT_INITIAL_PROFILE_ENABLED:
    (
        ANSYS_THERMAL_SURFACE_RADII,
        ANSYS_THERMAL_SURFACE_UY,
    ) = _read_ansys_directional_displacement_surface(
        ANSYS_THERMAL_DISPLACEMENT_PATH
    )
else:
    ANSYS_THERMAL_SURFACE_RADII = np.array([0.0, pad_radius], dtype=float)
    ANSYS_THERMAL_SURFACE_UY = np.array([0.0, 0.0], dtype=float)


def ansys_thermal_displacement_correction(radial_position):
    """Return imported ANSYS vertical thermal displacement in scaled units."""

    r = _as_array(radial_position)
    if not ANSYS_THERMAL_DISPLACEMENT_INITIAL_PROFILE_ENABLED:
        return np.zeros_like(r, dtype=float)

    return np.interp(
        np.clip(r, 0.0, pad_radius),
        ANSYS_THERMAL_SURFACE_RADII,
        ANSYS_THERMAL_SURFACE_UY,
        left=ANSYS_THERMAL_SURFACE_UY[0],
        right=ANSYS_THERMAL_SURFACE_UY[-1],
    )


def initial_profile_correction(radial_position):
    correction = raw_initial_profile_correction(radial_position)

    if INITIAL_PROFILE_VOLUME_NORMALIZE:
        correction = correction - INITIAL_PROFILE_AREA_AVERAGE_CORRECTION

    return correction


def initial_surface_height(radial_position):
    """Return the post-CMP Cu surface height in scaled coordinates."""

    nominal_height = pad_height - recess
    return (
        nominal_height
        + initial_thermal_expansion_shift
        + initial_profile_correction(radial_position)
        + ansys_thermal_displacement_correction(radial_position)
    )


def initial_phi(x):
    radial_position = x[0]
    vertical_position = x[1]
    surface_height = initial_surface_height(radial_position)

    return 0.5 * (
        1.0
        - np.tanh(
            (vertical_position - surface_height)
            / (
                np.sqrt(2.0)
                * interface_width
            )
        )
    )

u.sub(0).interpolate(initial_phi)

u.sub(1).interpolate(
    lambda x: np.zeros(
        x.shape[1],
        dtype=default_real_type,
    )
)

u.x.scatter_forward()

u_old.x.array[:] = u.x.array
u_old.x.scatter_forward()


# 7. PAPER EQUATIONS 5-6: PHYSICAL SCALES

phase_evolution_prefactor = 3.0 / (2.0 * np.sqrt(2.0))
surface_energy_prefactor = 3.0 / np.sqrt(2.0)
mechanical_energy_prefactor = 15.0 / 8.0

interface_width_m = interface_width_nm * 1.0e-9

# Scale the paper chemical potential by gamma*Omega/epsilon.
# The solved variable mu is therefore O(1), while time remains seconds.
chemical_potential_scale_J_atom = (
    GAMMA_CU_J_M2
    * OMEGA_CU_M3_ATOM
    / interface_width_m
)

# Coefficient multiplying the dimensionless divergence term in Eq. 5.
# Units: 1/s. This is what converts the transient solve to seconds.
kinetic_rate_per_s = (
    phase_evolution_prefactor
    * OMEGA_CU_M3_ATOM
    * M_gamma
    * chemical_potential_scale_J_atom
    / (interface_width_m * L0_m**2)
)

# Dimensionless coefficient in the normalized curvature term of Eq. 6.
epsilon_squared_over_L0_squared = (
    interface_width_m / L0_m
) ** 2

# 7A. IMPORT ANSYS STRAIN-ENERGY DENSITY FIELD W(r,z)
#
# Stage 12 reads a nodal ANSYS export produced from the stress/elastic-
# strain components. The expected CSV columns are:
#
#   X Location (m), Y Location (m), W_J_m3
#
# Coordinates in ANSYS are dimensional meters. The FEniCSx mesh coordinates
# are dimensionless and scaled by L0_m, so this importer divides the ANSYS
# coordinates by L0_m before interpolation.
#
# Negative nodal W values can appear from nodal stress/strain averaging even
# though stored elastic energy density is nonnegative. By default, they are
# clipped to zero before coupling into the phase-field chemical potential.

ANSYS_W_ENABLED = True
ANSYS_W_CSV_PATH = ANSYS_RESULTS_CSV
ANSYS_W_CLIP_NEGATIVE = True

# Important for the ANSYS -> FEniCSx one-way coupling:
# the ANSYS export may include both Cu and surrounding dielectric.
# The phase-field model only evolves the Cu pad, so dielectric W should
# not be used as a Cu diffusion driving force. This mask keeps only points
# inside the initial Cu pad geometry and zeroes W outside that region.
ANSYS_W_MASK_TO_INITIAL_CU_PAD = True
ANSYS_W_MASK_PADDING_NM = 0.0

# Coordinate alignment between ANSYS and FEniCSx:
# In the ANSYS model, dielectric occupies 0 < z < 500 nm and Cu occupies
# 500 < z < 996.5 nm. In the phase-field model, the Cu pad is placed at
# 0 < z < 496.5 nm. Therefore subtract 500 nm from the ANSYS z-coordinate
# before masking/interpolating W onto the FEniCSx mesh.
ANSYS_W_Z_ORIGIN_SHIFT_NM = 500.0

ANSYS_W_DISCARD_SOURCE_POINTS_OUTSIDE_MASK = True
ANSYS_W_ZERO_DOF_VALUES_OUTSIDE_MASK = True

# After source masking, use zero outside the remaining Cu interpolation
# domain. This avoids nearest-neighbor leakage from Cu or dielectric points
# into the void/dielectric region of the phase-field mesh.
ANSYS_W_OUTSIDE_DOMAIN_VALUE = 0.0
MECHANICAL_W_MULTIPLIER = 1.0


def _initial_cu_pad_mask_dimless(points_dimless: np.ndarray) -> np.ndarray:
    """
    Return True for coordinates that belong to the initial Cu pad.

    The mask is evaluated in the same scaled coordinates used by the
    FEniCSx mesh. It removes ANSYS points from the surrounding dielectric
    and void before W is interpolated. The top boundary follows the selected Stage 13C initial profile.
    """

    r = points_dimless[:, 0]
    z = points_dimless[:, 1]
    padding = ANSYS_W_MASK_PADDING_NM / L0_nm
    surface = initial_surface_height(r)

    return (
        (r >= -padding)
        & (r <= pad_radius + padding)
        & (z >= -padding)
        & (z <= surface + padding)
    )


def _read_ansys_w_csv(csv_path: Path):
    """
    Read ANSYS strain-energy density W from the combined result table.

    Supports the old column name W_J_m3 and the new combined-file column
    W Strain Energy Density (J/m^3). Coordinates are shifted so that the
    ANSYS Cu bottom maps to z=0 in the phase-field model.
    """

    df = _read_ansys_table(csv_path)
    names = list(df.columns)

    x_col = _find_column(names, ["X", "Location"])
    y_col = _find_column(names, ["Y", "Location"])

    W_col = None
    # Prefer explicit strain-energy-density columns.
    for required in (["W", "J"], ["Strain", "Energy", "Density"], ["Energy", "Density"]):
        W_col = _find_column(names, required)
        if W_col is not None:
            break

    # Backward compatible exact-style fallback after pandas import.
    if W_col is None and "W_J_m3" in names:
        W_col = "W_J_m3"

    missing = []
    if x_col is None:
        missing.append("X Location")
    if y_col is None:
        missing.append("Y Location")
    if W_col is None:
        missing.append("W strain energy density")
    if missing:
        raise ValueError(
            f"Missing required columns {missing} in {csv_path}. Available columns: {names}"
        )

    x_dimless = _numeric_column(df, x_col) / L0_m

    # Shift ANSYS z/y so that the Cu bottom in ANSYS maps to z=0 in FEniCSx.
    # ANSYS: dielectric 0-500 nm, Cu 500-996.5 nm.
    # FEniCSx: Cu 0-496.5 nm.
    y_dimless = (
        _numeric_column(df, y_col) - ANSYS_W_Z_ORIGIN_SHIFT_NM * 1.0e-9
    ) / L0_m
    W_values = _numeric_column(df, W_col)

    finite = (
        np.isfinite(x_dimless)
        & np.isfinite(y_dimless)
        & np.isfinite(W_values)
    )
    x_dimless = x_dimless[finite]
    y_dimless = y_dimless[finite]
    W_values = W_values[finite]

    if ANSYS_W_CLIP_NEGATIVE:
        W_values = np.maximum(W_values, 0.0)

    points = np.column_stack([x_dimless, y_dimless])

    if ANSYS_W_MASK_TO_INITIAL_CU_PAD and ANSYS_W_DISCARD_SOURCE_POINTS_OUTSIDE_MASK:
        source_mask = _initial_cu_pad_mask_dimless(points)
        points = points[source_mask]
        W_values = W_values[source_mask]

        if len(points) == 0:
            raise RuntimeError(
                "The ANSYS W mask removed every source point. Check that "
                "the ANSYS coordinate origin/units match the FEniCSx model, "
                "or temporarily set ANSYS_W_MASK_TO_INITIAL_CU_PAD = False."
            )

    # Average exact duplicate coordinate rows if present. This avoids
    # interpolation issues at coincident ANSYS nodes.
    unique_points, inverse = np.unique(points, axis=0, return_inverse=True)
    if len(unique_points) != len(points):
        sums = np.zeros(len(unique_points), dtype=float)
        counts = np.zeros(len(unique_points), dtype=float)
        np.add.at(sums, inverse, W_values)
        np.add.at(counts, inverse, 1.0)
        W_values = sums / counts
        points = unique_points

    print(f"Read ANSYS W from {csv_path}: {len(points)} masked source points using column '{W_col}'")
    return points, W_values


def _interpolate_ansys_w_to_function(csv_path: Path, V):
    points, W_values = _read_ansys_w_csv(csv_path)

    dof_coordinates = V.tabulate_dof_coordinates()[:, :2]

    try:
        from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

        linear_interp = LinearNDInterpolator(points, W_values, fill_value=np.nan)
        interpolated_values = linear_interp(
            dof_coordinates[:, 0],
            dof_coordinates[:, 1],
        )

        if np.any(~np.isfinite(interpolated_values)):
            if ANSYS_W_OUTSIDE_DOMAIN_VALUE == "nearest":
                nearest_interp = NearestNDInterpolator(points, W_values)
                missing = ~np.isfinite(interpolated_values)
                interpolated_values[missing] = nearest_interp(
                    dof_coordinates[missing, 0],
                    dof_coordinates[missing, 1],
                )
            else:
                interpolated_values[~np.isfinite(interpolated_values)] = float(
                    ANSYS_W_OUTSIDE_DOMAIN_VALUE
                )

    except Exception:
        # Fallback when SciPy is unavailable: nearest-neighbor mapping.
        # This is less smooth than linear interpolation but keeps the script
        # runnable and preserves the spatial W scale.
        distances_squared = (
            (dof_coordinates[:, None, 0] - points[None, :, 0]) ** 2
            + (dof_coordinates[:, None, 1] - points[None, :, 1]) ** 2
        )
        nearest = np.argmin(distances_squared, axis=1)
        interpolated_values = W_values[nearest]

    if ANSYS_W_MASK_TO_INITIAL_CU_PAD and ANSYS_W_ZERO_DOF_VALUES_OUTSIDE_MASK:
        dof_mask = _initial_cu_pad_mask_dimless(dof_coordinates)
        interpolated_values[~dof_mask] = 0.0

    return interpolated_values, W_values


V_W = functionspace(msh, P1)
W_field = Function(V_W, name="ansys_strain_energy_density_J_m3")

if ANSYS_W_ENABLED:
    W_dof_values, W_source_values = _interpolate_ansys_w_to_function(
        ANSYS_W_CSV_PATH,
        V_W,
    )
    W_field.x.array[:] = MECHANICAL_W_MULTIPLIER * W_dof_values
else:
    W_source_values = np.array([0.0], dtype=float)
    W_field.x.array[:] = 0.0

W_field.x.scatter_forward()

# Dimensionless mechanical energy appearing in normalized Eq. 6.
W_dimensionless = (
    interface_width_m / GAMMA_CU_J_M2
) * W_field

# Static wall free energy has no independent wall-relaxation time.


# 8. DIFFUSE INTERFACE DELTA FUNCTION
#
# Standard normalized localization profile. It vanishes in bulk Cu
# and void and peaks at one at phi = 0.5. The explicit 1/epsilon
# factor from paper Eq. 5 is already included in kinetic_rate_per_s.

delta_phi_old = 16.0 * phi_old**2 * (1.0 - phi_old) ** 2
delta_phi = 16.0 * phi**2 * (1.0 - phi) ** 2


# 9. DOUBLE-WELL ENERGY

phi_variable = ufl.variable(phi)

free_energy = (
    phi_variable**2
    * (1.0 - phi_variable) ** 2
)

df_dphi = ufl.diff(
    free_energy,
    phi_variable,
)


# 10. PHYSICAL TIME INTEGRATION

theta = 1.0

# Seconds. Reduce dt if Newton convergence or time-step convergence
# is poor. The paper plots evolution to approximately 7000 s.
dt = INITIAL_DT_S
final_time = 7000.0

# UFL uses this Constant so the time step can be changed without rebuilding
# the nonlinear form.
dt_constant = fem.Constant(msh, default_real_type(dt))

mu_mid = (1.0 - theta) * mu_old + theta * mu


# 11. AXISYMMETRIC MEASURES

radial_coordinate = ufl.SpatialCoordinate(msh)[0]

dx = ufl.Measure(
    "dx",
    domain=msh,
)

ds = ufl.Measure(
    "ds",
    domain=msh,
    subdomain_data=facet_tags,
)


# 11A. CONTACT-LINE / EDGE MOBILITY LIMITER
# The thermodynamic wall energy remains active, but atomic transport is
# slowed near the outer Cu/liner edge. Since the Cahn-Hilliard flux is
# already localized by delta_phi_old, this mainly limits motion of the
# moving Cu interface in the edge/contact-line region.
#
# This is intentionally a kinetic limiter, not a new equilibrium angle.
# It should change how much of the 112.5 degree equilibrium wetting state
# is reached during finite anneal time.

if EDGE_MOBILITY_LIMITER_ENABLED:
    contact_line_mobility_factor_ufl = ufl.conditional(
        ufl.gt(radial_coordinate, edge_mobility_start_radius),
        default_real_type(EDGE_CONTACT_LINE_MOBILITY_FACTOR),
        default_real_type(1.0),
    )
else:
    contact_line_mobility_factor_ufl = default_real_type(1.0)


# 12. WEAK FORM

F_phase = (
    (phi - phi_old)
    * q
    * radial_coordinate
    * dx

    + dt_constant
    * kinetic_rate_per_s
    * contact_line_mobility_factor_ufl
    * delta_phi_old
    * ufl.inner(
        ufl.grad(mu_mid),
        ufl.grad(q),
    )
    * radial_coordinate
    * dx
)

# Cubic wall interpolation p(phi), with p(0)=0 and p(1)=1.
# The wall free energy distinguishes liner exposed to void from liner
# covered by Cu. Its variation produces a natural wetting boundary
# condition; no interface slope is prescribed directly.
wall_interpolation = phi**2 * (3.0 - 2.0 * phi)
dwall_interpolation_dphi = 6.0 * phi * (1.0 - phi)

# Dimensionless wall-energy density consistent with the normalized
# bulk free energy used below. In Stage 13C the active wall term is full
# strength, while the edge mobility limiter controls how quickly the
# interface can approach this equilibrium wetting state.
wall_energy_density = (
    -0.5
    * interface_width
    * active_wall_cosine
    * wall_interpolation
)

dwall_energy_dphi = (
    -0.5
    * interface_width
    * active_wall_cosine
    * dwall_interpolation_dphi
)

# Normalized form of paper Eq. 6. The physical chemical
# potential is mu * chemical_potential_scale_J_atom.
F_mu = (
    mu
    * v
    * radial_coordinate
    * dx

    - surface_energy_prefactor
    * epsilon_squared_over_L0_squared
    * ufl.inner(
        ufl.grad(phi),
        ufl.grad(v),
    )
    * radial_coordinate
    * dx

    - surface_energy_prefactor
    * df_dphi
    * v
    * radial_coordinate
    * dx

    + mechanical_energy_prefactor
    * delta_phi
    * W_dimensionless
    * v
    * radial_coordinate
    * dx

    # Natural boundary contribution from the liner wall free energy.
    # This replaces the former fixed normal-gradient/contact-angle term.
    - dwall_energy_dphi
    * v
    * radial_coordinate
    * ds(wall_marker)

)

F = F_phase + F_mu


# 13. NONLINEAR SOLVER

SNES_RTOL = 1.0e-8
SNES_ATOL = 1.0e-10
SNES_STOL = 1.0e-10
SNES_MAX_IT = 120

petsc_options = {
    "snes_type": "newtonls",
    "snes_linesearch_type": "bt",
    "snes_rtol": SNES_RTOL,
    "snes_atol": SNES_ATOL,
    "snes_stol": SNES_STOL,
    "snes_max_it": SNES_MAX_IT,
    "ksp_type": "preonly",
    "pc_type": "lu",
}

problem = NonlinearProblem(
    F,
    u,
    petsc_options_prefix="paper_stage15_ansys_thermal_ic_ansys_W_ch_",
    petsc_options=petsc_options,
)


# 14. OUTPUT

def format_label_number(value):
    return f"{value:.2f}".replace("-", "m").replace(".", "p")


# Stage 15 clean output naming. Only include the primary parameters being swept.
# Example folder:
#   results_stage15_T250C_D500nm_R3p50nm_Wall1p00
wall_strength_label = format_label_number(WALL_ENERGY_STRENGTH)
recess_label = format_label_number(recess * L0_nm)
diameter_label = format_label_number(2.0 * pad_radius * L0_nm)
material_label = str(DIELECTRIC_MATERIAL_LABEL).replace(" ", "").replace("/", "-").replace("_", "-")
cd_label = format_label_number(DIFFUSIVITY_CALIBRATION_FACTOR)
profile_mode_label = str(INITIAL_PROFILE_MODE).replace("_", "-")
dish_label = format_label_number(CENTER_DISH_NM)
protrusion_label = format_label_number(CENTER_PROTRUSION_NM)
ero_label = format_label_number(EDGE_ROLLOFF_NM)

run_label = (
    f"T{int(anneal_temperature_C)}C"
    f"_D{diameter_label}nm"
    f"_R{recess_label}nm"
    f"_Wall{wall_strength_label}"
    f"_Mat{material_label}"
)

output_folder = Path(f"results_stage15_{run_label}")

output_folder.mkdir(
    parents=True,
    exist_ok=True,
)

V_phi, phi_dofs = ME.sub(0).collapse()

phi_output = Function(
    V_phi,
    name="phi",
)

xdmf_path = (
    output_folder
    / f"stage15_{run_label}_phi.xdmf"
)

xdmf = XDMFFile(
    MPI.COMM_WORLD,
    xdmf_path,
    "w",
)

xdmf.write_mesh(msh)

# Save the imported/interpolated ANSYS W field for ParaView inspection.
if ANSYS_W_ENABLED:
    W_xdmf_path = output_folder / f"stage15_{run_label}_interpolated_Cu_masked_ansys_W.xdmf"
    with XDMFFile(MPI.COMM_WORLD, W_xdmf_path, "w") as W_xdmf:
        W_xdmf.write_mesh(msh)
        W_xdmf.write_function(W_field, 0.0)
else:
    W_xdmf_path = None

# Save the actual selected initial CMP profile used to initialize phi.
if MPI.COMM_WORLD.rank == 0:
    initial_profile_path = output_folder / f"stage15_{run_label}_initial_profile.csv"
    profile_r_nm = np.linspace(0.0, pad_radius * L0_nm, 251)
    profile_z_nm = initial_surface_height(profile_r_nm / L0_nm) * L0_nm

    with open(initial_profile_path, "w", newline="") as profile_file:
        profile_writer = csv.writer(profile_file)
        profile_writer.writerow(["radius_nm", "surface_height_nm"])
        profile_writer.writerows(zip(profile_r_nm, profile_z_nm))


def update_output():
    phi_output.x.array[:] = u.x.array[phi_dofs]
    phi_output.x.scatter_forward()


def copper_mass():
    mass_form = fem.form(
        phi
        * radial_coordinate
        * dx
    )

    local_mass = fem.assemble_scalar(
        mass_form
    )

    return msh.comm.allreduce(
        local_mass,
        op=MPI.SUM,
    )


def phase_field_bounds():
    """Return global minimum and maximum of the current phase field."""
    update_output()
    local_min = float(np.min(phi_output.x.array))
    local_max = float(np.max(phi_output.x.array))
    global_min = msh.comm.allreduce(local_min, op=MPI.MIN)
    global_max = msh.comm.allreduce(local_max, op=MPI.MAX)
    return global_min, global_max


# 15. BULGE-OUT SAMPLING GRID

if MPI.COMM_WORLD.size != 1:
    raise RuntimeError(
        "Run this measurement script with one MPI rank."
    )

reference_height = pad_height - recess

number_of_radial_points = 201
number_of_vertical_points = 401

sampling_offset = 1.0e-10

radial_samples = np.linspace(
    sampling_offset,
    pad_radius - sampling_offset,
    number_of_radial_points,
)

vertical_samples = np.linspace(
    pad_height - 0.15,
    domain_height - sampling_offset,
    number_of_vertical_points,
)

sample_points = np.array(
    [
        [r_value, z_value, 0.0]
        for r_value in radial_samples
        for z_value in vertical_samples
    ],
    dtype=np.float64,
)

tree = geometry.bb_tree(
    msh,
    msh.topology.dim,
)

candidate_cells = geometry.compute_collisions_points(
    tree,
    sample_points,
)

colliding_cells = geometry.compute_colliding_cells(
    msh,
    candidate_cells,
    sample_points,
)

sample_cells = np.full(
    len(sample_points),
    -1,
    dtype=np.int32,
)

for point_index in range(len(sample_points)):
    cells = colliding_cells.links(point_index)

    if len(cells) > 0:
        sample_cells[point_index] = cells[0]

if np.any(sample_cells < 0):
    raise RuntimeError(
        "Some measurement points were not found in the mesh."
    )


# 16. BULGE-OUT MEASUREMENTS

def extract_interface_profile():
    """Return radius and height arrays for the upper phi=0.5 contour."""
    phi_output.x.scatter_forward()

    phase_values = phi_output.eval(
        sample_points,
        sample_cells,
    ).reshape(
        number_of_radial_points,
        number_of_vertical_points,
    )

    interface_heights = np.full(number_of_radial_points, np.nan)

    for radial_index in range(number_of_radial_points):
        phase_line = phase_values[radial_index]
        crossings = np.where(
            (phase_line[:-1] >= 0.5)
            & (phase_line[1:] < 0.5)
        )[0]

        if len(crossings) == 0:
            continue

        lower_index = crossings[-1]
        upper_index = lower_index + 1
        phi_lower = phase_line[lower_index]
        phi_upper = phase_line[upper_index]
        denominator = phi_upper - phi_lower

        if abs(denominator) < 1.0e-14:
            continue

        fraction = (0.5 - phi_lower) / denominator
        interface_heights[radial_index] = (
            vertical_samples[lower_index]
            + fraction
            * (vertical_samples[upper_index] - vertical_samples[lower_index])
        )

    valid = np.isfinite(interface_heights)
    if np.count_nonzero(valid) < 2:
        raise RuntimeError("Could not extract the phi = 0.5 interface.")

    return radial_samples[valid], interface_heights[valid]


def measure_bulge_out():
    r_valid, z_valid = extract_interface_profile()

    displacement = z_valid - reference_height
    initial_profile_valid = initial_surface_height(r_valid)
    evolution_from_initial = z_valid - initial_profile_valid

    # Estimate the apparent contact angle through the Cu phase by fitting
    # the final 10 nm of the phi=0.5 contour. A downward slope toward the
    # liner gives an angle greater than 90 degrees.
    angle_fit_width_nm = 10.0
    angle_fit_width = angle_fit_width_nm / L0_nm
    edge_fit_mask = r_valid >= (pad_radius - angle_fit_width)

    apparent_contact_angle_deg = np.nan
    edge_slope = np.nan

    if np.count_nonzero(edge_fit_mask) >= 3:
        edge_slope, _ = np.polyfit(
            r_valid[edge_fit_mask],
            z_valid[edge_fit_mask],
            1,
        )
        apparent_contact_angle_deg = (
            90.0 - np.degrees(np.arctan(edge_slope))
        )

    # Apply 0.01 nm classification tolerance.
    positive_boolean = (
        displacement > positive_tolerance
    )

    positive_mask = positive_boolean.astype(float)

    positive_displacement = np.where(
        positive_boolean,
        displacement,
        0.0,
    )

    maximum_bulge = np.max(displacement)

    # Positive-only area-equivalent bulge. Negative edge recession is
    # clipped to zero, matching the previous Stage 11 diagnostic.
    full_pad_average = (
        2.0
        / pad_radius**2
        * np.trapezoid(
            positive_displacement * r_valid,
            r_valid,
        )
    )

    # True signed area-average over the complete circular pad. Positive
    # center bulge contributes positively and edge recession contributes
    # negatively. The radial factor applies the correct axisymmetric area
    # weighting.
    signed_full_pad_average = (
        2.0
        / pad_radius**2
        * np.trapezoid(
            displacement * r_valid,
            r_valid,
        )
    )

    # Unweighted arithmetic average of sampled surface points. This is
    # included only as a secondary diagnostic because uniformly spaced
    # radius samples do not represent equal annular areas.
    signed_point_average = np.mean(displacement)

    positive_area_fraction = (
        2.0
        / pad_radius**2
        * np.trapezoid(
            positive_mask * r_valid,
            r_valid,
        )
    )

    if positive_area_fraction > 1.0e-14:
        positive_region_average = (
            full_pad_average
            / positive_area_fraction
        )
    else:
        positive_region_average = 0.0

    center_displacement = displacement[0]
    edge_displacement = displacement[-1]

    return {
        "maximum_bulge_nm":
            maximum_bulge * L0_nm,

        "full_pad_average_nm":
            full_pad_average * L0_nm,

        "signed_full_pad_average_nm":
            signed_full_pad_average * L0_nm,

        "signed_point_average_nm":
            signed_point_average * L0_nm,

        "positive_region_average_nm":
            positive_region_average * L0_nm,

        "center_displacement_nm":
            center_displacement * L0_nm,

        "edge_displacement_nm":
            edge_displacement * L0_nm,

        "maximum_change_from_initial_nm":
            np.max(evolution_from_initial) * L0_nm,

        "center_change_from_initial_nm":
            evolution_from_initial[0] * L0_nm,

        "edge_change_from_initial_nm":
            evolution_from_initial[-1] * L0_nm,

        "apparent_contact_angle_deg":
            apparent_contact_angle_deg,

        "edge_slope_dz_dr":
            edge_slope,

        "positive_area_fraction":
            positive_area_fraction,
    }


# 17. INITIAL OUTPUT AND CSV

update_output()

xdmf.write_function(
    phi_output,
    0.0,
)

initial_mass = copper_mass()

csv_path = output_folder / f"stage15_{run_label}_bulge_metrics.csv"

csv_file = open(
    csv_path,
    "w",
    newline="",
)

csv_writer = csv.writer(csv_file)

csv_writer.writerow(
    [
        "step",
        "time_s",
        "maximum_bulge_nm",
        "full_pad_average_positive_bulge_nm",
        "signed_full_pad_area_average_nm",
        "signed_unweighted_point_average_nm",
        "positive_region_average_bulge_nm",
        "center_displacement_nm",
        "edge_displacement_from_flat_reference_nm",
        "maximum_change_from_initial_profile_nm",
        "center_change_from_initial_profile_nm",
        "edge_change_from_initial_profile_nm",
        "apparent_contact_angle_deg",
        "edge_slope_dz_dr",
        "positive_area_fraction",
        "relative_mass_change",
        "accepted_dt_s",
        "newton_iterations",
        "snes_convergence_reason",
        "final_residual_norm",
        "phi_min",
        "phi_max",
    ]
)


def write_measurement(
    step_number,
    time_value,
    mass_change,
    accepted_dt=np.nan,
    newton_iterations=0,
    convergence_reason=0,
    residual_norm=np.nan,
):
    update_output()

    metrics = measure_bulge_out()
    phi_min, phi_max = phase_field_bounds()

    csv_writer.writerow(
        [
            step_number,
            time_value,
            metrics["maximum_bulge_nm"],
            metrics["full_pad_average_nm"],
            metrics["signed_full_pad_average_nm"],
            metrics["signed_point_average_nm"],
            metrics["positive_region_average_nm"],
            metrics["center_displacement_nm"],
            metrics["edge_displacement_nm"],
            metrics["maximum_change_from_initial_nm"],
            metrics["center_change_from_initial_nm"],
            metrics["edge_change_from_initial_nm"],
            metrics["apparent_contact_angle_deg"],
            metrics["edge_slope_dz_dr"],
            metrics["positive_area_fraction"],
            mass_change,
            accepted_dt,
            newton_iterations,
            convergence_reason,
            residual_norm,
            phi_min,
            phi_max,
        ]
    )

    csv_file.flush()

    return metrics


last_metrics = write_measurement(
    0,
    0.0,
    0.0,
)


# 18. RUN INFORMATION

if MPI.COMM_WORLD.rank == 0:
    print("\nStage 13C initial-profile-natural-angle + ANSYS-W surface-diffusion model")
    print("---------------------------------------------")
    print(f"Temperature: {anneal_temperature_C:.1f} C")
    print(f"Pad diameter: {2.0 * pad_radius * L0_nm:.1f} nm")
    print(f"Pad height: {pad_height * L0_nm:.1f} nm")
    print(f"Recess: {recess * L0_nm:.2f} nm")
    print(
        f"Reference wall-energy angle from paper: "
        f"{reference_wall_equilibrium_angle_deg:.2f} degrees"
    )
    print(
        f"Active wall energy strength: {WALL_ENERGY_STRENGTH:.3g}; "
        f"active wall term equivalent angle: {active_wall_equilibrium_angle_deg:.2f} degrees"
    )
    print(
        "Reference liner-void minus Cu-liner energy: "
        f"{liner_void_minus_cu_liner_J_m2:.6f} J/m^2"
    )
    print(
        "Active liner-void minus Cu-liner energy: "
        f"{active_wall_energy_difference_J_m2:.6f} J/m^2"
    )
    print("Finite wall mobility: disabled; Stage 15 uses ANSYS thermal-displacement initial condition")
    print(f"Edge mobility limiter enabled: {EDGE_MOBILITY_LIMITER_ENABLED}")
    print(f"Edge contact-line mobility factor: {EDGE_CONTACT_LINE_MOBILITY_FACTOR:.6g}")
    print(f"Edge mobility start radius: {EDGE_MOBILITY_START_RADIUS_NM:.3f} nm")
    print(f"Initial profile mode: {INITIAL_PROFILE_MODE}")
    print(f"Center dish amplitude: {CENTER_DISH_NM:.3f} nm")
    print(f"Edge roll-off amplitude: {EDGE_ROLLOFF_NM:.3f} nm")
    print(f"Edge roll-off start radius: {EDGE_ROLLOFF_START_RADIUS_NM:.3f} nm")
    print(f"Initial profile volume-normalized: {INITIAL_PROFILE_VOLUME_NORMALIZE}")
    print(
        f"Initial profile raw area-average correction: "
        f"{INITIAL_PROFILE_AREA_AVERAGE_CORRECTION * L0_nm:.6f} nm"
    )
    print(f"Initial thermal expansion shift: {INITIAL_THERMAL_EXPANSION_SHIFT_NM:.3f} nm")
    print(f"ANSYS thermal displacement IC enabled: {ANSYS_THERMAL_DISPLACEMENT_INITIAL_PROFILE_ENABLED}")
    if ANSYS_THERMAL_DISPLACEMENT_INITIAL_PROFILE_ENABLED:
        print(f"ANSYS thermal displacement path: {ANSYS_THERMAL_DISPLACEMENT_PATH}")
        print(f"ANSYS thermal top-surface tolerance: {ANSYS_THERMAL_DISPLACEMENT_SURFACE_TOLERANCE_NM:.3f} nm")
        print(f"ANSYS thermal zero reference: {ANSYS_THERMAL_DISPLACEMENT_ZERO_REFERENCE}")
        print(
            f"ANSYS thermal displacement min/mean/max: "
            f"{np.min(ANSYS_THERMAL_SURFACE_UY) * L0_nm:.6f} / "
            f"{np.mean(ANSYS_THERMAL_SURFACE_UY) * L0_nm:.6f} / "
            f"{np.max(ANSYS_THERMAL_SURFACE_UY) * L0_nm:.6f} nm"
        )
    print(f"Raw surface diffusivity D_s: {D_surface_raw:.6e} m^2/s")
    print(f"Diffusivity calibration C_D: {DIFFUSIVITY_CALIBRATION_FACTOR:.6g}")
    print(f"Effective surface diffusivity: {D_surface_effective:.6e} m^2/s")
    print(f"Effective surface atomic mobility M_Gamma: {M_gamma:.6e} 1/(J s)")
    print(f"Cu surface energy gamma: {GAMMA_CU_J_M2:.4f} J/m^2")
    print(f"Cu atomic volume Omega: {OMEGA_CU_M3_ATOM:.6e} m^3/atom")
    print(f"Surface layer thickness h: {SURFACE_LAYER_THICKNESS_M * 1e9:.3f} nm")
    print(f"Chemical-potential scale: {chemical_potential_scale_J_atom:.6e} J/atom")
    print(f"Kinetic rate coefficient: {kinetic_rate_per_s:.6e} 1/s")
    print(f"Positive tolerance: {positive_tolerance_nm:.3f} nm")
    print(f"Initial time step: {dt:.6e} s")
    print(f"Minimum time step: {MIN_DT_S:.6e} s")
    print(f"Maximum time step: {MAX_DT_S:.6e} s")
    print(f"Final time: {final_time:.6e} s")
    print(f"Mesh cells (r x z): {MESH_CELLS_RADIAL} x {MESH_CELLS_VERTICAL}")
    print(f"Interface width epsilon: {interface_width_nm:.3f} nm")
    print(f"SNES tolerances: rtol={SNES_RTOL:.1e}, atol={SNES_ATOL:.1e}, stol={SNES_STOL:.1e}")
    
    if ANSYS_W_ENABLED:
        W_local_min = float(np.min(W_field.x.array))
        W_local_max = float(np.max(W_field.x.array))
        W_local_mean = float(np.mean(W_field.x.array))
        W_min = msh.comm.allreduce(W_local_min, op=MPI.MIN)
        W_max = msh.comm.allreduce(W_local_max, op=MPI.MAX)
        W_mean = msh.comm.allreduce(W_local_mean, op=MPI.SUM) / msh.comm.size
        print(f"ANSYS W CSV: {ANSYS_W_CSV_PATH}")
        print(f"ANSYS W negative clipping: {ANSYS_W_CLIP_NEGATIVE}")
        print(f"ANSYS W mask to initial Cu pad: {ANSYS_W_MASK_TO_INITIAL_CU_PAD}")
        print(f"ANSYS W mask padding: {ANSYS_W_MASK_PADDING_NM:.6g} nm")
        print(f"ANSYS W outside-domain value: {ANSYS_W_OUTSIDE_DOMAIN_VALUE}")
        print(f"Mechanical W multiplier: {MECHANICAL_W_MULTIPLIER:.6g}")
        print(f"Interpolated W min/mean/max: {W_min:.6e} / {W_mean:.6e} / {W_max:.6e} J/m^3")
        print(f"W XDMF: {W_xdmf_path}")
    else:
        print("Mechanical strain energy W: disabled, zero field")
    print(f"CSV: {csv_path}")
    print(f"XDMF: {xdmf_path}\n")


# 19. TIME STEPPING

time = 0.0
step = 0

measure_interval_s = 50.0
save_interval_s = 100.0
next_measure_time = measure_interval_s
next_save_time = save_interval_s

last_measured_step = 0
last_saved_step = 0
easy_step_count = 0

while time < final_time - 1.0e-12:
    # Do not overshoot the requested final physical time.
    dt = min(dt, final_time - time)
    retry = 0

    while True:
        dt_constant.value = default_real_type(dt)

        # Always restart a retry from the last converged state. A failed
        # SNES solve may have overwritten u with a nonphysical iterate.
        u.x.array[:] = u_old.x.array
        u.x.scatter_forward()

        problem.solve()

        convergence_reason = problem.solver.getConvergedReason()
        nonlinear_iterations = problem.solver.getIterationNumber()
        residual_norm = problem.solver.getFunctionNorm()

        if convergence_reason > 0:
            break

        retry += 1
        proposed_dt = dt * DT_REDUCTION_FACTOR

        if MPI.COMM_WORLD.rank == 0:
            print(
                f"SNES retry at t={time:.6f} s: reason={convergence_reason}, "
                f"Newton={nonlinear_iterations}, dt {dt:.6g} -> {proposed_dt:.6g} s"
            )

        if retry > MAX_RETRIES_PER_STEP or proposed_dt < MIN_DT_S - 1.0e-15:
            xdmf.close()
            csv_file.close()
            raise RuntimeError(
                f"Solver failed after {retry} retries at accepted time "
                f"{time:.6e} s. Last reason {convergence_reason}; "
                f"minimum attempted dt={proposed_dt:.6e} s."
            )

        dt = max(proposed_dt, MIN_DT_S)

    # Accept the converged trial step.
    accepted_dt = dt
    time += accepted_dt
    step += 1

    u_old.x.array[:] = u.x.array
    u_old.x.scatter_forward()

    # Gradually recover the nominal step after a sequence of easy solves.
    if nonlinear_iterations <= EASY_NEWTON_ITERATIONS:
        easy_step_count += 1
    else:
        easy_step_count = 0

    if easy_step_count >= EASY_STEPS_BEFORE_GROWTH and dt < MAX_DT_S:
        dt = min(MAX_DT_S, dt * DT_GROWTH_FACTOR)
        easy_step_count = 0

    if time + 1.0e-10 >= next_measure_time:
        current_mass = copper_mass()
        relative_mass_change = (current_mass - initial_mass) / initial_mass

        last_metrics = write_measurement(
            step,
            time,
            relative_mass_change,
            accepted_dt=accepted_dt,
            newton_iterations=nonlinear_iterations,
            convergence_reason=convergence_reason,
            residual_norm=residual_norm,
        )
        last_measured_step = step

        if MPI.COMM_WORLD.rank == 0:
            print(
                f"Step {step:5d} | "
                f"time={time:.3f} s | dt={float(dt_constant.value):.5g} s | "
                f"max={last_metrics['maximum_bulge_nm']:.3f} nm | "
                f"positive avg={last_metrics['full_pad_average_nm']:.3f} nm | "
                f"signed avg={last_metrics['signed_full_pad_average_nm']:.3f} nm | "
                f"raised avg={last_metrics['positive_region_average_nm']:.3f} nm | "
                f"edge dInit={last_metrics['edge_change_from_initial_nm']:.3f} nm | "
                f"angle={last_metrics['apparent_contact_angle_deg']:.2f} deg | "
                f"area={last_metrics['positive_area_fraction']:.3f} | "
                f"Newton={nonlinear_iterations:2d} | "
                f"res={residual_norm:.3e} | "
                f"mass={relative_mass_change:.3e}"
            )

        while next_measure_time <= time + 1.0e-10:
            next_measure_time += measure_interval_s

    if time + 1.0e-10 >= next_save_time:
        update_output()
        xdmf.write_function(phi_output, time)
        last_saved_step = step

        while next_save_time <= time + 1.0e-10:
            next_save_time += save_interval_s




def moving_average_1d(values, window):
    """Centered moving average with edge padding. Used only for curvature diagnostics."""
    if window <= 1:
        return np.asarray(values, dtype=float).copy()
    if window % 2 == 0:
        raise ValueError("Curvature smoothing window must be odd.")
    pad = window // 2
    values_padded = np.pad(np.asarray(values, dtype=float), pad_width=pad, mode="edge")
    kernel = np.ones(window, dtype=float) / float(window)
    return np.convolve(values_padded, kernel, mode="valid")


def compute_physical_surface_curvature(radius_dimless, height_dimless, gamma_j_m2, smooth_window=7):
    """
    Compute physical axisymmetric curvature from extracted z(r), not from phi derivatives.

    Axisymmetric surface z(r):
        kappa = z''/(1 + z'^2)^(3/2) + z'/(r * sqrt(1 + z'^2))
    At r=0, z'/r is replaced by z'' for a smooth axisymmetric surface.

    Inputs are dimensionless FEniCSx coordinates scaled by L0.
    Curvature output is 1/m; gamma*kappa is J/m^3.
    """
    r_m = np.asarray(radius_dimless, dtype=float) * L0_m
    z_m_raw = np.asarray(height_dimless, dtype=float) * L0_m

    order = np.argsort(r_m)
    r_m = r_m[order]
    z_m_raw = z_m_raw[order]

    # Remove duplicate radii, if any, by averaging the corresponding heights.
    profile_df = pd.DataFrame({"r_m": r_m, "z_m": z_m_raw}).groupby("r_m", as_index=False).mean()
    r_m = profile_df["r_m"].to_numpy(dtype=float)
    z_m_raw = profile_df["z_m"].to_numpy(dtype=float)

    if len(r_m) < max(7, smooth_window + 2):
        raise RuntimeError(f"Need more surface points for curvature calculation. Got {len(r_m)}.")

    z_m_smooth = moving_average_1d(z_m_raw, smooth_window)

    dz_dr = np.gradient(z_m_smooth, r_m, edge_order=2)
    d2z_dr2 = np.gradient(dz_dr, r_m, edge_order=2)

    kappa_meridional = d2z_dr2 / ((1.0 + dz_dr**2) ** 1.5)

    kappa_azimuthal = np.empty_like(kappa_meridional)
    nonzero_radius = np.abs(r_m) > 1.0e-15
    kappa_azimuthal[nonzero_radius] = (
        dz_dr[nonzero_radius]
        / (r_m[nonzero_radius] * np.sqrt(1.0 + dz_dr[nonzero_radius]**2))
    )
    # Smooth limit at r=0: dz/dr ~ z'' r, so dz/dr/r -> z''.
    kappa_azimuthal[~nonzero_radius] = d2z_dr2[~nonzero_radius]

    kappa_axisymmetric = kappa_meridional + kappa_azimuthal
    gamma_kappa = gamma_j_m2 * kappa_axisymmetric

    curvature_df = pd.DataFrame(
        {
            "radius_nm": r_m / 1.0e-9,
            "z_raw_nm": z_m_raw / 1.0e-9,
            "z_smoothed_nm": z_m_smooth / 1.0e-9,
            "dz_dr_dimensionless": dz_dr,
            "d2z_dr2_1_per_m": d2z_dr2,
            "kappa_meridional_1_per_m": kappa_meridional,
            "kappa_azimuthal_1_per_m": kappa_azimuthal,
            "kappa_axisymmetric_1_per_m": kappa_axisymmetric,
            "gamma_kappa_axisymmetric_J_m3": gamma_kappa,
            "abs_gamma_kappa_axisymmetric_J_m3": np.abs(gamma_kappa),
        }
    )

    curvature_stats = {
        "curvature_smoothing_window_points": smooth_window,
        "curvature_points": int(len(curvature_df)),
        "kappa_min_1_m": float(np.nanmin(kappa_axisymmetric)),
        "kappa_mean_1_m": float(np.nanmean(kappa_axisymmetric)),
        "kappa_max_1_m": float(np.nanmax(kappa_axisymmetric)),
        "abs_kappa_mean_1_m": float(np.nanmean(np.abs(kappa_axisymmetric))),
        "abs_kappa_max_1_m": float(np.nanmax(np.abs(kappa_axisymmetric))),
        "gamma_kappa_min_J_m3": float(np.nanmin(gamma_kappa)),
        "gamma_kappa_mean_J_m3": float(np.nanmean(gamma_kappa)),
        "gamma_kappa_max_J_m3": float(np.nanmax(gamma_kappa)),
        "abs_gamma_kappa_mean_J_m3": float(np.nanmean(np.abs(gamma_kappa))),
        "abs_gamma_kappa_max_J_m3": float(np.nanmax(np.abs(gamma_kappa))),
    }

    return curvature_df, curvature_stats


# 20. FINAL OUTPUT

final_mass = copper_mass()

final_mass_change = (
    (final_mass - initial_mass)
    / initial_mass
)

if last_measured_step != step:
    last_metrics = write_measurement(
        step,
        time,
        final_mass_change,
        accepted_dt=accepted_dt,
        newton_iterations=nonlinear_iterations,
        convergence_reason=convergence_reason,
        residual_norm=residual_norm,
    )

if last_saved_step != step:
    update_output()

    xdmf.write_function(
        phi_output,
        time,
    )

xdmf.close()
csv_file.close()

# Export the final upper Cu surface for direct comparison with TEM/AFM data.
update_output()
final_r, final_z = extract_interface_profile()
final_profile_path = output_folder / f"stage15_{run_label}_final_phi_0p5_surface_profile.csv"

if MPI.COMM_WORLD.rank == 0:
    with open(final_profile_path, "w", newline="") as final_profile_file:
        writer = csv.writer(final_profile_file)
        writer.writerow([
            "radius_nm",
            "surface_height_nm",
            "displacement_from_flat_reference_nm",
            "change_from_initial_profile_nm",
        ])
        initial_z = initial_surface_height(final_r)
        for r_value, z_value, z_initial in zip(final_r, final_z, initial_z):
            writer.writerow([
                r_value * L0_nm,
                z_value * L0_nm,
                (z_value - reference_height) * L0_nm,
                (z_value - z_initial) * L0_nm,
            ])


curvature_path = output_folder / f"stage15_{run_label}_physical_curvature.csv"
curvature_stats = None

if MPI.COMM_WORLD.rank == 0:
    curvature_df, curvature_stats = compute_physical_surface_curvature(
        final_r,
        final_z,
        GAMMA_CU_J_M2,
        smooth_window=7,
    )
    curvature_df.to_csv(curvature_path, index=False)

summary_path = output_folder / f"stage15_{run_label}_summary.txt"

if MPI.COMM_WORLD.rank == 0:
    with open(summary_path, "w") as summary:
        summary.write(
            "STAGE 15 ANSYS-THERMAL-INITIALIZED + ANSYS-W SURFACE-DIFFUSION MODEL\n"
        )
        summary.write(
            "ANSYS mechanical field W is shifted, Cu-masked, and interpolated to the FEniCSx mesh.\n"
        )
        summary.write(
            "Stage 15 initializes phi from the ANSYS thermally displaced Cu top surface, then evolves surface diffusion with the physical wall energy and optional ANSYS W coupling.\n\n"
        )
        summary.write(
            f"Temperature_C: {anneal_temperature_C}\n"
        )
        summary.write(
            f"Dielectric_material_label: {DIELECTRIC_MATERIAL_LABEL}"
        )
        summary.write(
            f"Pad_diameter_nm: {2.0 * pad_radius * L0_nm}\n"
        )
        summary.write(
            f"Pad_height_nm: {pad_height * L0_nm}\n"
        )
        summary.write(
            f"Recess_nm: {recess * L0_nm}\n"
        )
        summary.write(
            f"Reference_wall_energy_equilibrium_angle_deg: "
            f"{reference_wall_equilibrium_angle_deg}\n"
        )
        summary.write(f"Wall_energy_strength: {WALL_ENERGY_STRENGTH}\n")
        summary.write(
            f"Active_wall_energy_equivalent_angle_deg: "
            f"{active_wall_equilibrium_angle_deg}\n"
        )
        summary.write(
            f"Reference_liner_void_minus_Cu_liner_J_m2: "
            f"{liner_void_minus_cu_liner_J_m2}\n"
        )
        summary.write(
            f"Active_liner_void_minus_Cu_liner_J_m2: "
            f"{active_wall_energy_difference_J_m2}\n"
        )
        summary.write(
            "Finite_wall_mobility_enabled: False\n"
        )
        summary.write(
            "Wall_relaxation_time_s: not_used\n"
        )
        summary.write(f"Edge_mobility_limiter_enabled: {EDGE_MOBILITY_LIMITER_ENABLED}\n")
        summary.write(f"Edge_contact_line_mobility_factor: {EDGE_CONTACT_LINE_MOBILITY_FACTOR}\n")
        summary.write(f"Edge_mobility_start_radius_nm: {EDGE_MOBILITY_START_RADIUS_NM}\n")
        summary.write(f"Initial_profile_mode: {INITIAL_PROFILE_MODE}\n")
        summary.write(f"Center_dish_nm: {CENTER_DISH_NM}\n")
        summary.write(f"Edge_rolloff_nm: {EDGE_ROLLOFF_NM}\n")
        summary.write(f"Edge_rolloff_start_radius_nm: {EDGE_ROLLOFF_START_RADIUS_NM}\n")
        summary.write(f"Initial_profile_volume_normalize: {INITIAL_PROFILE_VOLUME_NORMALIZE}\n")
        summary.write(
            f"Initial_profile_raw_area_average_correction_nm: "
            f"{INITIAL_PROFILE_AREA_AVERAGE_CORRECTION * L0_nm}\n"
        )
        summary.write(f"Initial_thermal_expansion_shift_nm: {INITIAL_THERMAL_EXPANSION_SHIFT_NM}\n")
        summary.write(f"ANSYS_thermal_displacement_initial_profile_enabled: {ANSYS_THERMAL_DISPLACEMENT_INITIAL_PROFILE_ENABLED}\n")
        summary.write(f"ANSYS_thermal_displacement_path: {ANSYS_THERMAL_DISPLACEMENT_PATH}\n")
        summary.write(f"ANSYS_thermal_displacement_surface_tolerance_nm: {ANSYS_THERMAL_DISPLACEMENT_SURFACE_TOLERANCE_NM}\n")
        summary.write(f"ANSYS_thermal_displacement_value_scale: {ANSYS_THERMAL_DISPLACEMENT_VALUE_SCALE}\n")
        summary.write(f"ANSYS_Cu_bottom_y_nm: {ANSYS_CU_BOTTOM_Y_NM}\n")
        summary.write(f"ANSYS_thermal_displacement_zero_reference: {ANSYS_THERMAL_DISPLACEMENT_ZERO_REFERENCE}\n")
        summary.write(f"ANSYS_thermal_displacement_min_nm: {np.min(ANSYS_THERMAL_SURFACE_UY) * L0_nm}\n")
        summary.write(f"ANSYS_thermal_displacement_mean_nm: {np.mean(ANSYS_THERMAL_SURFACE_UY) * L0_nm}\n")
        summary.write(f"ANSYS_thermal_displacement_max_nm: {np.max(ANSYS_THERMAL_SURFACE_UY) * L0_nm}\n")
        summary.write(f"Initial_center_surface_height_nm: {initial_surface_height(np.array([0.0]))[0] * L0_nm}\n")
        summary.write(f"Initial_edge_surface_height_nm: {initial_surface_height(np.array([pad_radius]))[0] * L0_nm}\n")
        summary.write(
            f"D_surface_raw_m2_s: {D_surface_raw}\n"
        )
        summary.write(
            f"Diffusivity_calibration_factor_C_D: "
            f"{DIFFUSIVITY_CALIBRATION_FACTOR}\n"
        )
        summary.write(
            f"D_surface_effective_m2_s: {D_surface_effective}\n"
        )
        summary.write(
            f"M_gamma_effective_1_J_s: {M_gamma}\n"
        )
        summary.write(
            f"Gamma_Cu_J_m2: {GAMMA_CU_J_M2}\n"
        )
        summary.write(
            f"Omega_Cu_m3_atom: {OMEGA_CU_M3_ATOM}\n"
        )
        summary.write(
            f"Surface_layer_h_nm: {SURFACE_LAYER_THICKNESS_M * 1e9}\n"
        )
        summary.write(
            f"Interface_width_epsilon_nm: {interface_width_nm}\n"
        )
        summary.write(
            f"initial_dt_s: {INITIAL_DT_S}\n"
            f"minimum_dt_s: {MIN_DT_S}\n"
            f"maximum_dt_s: {MAX_DT_S}\n"
        )
        summary.write(
            f"final_time_s: {final_time}\n"
        )
        summary.write(f"mesh_cells_radial: {MESH_CELLS_RADIAL}\n")
        summary.write(f"mesh_cells_vertical: {MESH_CELLS_VERTICAL}\n")
        summary.write(f"snes_rtol: {SNES_RTOL}\n")
        summary.write(f"snes_atol: {SNES_ATOL}\n")
        summary.write(f"snes_stol: {SNES_STOL}\n")
        summary.write(f"snes_max_it: {SNES_MAX_IT}\n")
        summary.write(
            f"Positive_tolerance_nm: {positive_tolerance_nm}\n"
        )
        summary.write(f"ANSYS_W_enabled: {ANSYS_W_ENABLED}\n")
        summary.write(f"ANSYS_W_csv_path: {ANSYS_W_CSV_PATH}\n")
        summary.write(f"ANSYS_W_clip_negative: {ANSYS_W_CLIP_NEGATIVE}\n")
        summary.write(f"ANSYS_W_mask_to_initial_Cu_pad: {ANSYS_W_MASK_TO_INITIAL_CU_PAD}\n")
        summary.write(f"ANSYS_W_mask_padding_nm: {ANSYS_W_MASK_PADDING_NM}\n")
        summary.write(f"ANSYS_W_discard_source_points_outside_mask: {ANSYS_W_DISCARD_SOURCE_POINTS_OUTSIDE_MASK}\n")
        summary.write(f"ANSYS_W_zero_dof_values_outside_mask: {ANSYS_W_ZERO_DOF_VALUES_OUTSIDE_MASK}\n")
        summary.write(f"ANSYS_W_outside_domain_value: {ANSYS_W_OUTSIDE_DOMAIN_VALUE}\n")
        summary.write(f"Mechanical_W_multiplier: {MECHANICAL_W_MULTIPLIER}\n")
        summary.write(f"Interpolated_W_min_J_m3: {float(np.min(W_field.x.array))}\n")
        summary.write(f"Interpolated_W_mean_J_m3: {float(np.mean(W_field.x.array))}\n")
        summary.write(f"Interpolated_W_max_J_m3: {float(np.max(W_field.x.array))}\n\n")
        if curvature_stats is not None:
            summary.write("# Physical curvature diagnostic from extracted final z(r) surface\n")
            summary.write(f"Physical_curvature_csv: {curvature_path}\n")
            for key, value in curvature_stats.items():
                summary.write(f"{key}: {value}\n")
            summary.write("\n")
        summary.write(
            f"Final_maximum_bulge_nm: "
            f"{last_metrics['maximum_bulge_nm']}\n"
        )
        summary.write(
            f"Final_full_pad_average_positive_bulge_nm: "
            f"{last_metrics['full_pad_average_nm']}\n"
        )
        summary.write(
            f"Final_signed_full_pad_area_average_nm: "
            f"{last_metrics['signed_full_pad_average_nm']}\n"
        )
        summary.write(
            f"Final_signed_unweighted_point_average_nm: "
            f"{last_metrics['signed_point_average_nm']}\n"
        )
        summary.write(
            f"Final_positive_region_average_nm: "
            f"{last_metrics['positive_region_average_nm']}\n"
        )
        summary.write(
            f"Final_center_displacement_nm: "
            f"{last_metrics['center_displacement_nm']}\n"
        )
        summary.write(
            f"Final_edge_displacement_from_flat_reference_nm: "
            f"{last_metrics['edge_displacement_nm']}\n"
        )
        summary.write(
            f"Final_maximum_change_from_initial_profile_nm: "
            f"{last_metrics['maximum_change_from_initial_nm']}\n"
        )
        summary.write(
            f"Final_center_change_from_initial_profile_nm: "
            f"{last_metrics['center_change_from_initial_nm']}\n"
        )
        summary.write(
            f"Final_edge_change_from_initial_profile_nm: "
            f"{last_metrics['edge_change_from_initial_nm']}\n"
        )
        summary.write(
            f"Final_apparent_contact_angle_deg: "
            f"{last_metrics['apparent_contact_angle_deg']}\n"
        )
        summary.write(
            f"Final_positive_area_fraction: "
            f"{last_metrics['positive_area_fraction']}\n"
        )
        summary.write(
            f"Final_relative_mass_change: "
            f"{final_mass_change}\n"
        )

    print("\nStage 15 simulation complete.")
    print(
        f"Final maximum bulge: "
        f"{last_metrics['maximum_bulge_nm']:.4f} nm"
    )
    print(
        f"Final positive-only full-pad equivalent bulge: "
        f"{last_metrics['full_pad_average_nm']:.4f} nm"
    )
    print(
        f"Final signed full-pad area average: "
        f"{last_metrics['signed_full_pad_average_nm']:.4f} nm"
    )
    print(
        f"Final signed unweighted point average: "
        f"{last_metrics['signed_point_average_nm']:.4f} nm"
    )
    print(
        f"Final raised-region average: "
        f"{last_metrics['positive_region_average_nm']:.4f} nm"
    )
    print(
        f"Final mass change: {final_mass_change:.3e}"
    )
    print(f"Summary: {summary_path}")
    print(f"Final surface profile: {final_profile_path}")
    print(f"Physical curvature: {curvature_path}")
    if curvature_stats is not None:
        print(
            f"Mean |gamma*kappa|: "
            f"{curvature_stats['abs_gamma_kappa_mean_J_m3']:.6e} J/m^3"
        )
        print(
            f"Max |gamma*kappa|: "
            f"{curvature_stats['abs_gamma_kappa_max_J_m3']:.6e} J/m^3"
        )
    print(
        "Open ParaView with:\n"
        f"open -a ParaView {xdmf_path}"
    )