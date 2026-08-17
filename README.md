# Coupled ANSYS–FEniCSx Modeling of Cu Pad Deformation in Hybrid Bonding

## Overview

This project models the evolution of nanoscale copper (Cu) pads during **Cu–Cu hybrid-bonding anneal**. The goal is to understand how process and material parameters influence Cu bulge-out and the surface morphology that ultimately determines whether opposing Cu pads can make contact.

The model combines:

1. **ANSYS finite element analysis (FEA)** — calculates thermoelastic deformation and strain-energy density.
2. **FEniCSx phase-field modeling** — simulates Cu surface evolution driven by surface diffusion, surface energy, and the imported mechanical fields.

The baseline geometry is a **500 nm-diameter Cu pad with a 3.5 nm recess**. Parameter studies include:

- Annealing temperature
- Cu surface diffusivity
- Dielectric material/properties
- Initial Cu surface profile, including dishing, protrusion, and edge roll-off

---

## Model Workflow

```text
ANSYS Transient Thermal
        ↓
ANSYS Static/Transient Structural
        ↓
Export displacement + strain-energy density CSV
        ↓
Stage 15 FEniCSx phase-field model
        ↓
Surface diffusion / Cu morphology evolution
        ↓
CSV + XDMF + summary outputs
        ↓
ParaView / plotting / post-processing
```

The ANSYS result provides the **thermoelastic initial condition** and mechanical strain-energy density \(W\). FEniCSx then evolves the Cu surface during the anneal.

---

# 1. ANSYS Simulation

## Geometry

The model uses a **2D axisymmetric Cu pad embedded in a dielectric**.

Baseline Cu dimensions:

- Diameter: **500 nm**
- Radius: **250 nm**
- Height: **500 nm**
- Initial recess: **3.5 nm**

The left boundary represents the axis of symmetry.

## Analysis

The ANSYS workflow consists of:

### Transient Thermal

Apply the desired annealing temperature and calculate the temperature field.

### Structural Analysis

Import the thermal solution into the structural model and calculate the resulting thermoelastic deformation and stress.

The quantities required by Stage 15 are:

- Node position
- **Directional deformation in the vertical direction**
- Stress/elastic strain information used to obtain **strain-energy density \(W\)**

where

$$
W = \frac{1}{2}\boldsymbol{\sigma}:\boldsymbol{\varepsilon}^{el}
$$

The vertical displacement is used to modify the initial Cu surface, while \(W\) provides the mechanical contribution to the phase-field model.

## Exporting ANSYS Results

Export the required nodal results into a combined CSV and place it in:

```text
fenicsx/
└── ansys-results/
    └── your_result_file.csv
```

For example:

```text
ansys-results/250_500nm_3.5_SiO2.csv
```

Stage 15 must be pointed to the corresponding CSV before running.

---

# 2. Installing FEniCSx

FEniCSx is most easily installed using **Conda/Mamba**.

Install Miniconda or Miniforge if a Conda installation is not already available.

Create a FEniCSx environment:

```bash
conda create -n fenicsx-env
conda activate fenicsx-env
conda install -c conda-forge fenics-dolfinx mpich pyvista pandas scipy matplotlib
```

Activate this environment each time before running the model:

```bash
conda activate fenicsx-env
```

Confirm that FEniCSx loads correctly:

```bash
python -c "import dolfinx; print(dolfinx.__version__)"
```

---

# 3. Running Stage 15

Navigate to the directory containing the Stage 15 Python script:

```bash
cd path/to/hybrid-bonding/fenicsx
```

Activate the environment:

```bash
conda activate fenicsx-env
```

Run the model:

```bash
python ch_pad_stage15_with_curvature.py
```

The exact filename may differ depending on the version stored in the repository.

Stage 15:

1. Reads the ANSYS result CSV.
2. Extracts the Cu top-surface vertical thermal displacement.
3. Interpolates the Cu-masked ANSYS \(W\) field onto the FEniCSx mesh.
4. Constructs the initial Cu geometry.
5. Runs the phase-field surface-diffusion simulation.
6. Tracks Cu surface evolution over the annealing time.
7. Calculates bulge-out, surface profiles, curvature, and other summary metrics.

---

# 4. Important Stage 15 Parameters

Most simulations can be changed through the parameter/settings section near the beginning of the Stage 15 script.

## Geometry

Typical baseline:

```python
pad_radius = 250e-9
pad_height = 500e-9
recess = 3.5e-9
```

Change these parameters to investigate pad diameter, height, or CMP recess.

---

## Annealing Temperature

Set the desired temperature using the temperature variable in Stage 15, for example:

```python
anneal_temperature_C = 250
```

The temperature affects the calculated Cu surface diffusivity through its temperature dependence.

When changing temperature, use the corresponding ANSYS result file.

Example:

```text
200°C → ANSYS 200°C result
250°C → ANSYS 250°C result
300°C → ANSYS 300°C result
```

---

## Surface Diffusivity

Stage 15 calculates the raw Cu surface diffusivity and applies a calibration factor:

\[
D_s^{eff}=C_DD_s.
\]

The baseline simulations use:

```python
C_D = 0.1
```

At 250°C, this corresponds to approximately:

```text
Raw Ds       = 5.86 × 10^-14 m²/s
Effective Ds = 5.86 × 10^-15 m²/s
```

Change `C_D` to perform the surface-diffusivity sweep.

Values investigated in this project include:

```text
0.01
0.03
0.10
0.30
0.50
1.00
```

---

## ANSYS Result Files

Update the paths to match the desired ANSYS simulation:

```python
ANSYS_THERMAL_DISPLACEMENT_PATH = "ansys-results/your_file.csv"
ANSYS_W_CSV_PATH = "ansys-results/your_file.csv"
```

The same combined ANSYS CSV can contain both displacement and \(W\).

---

## Dielectric Label

Set the material/property label used in output filenames and summaries:

```python
DIELECTRIC_MATERIAL_LABEL = "SiO2"
```

Examples used in the material study include:

```text
SiO2
SiCN
SiCOH
BCB
Diamond
```

The actual dielectric mechanical response is calculated in ANSYS. The label in Stage 15 primarily identifies which ANSYS dataset is being used.

---

## Initial Surface Profile

Stage 15 can initialize different CMP surface morphologies.

The profile mode is controlled using:

```python
INITIAL_PROFILE_MODE
```

Additional parameters control the amount of:

```python
CENTER_DISH_NM
CENTER_PROTRUSION_NM
EDGE_ROLLOFF_NM
```

Examples include:

```text
Flat
3 nm dish
5 nm dish
5 nm dish + 3 nm ERO
3 nm protrusion
3 nm protrusion + 3 nm ERO
```

For the 500 nm-diameter pad, the current edge-roll-off region begins near:

```text
r = 225 nm
```

and extends to the pad edge at:

```text
r = 250 nm.
```

---

## Simulation Time

The current simulations typically use:

```python
final_time_s = 7000
```

with adaptive time stepping.

Important time-step settings include:

```python
initial_dt_s
minimum_dt_s
maximum_dt_s
```

Reducing the time step may improve convergence but increases simulation time.

---

## Mesh

The baseline mesh is approximately:

```python
Nr = 100
Nz = 220
```

Increasing mesh resolution improves spatial resolution but substantially increases computational cost.

---

# 5. Output Files

Each Stage 15 simulation creates a labeled results directory so different parameter combinations can be compared.

Important outputs include:

### Summary file

Contains the major simulation parameters and final measurements, including quantities such as:

```text
Final center displacement
Final center change from initial profile
Maximum surface displacement
Positive-region average bulge
Surface diffusivity
Curvature statistics
ANSYS displacement statistics
ANSYS strain-energy-density statistics
```

### Surface-profile CSV

Contains the extracted \(\phi=0.5\) Cu surface:

```text
r
z(r)
```

This is useful for comparing final surface shapes between simulations.

### Time-history CSV

Tracks quantities such as bulge-out as the simulation progresses.

### XDMF Results

The XDMF files can be opened in **ParaView** to visualize the phase field and Cu surface evolution.

---

# 6. Running a Parameter Sweep

Only change the parameter being investigated while keeping the other baseline conditions fixed.

For example, a temperature sweep could use:

```text
200°C
250°C
300°C
350°C
```

with the same geometry, dielectric, initial profile, and diffusivity model.

A surface-diffusivity sweep could instead hold the temperature at 250°C and change only:

```text
C_D = 0.01
C_D = 0.03
C_D = 0.10
C_D = 0.30
C_D = 0.50
C_D = 1.00
```

This controlled approach allows the effect of each parameter to be isolated.

---

# 7. Current Model Interpretation

The simulations indicate that, at the 500 nm scale studied here:

- **Temperature strongly affects Cu evolution**, primarily through its effect on surface-diffusion kinetics.
- **Surface diffusivity controls how quickly the Cu approaches its preferred morphology**, with diminishing changes once the surface approaches equilibrium.
- **Initial CMP profile affects both the starting contact margin and the amount of additional diffusion required.**
- **Dielectric mechanical properties alter thermoelastic displacement and strain energy**, although surface-energy/curvature effects are stronger drivers of the subsequent surface evolution in the current model.
- **Edge recession develops as the Cu surface evolves toward the imposed contact-angle condition.**

These results motivate using temperature, anneal time, recess, initial profile, and surface mobility together as a hybrid-bonding **process window** rather than treating each parameter independently.

---

# 8. Current Limitations

The current model primarily describes **single-pad evolution before Cu–Cu contact**.

It does not yet explicitly model:

- Two opposing Cu pads
- Mechanical contact pressure
- Plastic deformation or creep after contact
- Adhesion/bond formation
- Post-contact stress redistribution
- Trapped edge-void closure
- Oxide or contamination effects

Therefore, the model can predict the morphology leading up to contact but cannot yet determine whether an edge recession will close after the opposing Cu pads touch.

Future work will extend the model toward **Cu–Cu contact, contact-area growth, and residual void prediction**.

---

## Software

- **ANSYS Mechanical / Workbench** — thermoelastic FEA
- **FEniCSx** — phase-field simulation
- **Python** — simulation and post-processing
- **ParaView** — visualization of XDMF results
