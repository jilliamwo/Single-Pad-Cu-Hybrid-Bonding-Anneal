# Running the model

The active entry point is `run_model.py`. Mechanics lives in
`new/thermoelastic_solver.py` and `new/plasticity.py`; the surface-diffusion equations
live in `new/phase_field_solver.py`. `new/static_coupling.py` connects them.
The old script is retained under `old/` as an archived reference. Active commands,
output names, and solver names use **the model**.

## Start here

Run from a terminal:

```sh
cd "/Users/williamjo/VSCode/Single-Pad-Cu-Hybrid-Bonding-Anneal"
conda activate fenicsx-env
python run_model.py doctor
python run_model.py run --preset smoke --output runs/smoke
```

`doctor` reports the Python environment and verifies scientific imports.
SciPy was installed in `fenicsx-env` to support ANSYS linear interpolation.
If `conda activate` is unavailable in your shell, replace `python` in these
commands with `/opt/anaconda3/envs/fenicsx-env/bin/python`.
The workflow targets DOLFINx 0.10.0 and runs on **one MPI rank**; do not use `mpirun`.

The smoke preset uses a coarse mesh, four thermal increments, and only 0.1 s
of phase-field evolution. It checks installation and data flow; its shapes are
not converged physical predictions. All output directories must be new or empty.
Choose a new name when repeating a run; existing results are never overwritten.

## First simulations to send for analysis

Use the same phase-field settings with two sources of mechanical input:

```sh
python run_model.py run --config configs/model_baseline.json --duration-s 100 --output runs/native_100s
python run_model.py run --config configs/model_baseline.json --source ansys --duration-s 100 --output runs/ansys_100s
python run_model.py compare runs/native_100s runs/ansys_100s --output runs/comparison_100s
```

The first command solves native mechanics with Cu plasticity, then evolves the
surface for 100 s. The second imports displacement and elastic energy from
`old/ansys-results/250_500nm_3.5_SiO2.csv`, then runs the same phase-field model.
The comparison reports initial and final surface differences, plus differences
in the surface *change* from each run's own initial condition.

Send these three files, identified by their parent folders:

- `runs/native_100s/analysis.zip`
- `runs/ansys_100s/analysis.zip`
- `runs/comparison_100s/analysis.zip`

These contain settings, source snapshots, mechanics diagnostics, surface profiles,
solver histories, mass conservation, and logs. The comparison alone is not enough
for diagnosing a mismatch. Large visualization files are excluded from the ZIPs.

I will analyze transfer consistency, initial/final shape differences, mass drift,
phase overshoot, nonlinear convergence, and sensitivity needed before longer runs.

## Full-duration runs

After reviewing the short comparison, the equivalent 7000 s commands are:

```sh
python run_model.py run --config configs/model_baseline.json --output runs/native_7000s
python run_model.py run --config configs/model_baseline.json --source ansys --output runs/ansys_7000s
python run_model.py compare runs/native_7000s runs/ansys_7000s --output runs/comparison_7000s
```

These start fresh; the 100 s runs are not restart checkpoints. A terminal process
must remain running. ParaView can open `phi.xdmf` for time-dependent phase fields
and `static_energy.xdmf` for the transferred energy. Keep each XDMF file together
with its HDF5 file. Add `--no-fields` if only CSV/JSON diagnostics are needed.

## Controls and assumptions

`configs/model_baseline.json` contains editable settings. Copy it to a new file
for a controlled experiment and pass `--config your_file.json`. Explicit command
options override the JSON, which overrides the selected preset. Each run records
the fully resolved settings in its own `config.json`.

The baseline freezes the previously studied SiO2 case:

- Cu radius 250 nm; dielectric outer radius 500 nm; base thickness 500 nm;
  Cu top at 996.5 nm; dielectric top at 1000 nm.
- Axis `ur=0`, bottom `ur=uz=0`, outer boundary `ur=0`.
- Anneal temperature 250 °C and stress-free reference 22 °C.
- Cu J2 plasticity: initial yield 172.3 MPa and tangent modulus 517 MPa.
  Dielectric remains elastic; no unmeasured dielectric plasticity is introduced.
- Mechanics: 80 radial subdivisions and 80 per main vertical layer, plus recess
  rows; 40 quasi-static thermal load increments.
- Phase field: 100 × 220 triangular-cell subdivisions over 250 × 550 nm
  (nominal spacings 2.5 × 2.5 nm); interface-width parameter 5 nm.
- Initial/maximum phase time step 1 s, minimum 0.03125 s; failed nonlinear steps
  restore the last accepted state and reduce the time step.
- Surface diagnostics every 50 s; visualization every 100 s; final output always
  saved. Adaptive stepping can place intermediate outputs just after these times.
- Existing wall energy −0.651 J/m² and diffusivity factor 0.1 are retained.
  The mechanical calibration does **not** validate these phase-field parameters.

This is **static one-way coupling**. The model first ramps prescribed uniform
thermal loading, then transfers displacement and elastic energy once. Phase time
starts after that mechanical initialization. It does not yet solve transient
heat conduction, update mechanical equilibrium as the surface changes, introduce
contact, or model creep. It therefore does not reproduce every aspect of an
ANSYS transient thermal + transient structural analysis.

Mechanical coordinates use the bottom of the dielectric as z=0; phase coordinates
use the Cu bottom as z=0. Transfer subtracts 500 nm. The native energy is sampled
on the Cu side of material interfaces and set to zero outside reference Cu, then
interpolated into the phase-field P1 space. This is a frozen reference-coordinate
field, not a conservative remap or an energy recomputed on the evolving shape.
P1 interpolation smooths DG0 values near jumps. Its sensitivity remains to be
checked. The ANSYS path retains the archived initial-surface mask and linear
interpolation convention. `--ansys-interpolation nearest` is available as an
explicit sensitivity experiment, not an automatic fallback.

Reported center height relative to the dielectric subtracts the **deformed
inner edge of the dielectric surface**, not its nominal height. Mechanical
vertical displacement is already included in the initial phase surface and is
not added a second time. The extracted 0.5 contour can differ slightly from the
analytic initial shape due to mesh resolution; both are saved in
`initial_surface.csv`. Surface extraction scans phase z=425–550 nm; missing
crossings are recorded and prevent a misleading automated comparison.

A CSV does not encode its thermal history or constitutive settings. Changing
`temperature_C` does not change imported ANSYS mechanics. Supply a matching CSV
and verify its provenance before doing other-temperature comparisons. These
commands deliberately start with the calibrated SiO2 case; they do not establish
validation across all materials or temperatures.

For a controlled mechanics ablation, use `--no-plasticity`. `--source none` starts
from the nominal recessed shape with zero mechanical energy. These change the
physical input and should be treated as separate experiments.

## Interrupted or failed runs

A Python-level failure records `status.json`, a traceback in `run.log`, and an
archive. The last accepted step is recorded in `step_history.csv`; profiles are
available through the last measurement. After a force quit or machine restart,
rebuild the archive with:

```sh
python run_model.py package runs/native_100s
```

Send that archive even if the run failed. Packaging only needs standard Python.
A process killed outside Python may leave status marked `running`; the logs and
step history help identify where it stopped. Restarting into the same directory
is intentionally refused; use a fresh output directory.

## Implementation verification

```sh
python new/validate_model.py
```

This checks the coordinate translation, Cu-side sampling and void mask, and
compares the residual and Jacobian with the archived weak form. Short native
plasticity and ANSYS-driven runs were also exercised, including comparison and
ZIP creation. Full 100 s and 7000 s simulations are left for you to run.
The next numerical study after the first comparison is mesh/time-step sensitivity;
full feedback coupling comes after that baseline is understood.
