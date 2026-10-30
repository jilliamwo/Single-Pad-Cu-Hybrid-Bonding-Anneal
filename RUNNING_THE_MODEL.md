# Running the model

## Runtime limit for future work

The completed 400 × 880 (0.625 nm) benchmark took multiple hours. At the user's
request, do not repeat it or propose still finer uniform meshes. Its configuration
`configs/model_mesh_finer_dt025.json` is retained for reproducibility only.
Use 100 × 220 for routine screening and at most 200 × 440 for selected accuracy
checks; reuse the completed benchmark to quantify sensitivity. Prefer saved-field
analysis and short tests before adding simulations. The finer working mesh is
still more costly, so this limit does not justify automatically launching long runs.


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
- Anneal temperature 250 °C. The 22 °C stress-free reference is an
  **unverified calibration assumption**, not an annealing run or a confirmed
  ANSYS setting. Change it with `--reference-c VALUE` once verified; the model
  records the value used. Changing it requires revisiting mechanical calibration.
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
archive. A PETSc exception first restores the accepted state and rebuilds the
solver once at the same time step. A repeated exception stops the run; ordinary
nonlinear non-convergence instead reduces the time step. Full PETSc messages
are kept in `solver_errors.json`, and `step_history.csv` records error codes.
On a rejected step the accepted surface and mass diagnostics are saved immediately.
At terminal solver failure, `last_accepted_state.npz` also saves mixed field values
for local diagnosis; it is not a supported restart file and is excluded from the
small ZIP. The last accepted step is recorded in `step_history.csv`; profiles are
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
python new/validate_model_recovery.py
```

This checks the coordinate translation, Cu-side sampling and void mask, and
compares the residual and Jacobian with the archived weak form. Short native
plasticity and ANSYS-driven runs were also exercised, including comparison and
ZIP creation. Full 100 s and 7000 s simulations are left for you to run.
The next numerical study after the first comparison is mesh/time-step sensitivity;
full feedback coupling comes after that baseline is understood.

## Rerun after the saved 44-second failure

The existing `runs/native_100s` stopped at 44 s with PETSc error 73 (a vector
read lock), not a reported SNES convergence failure. The original retry logic
reused the same solver while reducing dt and repeated the error. Its surface
CSV only reached t=0 because the measurement interval was 50 s; do not interpret
that file as the 44 s shape. Those original results have been preserved.

The new recovery path passed short fault-injection tests. The actual error was
not reproduced on a coarse 60-step test, so recovery of the full-size case still
needs this rerun. Use fresh directories:

```sh
python run_model.py run --duration-s 100 --output runs/native_100s_retry
python run_model.py run --source ansys --duration-s 100 --output runs/ansys_100s_retry
python run_model.py compare runs/native_100s_retry runs/ansys_100s_retry --output runs/comparison_100s_retry
```

Send the three `analysis.zip` files. If either simulation fails, send its ZIP
immediately; comparison requires two completed runs.

## Phase-field sensitivity study after the successful 100 s comparison

Use `runs/native_100s_retry2` as the existing baseline (dt=1 s, phase mesh
100 × 220). The study configurations keep mechanics at 80 radial subdivisions,
80 per vertical main layer, and 40 thermal increments; material properties,
interface width, wall energy, diffusivity, and the assumed reference temperature
are unchanged. Each run recomputes the same mechanical problem. No long study
simulations were run automatically.

First check time steps on the same phase mesh:

```sh
conda activate fenicsx-env
python run_model.py run --config configs/model_dt05.json --output runs/native_dt05_100s
python run_model.py run --config configs/model_dt025.json --output runs/native_dt025_100s
python run_model.py compare runs/native_100s_retry2 runs/native_dt05_100s --study time --output runs/time_1_vs_05
python run_model.py compare runs/native_dt05_100s runs/native_dt025_100s --study time --output runs/time_05_vs_025
```

Both initial and maximum dt are capped at the specified value, so automatic
step growth cannot undo refinement. Failed nonlinear steps may still reduce dt;
we will inspect accepted-step histories when interpreting the study.
Send the two new run archives and two comparison archives, or tell me when the
workspace results are ready. We will check whether the surface differences shrink
from the first pair to the second, especially at the Cu edge. Pairwise agreement
alone is not a proof of convergence.

Then check phase-mesh sensitivity at dt=0.5 s:

```sh
python run_model.py run --config configs/model_mesh_fine_dt05.json --output runs/native_meshfine_dt05_100s
python run_model.py compare runs/native_dt05_100s runs/native_meshfine_dt05_100s --study mesh --output runs/mesh_100_vs_200
```

This doubles subdivisions to 200 × 440 (nominal spacing 1.25 nm), with roughly
four times as many phase cells. It leaves the physical interface width at 5 nm.
If the time-step study shows 0.5 s is still too coarse, postpone this run and
prepare both meshes at the smaller step instead.

Mesh refinement changes both the phase discretization and the P1 representation
of the initial surface and frozen energy. Compare initial profiles as well as
final profiles and changes from each initial state. Mechanics itself is not
refined in this experiment. The comparisons report sensitivity, not an automatic
pass/fail or claimed convergence order. We want numerical uncertainty comfortably
below the approximately 0.052 nm native/ANSYS final-surface RMS difference before
interpreting that difference quantitatively.

`--study time` permits only initial/maximum dt differences; `--study mesh` permits
only phase mesh differences (visualization output may differ in either study).
The default `--study inputs` retains the existing mechanical-input comparison.
All studies require completed runs reaching the configured final time.

Validation of comparison isolation, including rejection of changed reference
temperatures and mixed mesh/time changes:

```sh
python new/validate_sensitivity.py
```

## Routine-mesh mechanical-energy ablation

Saved-field analysis shows the corner mismatch persists under denser sampling.
To assess the phase evolution's sensitivity to energy forcing, run:

```sh
python run_model.py run --config configs/model_energy_off_dt05.json --output runs/native_energyoff_dt05_100s
python run_model.py compare runs/native_dt05_100s runs/native_energyoff_dt05_100s --study energy --output runs/energy_on_vs_off
```

This uses the routine 100 × 220 mesh, dt capped at 0.5 s, and 100 s duration.
It retains native plastic mechanics and the same initial deformed shape. Only
energy forcing is disabled. `--source none` is not equivalent. The energy study
mode rejects changes to other physical/numerical settings. Its difference reports
the total effect of energy forcing; it does not directly measure the corner's
contribution or establish which source field is more accurate.

Optional small developer checks for the experimental mass-conserving transfer
(run inside `fenicsx-env`; these do not launch production simulations):

```sh
python new/validate_column_mass_transfer.py
python new/validate_feedback_loop.py --occupancy-mode fractional --transfer-mode column_mass
```

Reports are written under `new/verification_results/`. This transfers diffuse
phase-mass changes to a mass-equivalent mechanics surface, which can differ from
the plotted phi=0.5 contour. It is not enabled by `run_model.py run` yet.

## Short experimental two-way pilot

The separate pilot runner now exposes the fractional, mass-conserving feedback
path for a controlled interval comparison. The ordinary `run_model.py run`
command still uses static mechanics. Run from the repository root:

```sh
conda activate fenicsx-env
python run_feedback_pilot.py --output "runs/feedback_pilot_$(date +%Y%m%d_%H%M%S)"
```

This creates a unique output directory and runs two cases sequentially:

- Physical duration: 2 seconds each; fixed phase timestep: 0.25 seconds.
- Phase mesh: existing practical 100x220; mechanics: existing 80x80 per layer.
- Case `interval_1`: mechanics updates every phase step (0.25 seconds).
- Case `interval_2`: mechanics updates every two phase steps (0.5 seconds).
- Both cases start independently from the same native elastoplastic solution.
- Fractional occupancy, column-mass transfer, and explicit homogeneous
  virgin-at-reference mixing of new Cu history are enabled.

This is experimental. The 22 C stress-free temperature remains an unverified
ANSYS-reference assumption. The study measures update-interval sensitivity, not
physical validation. A two-second duration bounds the number of solves; runtime
still depends on mechanical convergence. It does not use the expensive fine
phase mesh.

On success, `comparison.json` reports final contour differences, sampled energy
differences, plastic-history differences, volume drift, and timings. Both cases'
full diagnostics and a code snapshot are included in the single `analysis.zip`
whose path is printed at completion. Send that archive for analysis. Large field
files are omitted in the routine pilot. On failure or Ctrl-C the runner saves a
failure status and packages available diagnostics; send that archive rather than
increasing mesh resolution or repeatedly changing numerical parameters.

This study intentionally disables timestep reductions: a failed step stops the
pilot so the cases cannot silently use different time discretizations. Existing
nonempty output folders are never overwritten.

For an optional tiny installation/workflow check only:

```sh
python run_feedback_pilot.py --smoke --output "runs/feedback_smoke_$(date +%Y%m%d_%H%M%S)"
```

The smoke mode uses a coarse mesh and 0.1 seconds, and is not a physical result.

The completed 2 s pilot showed small update-interval sensitivity. To extend the
same comparison to 10 s without changing either mesh or timestep:

```sh
conda activate fenicsx-env
python run_feedback_pilot.py --duration-s 10 --output "runs/feedback_10s_$(date +%Y%m%d_%H%M%S)"
```

This runs 40 phase steps per case and 40 versus 20 mechanics updates. Share the
single printed `analysis.zip`. The duration must be a positive multiple of 0.5 s
so both interval cases update mechanics at the final phase time. A 10 s comparison
is the next check; the 2 s result alone does not establish long-time sensitivity.

## Measure the effect of mechanics feedback

The 10 s interval pilot completed without rejected steps and showed only
0.000090 nm maximum surface difference between 0.25 and 0.5 s updates. The next
controlled comparison uses the provisional 0.5 s interval versus frozen initial
energy, with the *same* initial energy mapping in both cases:

```sh
conda activate fenicsx-env
python run_feedback_pilot.py --study feedback_strength --duration-s 10 --output "runs/feedback_effect_10s_$(date +%Y%m%d_%H%M%S)"
```

This runs two 10 s cases on the same practical meshes and fixed 0.25 s phase
steps: `interval_2` updates mechanics 20 times; `frozen` never updates mechanics.
The frozen case still uses the experimental initial column-coordinate energy
mapping, so it is not identical to the legacy static CLI's transfer convention.
Initial surface and energy arrays are checked for exact equality. The comparison
is updated minus frozen; energy relative L2 is normalized by the updated case.
Frozen-case plastic history remains at initialization, and its contour-to-mass
geometry discrepancy is not evaluated (reported null). Both cases and source
snapshots are included in the printed analysis.zip. Send that archive afterward.

The completed 10 s feedback-effect study found a 0.01056 nm edge difference
(0.001763 nm RMS), larger than update-interval sensitivity but small compared
with total edge motion. Extend the same controlled comparison to 100 s:

```sh
conda activate fenicsx-env
python run_feedback_pilot.py --study feedback_strength --duration-s 100 --output "runs/feedback_effect_100s_$(date +%Y%m%d_%H%M%S)"
```

This uses the same 100x220 phase mesh and 80x80 mechanics mesh, 400 fixed phase
steps per case, and 200 mechanics updates in the feedback case. No refinement.
A rough estimate from the 10 s run is 18-20 minutes total; this is not a runtime
guarantee. Send the printed analysis.zip on completion or failure.

Saved paired studies can also be reviewed without running any simulation:

```sh
/opt/anaconda3/bin/python new/review_feedback_study.py runs/feedback_effect_10s_20261029_194735
```

This rechecks equal initial arrays and writes review.json and feedback_review.png.
It needs numpy and matplotlib; the workspace's base Anaconda Python has both.
