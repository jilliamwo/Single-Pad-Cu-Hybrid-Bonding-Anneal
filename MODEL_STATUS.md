# Model status and next work

## Completed

- Native axisymmetric thermoelastoplastic mechanics: Cu J2 plasticity, elastic
  SiO2, stated geometry and constraints. Mechanical outputs calibrated against
  available ANSYS exports; this is not universal validation for all materials.
- Static one-way mechanics → phase-field coupling, with the archived phase-field
  equations preserved and equation/transfer tests passing.
- Paired 100 s native/ANSYS-input runs: final surface RMS difference 0.0520 nm
  on the routine mesh. Both paths shared that mesh's discretization error.
- Time sensitivity: approximately first-order on the coarse mesh; halving dt
  from 0.5 to 0.25 s changes the fine-mesh surface by 0.00109 nm RMS.
- Spatial sensitivity: at matched dt=0.25 s, surface RMS differences shrink from
  0.06544 to 0.01641 nm over three meshes, consistent with approximately second
  order for this output. This is not convergence of every mechanical field.
- Saved-field energy audit: mismatch is strongly concentrated near the Cu outer
  corner. Changing the ANSYS source mask does not explain the sampled discrepancy.

## Current stage

The implementation is working and has useful numerical verification. We are
assessing quantitative accuracy of the coupled model, not declaring the complete
ANSYS transient workflow replaced or the bonding physics fully validated.

## Open

- Separate corner mechanical-field differences from nodal/cellwise energy
  representation and interpolation differences.
- Verify the ANSYS stress-free reference temperature; 22 C remains an assumption.
- Confirm quantitative phase kinetics and wall parameters against physical data.
- Longer-duration behavior and generalization to other materials/temperatures.
- Two-way mechanical feedback, transient heat conduction, contact and creep are
  not implemented in the active static coupling.

## Runtime and workflow constraints

The user runs production simulations; the assistant analyzes and develops code.
Do not request or launch another 400 × 880 or finer uniform phase mesh: the user
reported multiple-hour runtime and explicitly ruled out repeating such fine
meshes. Reuse the benchmark. Use 100 × 220 for routine screening and reserve
200 × 440 for selected checks, considering runtime. Prefer saved-field analysis
and short tests; do not automatically escalate refinement or duration.

## Next

Continue corner-field and transfer diagnostics using saved data. No new user
simulation is required now. Preserve current physics until a specific discrepancy
is isolated. A later short energy-off experiment, if needed, can test whether the
remaining energy mismatch materially changes evolution without a fine-mesh run.

Saved corner sampling now confirms that the mismatch persists across all three
existing grids, with approximately 99% of upper-layer squared mismatch near the
outer edge. No extra fine simulation was performed. A routine-mesh energy-off
configuration and `compare --study energy` are ready to assess whether the whole
energy term materially affects the 100 s trajectory while retaining initial shape.


## Energy ablation result (2026-10-20)

The completed routine-mesh 100 s ablation preserved the initial mechanical shape
and set energy forcing exactly to zero. On/off final surface difference is
0.01595 nm RMS, 0.05989 nm at the edge, and 0.00089 nm at the center. The edge
effect is about 0.83% of total energy-off descent. This supports treating frozen
energy forcing as a modest short-time correction in this case, while retaining
the term. It does not validate the energy field or bound long-time effects.
No additional fine-mesh run is needed. Prioritize reference-temperature provenance
and physical validation of shared phase-field kinetics/wall inputs over tuning
corner energy to improve this short-time surface agreement.


## Two-way coupling development

The first prerequisite is implemented: fixed-mesh mechanical solves can continue
from copied plastic history and dimensionless displacement. Small-mesh tests pass
for split heating/cooling, same-load re-equilibration, and failed-call rollback.
The active phase runner still uses static one-way coupling. Evolving-geometry
history transport, coordinate consistency and transactional feedback remain to be
implemented; see `new/TWO_WAY_COUPLING.md`. No production or fine-mesh run was used.


The phase/reference mapping prerequisite is now implemented and verified on small
meshes. It follows the existing vertical-only initial-shape convention, preserves
subsequent phase changes, and rejects unsupported domain crossings. It is not a
conservative material transport map. The active simulation still uses one-way
coupling: diffuse occupancy versus sharp mechanical initialization and plastic
history treatment in newly occupied/vacated Cu remain the next implementation work.

Experimental mechanics feedback now has a working propose/commit API: a changed
binary Cu occupancy can trigger a mechanical solve with retained plastic history,
and failed/discarded trials leave accepted mechanics untouched. New Cu is rejected
by default; an explicit virgin-at-reference assumption exists only for experiments.
Small-mesh erosion/growth and rollback tests pass. The phase loop is still not
connected. Binary cell-centroid occupancy is nonconservative and can miss subcell
surface movement; reference-volume diagnostics expose this limitation. Do not
present this prototype as a validated two-way phase-field simulation.

## Experimental phase-loop connection

An optional programmatic feedback adapter now connects phase evolution, trial
mechanics, and updated energy. Small end-to-end rollback tests pass; a separate
synthetic growth test confirms energy is transferred in newly occupied Cu.
The tiny evolving-phase test itself changes no binary cells, so it is a control-
flow test, not physical validation. Production `run_model.py` remains one-way.
The unresolved binary/nonconservative transfer and new-Cu history assumptions
must be addressed before asking for meaningful long feedback simulations.
No user simulation is needed at this point, and no fine-mesh run is proposed.

## Conservative geometry transfer (2026-10-27)

An axisymmetric clipped-cell volume-fraction map is implemented and tested. It
preserves the supplied reference graph to roundoff and detects subcell motion
without mesh refinement. Experimental trials now carry the fraction map and
record binary-versus-graph volume error, but their constitutive solves still use
binary occupancy. Fractional plasticity/history rules remain to be implemented.
The saved-run audit also finds ~0.00927% reconstructed-graph volume drift despite
conserved diffuse phase mass; this is reported separately rather than corrected
by an artificial height shift. No new user simulation is needed yet.


### Fractional-cell mechanics update

The experimental two-way Python API now supports `occupancy_mode="fractional"`.
Small surface changes affect equilibrium without binary cell jumps. Cu uses its
intrinsic plasticity law; fraction weights stress and elastic energy once.
Partial-cell history uses an explicit homogeneous virgin-Cu mixing approximation;
growth is still rejected unless that policy is selected. Small mechanics and
four-step coupled rollback/retry checks pass. Production CLI remains static.
See `new/TWO_WAY_COUPLING.md` for equations, limitations, and verification commands.
Next: reconcile conserved diffuse phase mass with the mechanics surface-volume
map. No user production simulation or finer mesh is needed yet.

### Conservative phase-mass change transfer

Implemented opt-in `MechanicalFeedback(..., transfer_mode="column_mass")` for
fractional mechanics. Radial weighted projection converts changes in diffuse
phase mass to a mass-equivalent height graph, anchored to initial sharp geometry.
Exact integration and a pre-commit mechanics/phase balance check preserve mass
changes without a global surface offset. Analytic, independent FE integration,
and four-step coupled rollback/retry tests pass; maximum mechanics-volume drift
in that small test is 3.4e-16. Existing binary feedback and static weak-form
regressions pass.

The mass-equivalent surface is not the phi=0.5 contour: their difference reaches
0.882 nm on the tiny workflow mesh. Next is practical-mesh comparison and
transfer/update-interval sensitivity, not a new fine-mesh production run.
Constitutive/history mixture assumptions remain experimental. Production CLI
remains static and no user run is needed yet.

### Routine-mesh saved-field mass-transfer audit

Reused `runs/native_dt05_100s/phi.h5` on the existing 100x220 phase mesh; no new
simulation. Mass-equivalent versus contour height-change difference is 0.06356 nm
RMS and 0.18848 nm maximum, largest near the edge. At 101, 201, and 401 radial
projection nodes the profiles agree within 7e-10 nm. Keep existing 201-node
sampling; extra projection refinement is unnecessary for this saved field.
Mass changes agree with independent P1 integration to roundoff, and all reference
heights remain inside the mechanical domain. See `runs/column_mass_audit/ANALYSIS.md`
and its comparison plot. This verifies transfer resolution, not physical accuracy
or two-way evolution. Only endpoint fields exist: next is a short controlled
routine-mesh feedback/update-interval pilot, not a long or fine-mesh simulation.

### User-run two-way pilot ready

Added `run_feedback_pilot.py`: a separate experimental runner comparing mechanics
updates every 0.25 versus 0.5 s during two independent 2 s cases. Both use phase
100x220, mechanics80x80 per layer, identical initial native mechanics, fractional
occupancy, conservative column-mass transfer, and explicit virgin-Cu history
mixing. Phase dt is fixed at0.25 s to isolate interval sensitivity. Failure stops
and packages diagnostics rather than changing timestep. No production defaults
changed. The paired smoke test and archive-content checks passed; overwrite
protection was exercised. Routine pilot has not been run by the assistant.
Next action: user runs the command in RUNNING_THE_MODEL.md and provides the one
analysis.zip for interval-sensitivity review before any longer simulation.

### Routine two-way pilot analyzed (2 s)

`runs/feedback_pilot_20261028_170827` completed both cases without rejection.
Initial surfaces and energy arrays match exactly. Mechanics updates every 0.25 s
versus0.5 s yield surface differences 2.636e-5 nm RMS and 1.901e-4 nm maximum;
sampled energy relative L2 difference is 0.004845%. Phase mass and mechanics volume
remain conserved to roundoff; final intrinsic yield residual is below 5.1e-7 Pa.
Total runtime 82.76 s on the existing 100x220 phase/80x80 mechanics setup.
Mass-equivalent versus contour geometry discrepancy reaches 0.24056 nm, larger
than update-interval differences. This remains an experimental representation
assumption. Review saved in the run's ANALYSIS.md and repackaged archive.
Next user run: same paired study extended to 10 s using `--duration-s 10`. No mesh refinement or long run requested.

### Ten-second interval pilot reviewed; feedback-effect study ready

`runs/feedback_10s_20261028_175051` completed 40 steps per case without rejection.
Surface interval differences: 1.497e-5 nm RMS, 9.031e-5 nm maximum; energy relative
L2 difference 0.003606%. Mass is conserved to roundoff; final yield residuals
below 5.4e-7 Pa. Total runtime 213 s; 0.5 s mechanics updates save about 34% of
case runtime versus 0.25 s. Use 0.5 s provisionally for the next short check.

Added `--study feedback_strength`: updated versus frozen energy with identical
initial mapped energy, initial geometry, meshes, and phase timesteps. Both the
new smoke comparison and existing interval smoke regression pass. Next user run
is the 10 s feedback-effect command in RUNNING_THE_MODEL.md; it isolates whether
feedback materially changes surface evolution. No mesh refinement requested.

### Ten-second feedback effect measured

Analyzed `runs/feedback_effect_10s_20261029_194735`. Same initial mapped energy
and geometry confirmed. Updated versus frozen surface difference: 0.001763 nm
RMS, 0.010560 nm maximum at the edge, approximately 0.27% of frozen edge movement.
Sampled energy differs by 1.966%; center surface difference is negligible.
Both cases completed without rejection and conserved phase mass to roundoff.
The effect exceeds measured update-interval sensitivity by over 100x, but remains
small and subject to unvalidated surface/history assumptions. The geometry
representation discrepancy reaches 0.24013 nm; it is not a quantified error bound.
Saved review.json, ANALYSIS.md and feedback_review.png in the run directory.
Next user run: same feedback-strength pair for 100 s, same practical mesh, dt0.25 s
and 0.5 s mechanics updates. This tests effect accumulation over the established
100 s horizon, not a new physical feature or fine-mesh study.
