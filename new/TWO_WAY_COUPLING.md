# Two-way coupling implementation

## Implemented: persistent mechanics history on a fixed mesh

`thermoelastic_solver.solve` accepts `initial_plastic_state` and
`initial_displacement` from a prior result, plus explicit absolute `load_factors`.
Use `previous.u_dimless_field`, not the displacement field in metres.

```python
heated = solve(
    None, mesh, cu, 523.15, g.L0_m, .5, .993,
    geometry=g, material_dielectric=dielectric,
    plasticity_cu=law, load_factors=[.25, .5, .75, 1.],
)
reequilibrated = solve(
    None, mesh, cu, 523.15, g.L0_m, .5, .993,
    geometry=g, material_dielectric=dielectric,
    plasticity_cu=law, load_factors=[1.],
    initial_plastic_state=heated.plastic_state,
    initial_displacement=heated.u_dimless_field,
)
```

Factors are absolute fractions of the temperature difference from the material
reference, not fractions of the remaining temperature rise. A factor of 1 at the
same target temperature re-equilibrates at that temperature; it does not heat
again from room temperature. Explicit factors are required for continuation so
the default virgin heating ramp cannot accidentally be replayed.

Continuation copies the integration-point plastic tensor, accumulated equivalent
plastic strain, and displacement. The supplied state is never mutated, even if a
later increment fails. Each Newton iteration uses fixed committed history; both
plastic fields are evaluated before either is committed. Callers accept a returned
result only after the whole requested path succeeds. Different meshes or
incompatible field spaces are rejected. Keep geometry, material laws, normalization,
and constraints unchanged across this initial API; it does not validate or remap
history for those changes. This is in-memory continuation, not file restart.

Verification:

```sh
python new/validate_history_continuation.py
```

On an 8 × 8-per-layer mechanics mesh, tests check split versus continuous heating,
same-load re-equilibration, split cooling, input-state immutability, injected failure
after an accepted increment, and rejection of a different mesh. This is a small
implementation test, not a production calibration or a refined-mesh study.

## Not yet enabled: evolving-geometry feedback

The active terminal model remains static one-way coupling. Passing phi together
with plasticity still raises an explicit unsupported-operation error. Do not
remove that guard just to put a mechanics solve inside the phase loop.

The remaining work is:

1. Define the map between phase geometry and reference mechanical material
   occupancy. The phase surface already contains thermal displacement, so treating
   its raw height as a new stress-free geometry would count deformation twice.
2. Define plastic-history behavior for newly occupied and vacated material.
   A quadrature state tied to a fixed position is not automatically a transported
   material history. Remeshing, extending the mechanics domain, or changing the
   Cu/void constitutive interpolation also requires an explicit policy. The
   current mechanics top is 1000 nm while the phase domain reaches 1050 nm.
3. Implement a mechanical update using accepted phase states and preserved history,
   with transaction rollback of phase, displacement, plastic state and energy if
   either subsolve fails. Trial coupling iterations must not irreversibly commit
   plastic strain.
4. Transfer the updated elastic energy back without adding displacement a second
   time to the evolved surface. Test no-change and zero-feedback limits, then
   mechanical-update interval sensitivity on small/routine meshes.

No claim of a completed two-way coupled simulation is made at this stage. No
400 × 880 or finer uniform simulation is part of this work. The user runs
production simulations; implementation tests remain small.

## Implemented: phase/reference coordinate map

`phase_geometry_map.SurfaceReferenceMap` stores a copy of the initial Cu vertical
displacement profile. The active model initializes the phase surface at reference
radii and adds only vertical displacement; the new map follows that convention:

```
phase_z_query = reference_z - base_z + initial_uz(r)
reference_surface = phase_surface + base_z - initial_uz(r)
```

All quantities are normalized by L0. The original thermally displaced analytic
surface maps back to the reference Cu top. Subsequent phase-height changes are
preserved. This frozen column shift is an explicit approximation compatible with
the existing initializer, not a full inverse deformation map: it does not include
radial deformation, conserve transported mass, or advect plastic strain.

`pull_back_phase_pair` samples both initial and current scalar phase fields on a
mechanics scalar space and returns those fields plus their difference. Identical
input fields give exactly zero change. Dielectric/below-pad target nodes receive
zero Cu phase contribution; dielectric identity remains a separate mechanics mask.
Small coordinate roundoff at the column boundary is tolerated. Missing phase
queries and surfaces outside the mechanical domain are rejected, not clipped.
This transfer currently supports one MPI rank.

Run `python new/validate_phase_geometry_map.py` for inexpensive map checks. Tests
cover thermal-shift removal, retained surface motion, immutable initial profile,
exact affine FE transfer, zero-change behavior and out-of-domain rejection.

The initial mapped phase is diffuse, whereas the calibrated plastic mechanics
uses a sharp material map. They must not be substituted silently. Neither the
transferred difference nor the full mapped phase is currently fed into plastic
mechanics. The next implementation must choose and test consistent occupancy and
history treatment before enabling the existing evolving-plasticity guard.

## Implemented: experimental occupancy/history mechanics updates

`evolving_mechanics.EvolvingMechanics` now supports proposing and committing a
mechanical response to a changed reference surface on the existing mesh.
Initialize it from the original sharp-geometry thermal/plastic result. Geometry,
materials, temperature, normalization, and constraints must match that result.

- A binary DG0 Cu occupancy is sampled at cell centroids. Fixed dielectric
  identity is retained. No diffuse-mixture plastic law is introduced.
- Retained Cu keeps its integration-point plastic tensor and equivalent strain.
- Removed Cu has those histories cleared in the trial copy.
- New Cu is rejected with the default `new_cu_history='reject_growth'`.
- Explicit `new_cu_history='virgin_at_reference'` initializes new Cu's plastic
  history to zero, then equilibrates it with the existing thermal strain relative
  to the original reference temperature. This is **not** deposition stress-free
  at the current anneal temperature and is not a validated diffusion-history law.
- `propose(radii, phase_surface)` maps the surface back, prepares a copied history,
  and actually solves mechanical equilibrium at load factor 1. It returns updated
  displacement, stress, elastic energy, occupancy and transition diagnostics.
- `commit(trial)` accepts a successful proposal. Discarding it leaves accepted
  mechanics unchanged. Foreign/stale proposals are rejected. No trial history is
  committed during a mechanical Newton iteration.

The diagnostic reference Cu volumes omit the common 2π factor. Binary centroid
sampling is not conservative: it can change material volume in discrete jumps,
and small surface movements can produce no changed cells. The volume change is
reported, not corrected. The synthetic erosion/growth tests deliberately use
movements large enough to resolve on the tiny test mesh; they are not a physical
surface-diffusion simulation or a demonstration of mass-conserving coupling.

The low-level `solve(..., cu_occupancy=...)` entry point requires binary DG0 on the
same mesh and rejects simultaneous phi input. Correct history preparation is the
caller's responsibility; prefer the transactional wrapper. Dielectric plasticity
is intentionally unsupported for these occupancy transitions.

Verification:

```sh
python new/validate_evolving_mechanics.py
```

The small-mesh tests confirm unchanged geometry reproduces the original sharp
solve, erosion changes equilibrium and W, accepted histories are not mutated by
proposals or failure, stale commits fail, and explicit growth reaches a finite
mechanical solution. Existing history-continuation checks also pass.

## Remaining before a coupled phase-loop prototype

The active `run_model.py` remains one-way. The next integration must transfer the
trial energy using the **new** occupancy (the static sampler's original Cu-top
cutoff is inappropriate), coordinate phase/energy/mechanics rollback, and monitor
reference occupancy volume error. A production-quality scheme additionally needs
an accepted new-material history model and a conservative, sufficiently resolved
phase-to-material transfer. The frozen-column map and binary occupancy are
explicit prototype approximations, not settled coupling physics.

## Implemented: experimental staggered phase-loop integration

`run_phase_field(..., feedback=MechanicalFeedback(engine, ...))` now connects the
phase solver and transactional mechanics engine. Omitting `feedback` preserves
the existing production one-way path; no terminal default or production config
has been changed.

For each selected accepted phase step:

1. Solve phase evolution using the previous committed energy.
2. Extract the tentative surface. Subtract its sampled initial profile and add
   the analytic initial profile, so contour discretization error is not mistaken
   for a physical geometry change on the first update.
3. Propose mechanics at the corrected surface, retaining history according to
   the explicit occupancy policy.
4. Sample the new Cu-localized energy through the inverse frozen-column map.
   The old static Cu-top cutoff is not applied; newly occupied cells can supply
   energy. Initial experimental W uses this same map for consistency.
5. After all trial calculations succeed, commit mechanics, W and the phase
   state. The next phase step uses updated W. There are no inner coupling
   iterations and no second addition of mechanical displacement to phi.

A failed mechanics solve rejects the tentative phase step and reduces dt; its
uncommitted mechanics/energy trial is discarded. Invalid geometry/history policy
or excessive reference-volume change stops with the last accepted state intact.
The configurable cumulative reference-volume limit defaults to 1%. It is a
prototype rejection threshold, **not** a claim of acceptable physical mass error.
Phase-field conservation and binary mechanical occupancy volume are distinct.

`coupling_history.json` reports update times, changed cells, history transitions
and reference-volume changes. `latest_energy_samples.csv` stores the latest
committed phase-grid W; `coupling_samples.csv` stores the initial W. The existing
`static_energy.xdmf` file still contains only the initial field. For this prototype,
reported height relative to the dielectric uses its initial deformed edge,
explicitly labeled in the summary; it is not a dynamically updated gap metric.

The experimental initialization changes the energy-coordinate mapping relative
to the old static sampler. A direct static/prototype comparison therefore does
not isolate feedback alone. A matched map with updates disabled is needed for
that experiment.

Verification:

```sh
python new/validate_feedback_loop.py
```

Four tiny real coupled steps pass, with a retry after an injected post-mechanics
failure matching an uninterrupted smaller-step run. Input mechanics remains
unchanged. Separate synthetic growth verifies nonzero energy in newly occupied
regions, and a strict volume limit rejects the trial without commit. The tiny
phase trajectory itself changes no binary cells, so that test verifies control
flow and rollback, not a resolved geometry-feedback effect. The coarse test phi
has overshoot; it is not an accuracy benchmark.

## Current remaining work

No user production run is requested yet. The key remaining limitation is the
binary/nonconservative geometry transfer: it can miss small surface motion and
create abrupt volume changes. Address that limitation and the history policy
before interpreting a longer experimental feedback trajectory physically. The
production terminal workflow remains static. No fine uniform mesh is needed to
continue this implementation work.

## Implemented: conservative axisymmetric graph-to-cell fractions

`conservative_occupancy.py` clips each affine triangular mechanics cell against
the piecewise-linear reference surface. The DG0 fraction is
`integral_Cu(r dA) / integral_cell(r dA)`, with the common 2π omitted. Planar area
fractions are not substituted for axisymmetric volume fractions.

The fractions preserve the supplied graph's volume to roundoff and respond to
subcell movements that leave centroid binary occupancy unchanged. They reproduce
the original sharp baseline on its conforming mesh. Inputs outside the supported
domain are rejected, and a global coverage check catches incomplete mesh coverage.
The implementation currently requires serial affine triangles.

Each experimental mechanical proposal now computes this fraction map and records
its geometry diagnostics plus the binary-map volume error. The trial exposes
`reference_fraction` for later use. The default remains binary occupancy. The opt-in fractional constitutive
branch described below now uses these fractions in the mechanics solve.

Verification and inexpensive saved-result audit:

```sh
python new/validate_conservative_occupancy.py
python new/audit_geometry_transfer.py
```

Geometry tests include an analytic axisymmetric cut triangle, orientation
invariance, known subcell volume changes, a zero-volume shape redistribution,
multiple surface segments within one cell, and independent FE radial integration.
No new production simulation or finer mesh is used by the audit.

For `runs/native_dt05_100s`, the final graph map has relative integration error
about 1.1e-16. There are 86 partially occupied cells; fractions change in 97 cells,
whereas binary occupancy changes in only 26 cells. Binary final volume differs
from the reconstructed graph by -0.009157%.

A separate issue remains: the reference surface graph itself changes volume by
-0.009270%, whereas the saved diffuse phase-field mass drift is about -3.7e-15.
Thus conserving the extracted graph is **not** yet conserving transferred phase
mass. The audit deliberately reports both. No global surface offset, arbitrary
rescaling, or correction to the stored simulation has been applied.

Next: reconcile (or quantify) diffuse-phase mass versus reconstructed graph
volume before claiming conservative two-way coupling. The terminal model remains unchanged.


## Experimental fractional mechanics (implemented)

`EvolvingMechanics(..., occupancy_mode="fractional")` opts into the conservative
reference-graph fraction map. Default `binary` behavior and the static terminal
runner are preserved. This is a uniform-strain cell mixture approximation; it
is not an interface-resolved traction-free cut-cell formulation.

Cu return mapping uses intrinsic Cu elastic constants, thermal strain, yield
stress and hardening. Only the resulting stress and stored elastic energy are
weighted by Cu fraction. Elastic dielectric and residual void contributions are
calculated separately. Cu energy is weighted once, not squared by a second mask.
Exported partial-cell stress is mixture stress; plastic strain is intrinsic Cu
history. No plastic dissipation is included in the stored elastic-energy field.

History policy is explicit: partial loss retains the remaining Cu's intrinsic
history, complete loss clears it, and growth is rejected by default. With
`new_cu_history="virgin_at_reference"`, newly added Cu is mixed with retained Cu:
`new_history = old_history * min(old_fraction, new_fraction) / new_fraction`.
This applies to plastic strain and accumulated equivalent plastic strain, with
zero for empty cells. Dilution can decrease intrinsic accumulated strain at a
geometry transition; subsequent mechanical return mapping remains monotone.
This homogeneous mixing assumption does not transport plastic-history
populations or guarantee dissipation/energy conservation during shape changes.

Small verification commands (inside `fenicsx-env`):

```sh
python new/validate_fractional_mechanics.py
python new/validate_feedback_loop.py --occupancy-mode fractional
```

Verified: binary limit reproduces the original solution; 0.005 nm growth and
loss change equilibrium and energy on an 8-by-8-per-layer mesh; Cu energy is
weighted exactly once; yield consistency, history dilution/removal, caller-state
immutability, and four real phase/mechanics updates pass. The end-to-end test also
checks rollback after an injected failure and agreement with the uninterrupted
smaller-step run. These are implementation checks, not physical validation.
The tiny coupled test changes 13 fraction cells per update and accumulates about
-0.0232% reference-graph volume change; the volume guard accepts this within its
1% test limit. Graph-volume drift remains an outstanding transfer issue.

No new production or fine-mesh run is requested. Next work is a mass-consistent
phase-to-mechanics transfer, followed by a practical-mesh sensitivity study before
exposing two-way coupling through the production CLI.

## Mass-equivalent column transfer (implemented, opt-in)

The experimental adapter now accepts `transfer_mode="column_mass"` together with
`EvolvingMechanics(..., occupancy_mode="fractional")`. The default remains
`transfer_mode="contour"`; the production terminal runner still has no feedback.

`column_mass_transfer.py` projects the full P1 phase change relative to the saved
initial phase onto a radial piecewise-linear height change. For every radial
basis function N_i it solves

    integral N_i(r) * dh(r) * r dr
      = integral N_i(r) * (phi(r,z) - phi_initial(r,z)) * r dr dz.

Because the radial basis sums to one, the resulting graph-volume change equals
the diffuse phase-mass change. The reference mechanics height is its initial
sharp height plus dh. This anchors the original geometry instead of interpreting
initial diffuse-interface mass as a different Cu volume. It matches the existing
frozen vertical-column coordinate convention, whose axisymmetric integration
weight is unchanged by that shift; it is not a finite-deformation mass balance.

The transfer matrix is assembled once by clipping phase triangles at radial
knots and integrating the cubic polynomial exactly. It is reused at each update.
There is no global height adjustment, phase clipping, or fine-mesh requirement.
Fractional cell integration then preserves the graph volume. Each proposal
checks the final mechanics-volume change against the phase-mass change before
commit. Out-of-domain graphs are rejected rather than clipped.

This is a **mass-equivalent graph**, not the phi=0.5 contour. Vertical diffuse
redistribution, bulk phase deviations, and numerical overshoot all contribute
to the column integral. Local topology, overhangs, and separate voids cannot be
represented by a single graph. Radial L2 projection can also overshoot locally;
only the existing domain bounds constrain it. Fractional-cell stress and plastic
history still use the documented mixture approximation. Mass conservation alone
does not validate those physical assumptions.

Verification commands in `fenicsx-env`:

```sh
python new/validate_column_mass_transfer.py
python new/validate_feedback_loop.py --occupancy-mode fractional --transfer-mode column_mass
```

Analytic linear radial redistribution, zero change, nonaligned triangle/radial
knots, and independent FEniCS mass integration pass. The fractional mechanics map
also preserves a prescribed nonzero phase-mass change. The four-step coupled test
passes rollback/retry checks and updates 17 fraction cells at its final step.
Its maximum cumulative relative mechanics-volume drift is 3.4e-16, compared with
about 2.32e-4 for the contour-based fractional test. These are tiny workflow
checks, not converged simulations. The maximum mass-equivalent versus contour
height difference reaches about 0.882 nm on this deliberately coarse test; it
must be assessed on the routine mesh before production use.

Next: compare the two transferred surfaces and mechanics response on the routine
mesh, including sensitivity to radial projection resolution and update interval.
Keep phase meshes at the established practical sizes; no new fine-mesh run is
requested. Do not enable production two-way coupling solely on these tests.

## Routine-mesh saved-field audit

`audit_column_mass_transfer.py` reuses XDMF/HDF5 nodal snapshots with the same
projection implementation used by the live adapter. It requires numpy, scipy,
h5py and matplotlib, but no FEniCS. In this workspace the base Anaconda Python
has those dependencies:

```sh
/opt/anaconda3/bin/python new/audit_column_mass_transfer.py \
  --run runs/native_dt05_100s --output runs/column_mass_audit_repeat
```

The output directory must be empty. The script verifies the saved final time,
checks phase mass independently with exact triangle integration, and writes JSON,
CSV, and a comparison plot. It performs no phase/mechanics solves.

On the saved 100x220, 100 s trajectory, the contour and mass-equivalent changes
differ by 0.06356 nm RMS (uniform radial weighting), at most 0.18848 nm. The
101/201/401-node radial projections differ by less than 7e-10 nm for this field.
Retain 201 nodes; this is projection-resolution evidence, not a convergence claim
for the phase mesh or fractional constitutive law. All reference heights pass
the domain bounds. Mass is preserved to roundoff versus independent integration.
The routine-mesh difference is smaller than the 0.882 nm tiny workflow-test
result; these different trajectories are not a controlled mesh convergence pair.

Only two phase snapshots are saved. They cannot establish feedback interval or
plastic-path sensitivity. The next study needs a short two-way pilot on the
routine mesh with controlled update intervals, before a long-duration run.

## Controlled pilot entry point

`run_feedback_pilot.py` provides a separate experimental entry point for the next
study. It holds phase dt fixed, starts both interval cases from the same initial
mechanical state, checks accepted-step/update counts and final mechanical time,
and archives both histories and the source. The 2 s routine preset uses the
existing 100x220 phase and80x80 mechanics meshes, with update intervals1 and2
steps atdt0.25 s. A failure halts the study; there is no automatic timestep change
or inference of convergence from two intervals. `--smoke` verifies the workflow
on the small mesh. See RUNNING_THE_MODEL.md for the user command and artifacts.
