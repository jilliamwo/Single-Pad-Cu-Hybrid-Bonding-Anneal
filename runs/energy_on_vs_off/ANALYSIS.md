# Mechanical-energy ablation at 100 s

The comparison isolates mechanical-energy forcing on the routine 100 × 220 mesh
at dt=0.5 s. Mechanical surface profiles match within 3.81e-12 nm; the initial
phase surfaces match to the same precision. All saved energy-off W values are
exactly zero. Native thermal mechanics and Cu plasticity remain enabled.

Both runs completed 200 accepted steps without rejection or solver rebuild;
at most three Newton iterations were needed. All sampled surface crossings were
found. Energy-off mass drift remained below 8.3e-15 at recorded times. Small
nodal phi overshoot persists in both runs, as in earlier routine-mesh tests.

| At 100 s: energy on minus energy off | Difference |
|---|---:|
| Surface RMS | 0.015950 nm |
| Maximum absolute difference | 0.059891 nm |
| Center | -0.000890 nm |
| Outer edge | +0.059891 nm |

With energy on, the edge moves down 7.167025 nm; with it off, the edge moves down
7.226916 nm. The energy term therefore reduces edge descent by 0.059891 nm,
approximately 0.83% of the energy-off descent. About 69% of the squared surface
difference lies in the outermost 25 nm; excluding that region, RMS is 0.009359 nm.
Center motion is -0.190143 nm versus -0.189253 nm.

Interpretation: on this mesh and over 100 s, the frozen energy term makes a modest
correction, concentrated toward the edge. Most observed evolution remains when
that term is disabled. This does not mean mechanics is irrelevant: the deformed
initial shape is retained. Nor does it establish which mechanical-energy field is
correct or specifically isolate the corner's contribution.

The center effect (0.000890 nm) is smaller than the previously measured center
change from dt=0.5 to 0.25 s (~0.001176 nm); do not overinterpret its magnitude
without a matched ablation refinement. The full ablation RMS (0.015950 nm) is
smaller than the coarse-to-fine mesh sensitivity (~0.06544 nm). Correlated numerical
errors can cancel in on/off differences, so this comparison does not prove the
energy effect is unresolved, but its convergence has not been demonstrated.

The result supports lowering the priority of corner-energy tuning for short-time
routine-mesh screening. Keep the current energy term; do not remove or empirically
rescale it based on this test. The large local energy discrepancy alone does not
explain the earlier native/ANSYS center offset, which was mostly present initially.
Disabling all energy is not a rigorous upper bound on swapping two different fields.

No further fine-mesh run is proposed. Next prioritize confirming ANSYS's reference
temperature and assessing the shared phase-field parameters against physical
observables (surface profiles, edge shape, or evolution times). A longer routine-mesh
run can later be used as exploratory screening, not as a spatially converged
prediction. This ablation does not establish behavior at 7000 s or with feedback.
