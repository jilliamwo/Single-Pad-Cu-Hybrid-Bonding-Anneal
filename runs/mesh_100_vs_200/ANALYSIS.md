# Phase-mesh sensitivity at 100 s

Both runs completed 200 accepted steps at dt=0.5 s, with no rejected steps or
solver rebuilds and at most three Newton iterations per step. All surface
crossings were found. Mass drift magnitudes were below 4e-15. Mechanics resolution,
thermal load increments, material properties, and physical interface width were
held fixed.

| Metric | Coarse versus fine difference |
|---|---:|
| Initial surface RMS | 0.007665 nm |
| Final surface RMS (uniform radial samples) | 0.065416 nm |
| Final surface radial-area-weighted RMS | 0.049086 nm |
| RMS difference in change from initial surface | 0.061902 nm |
| Maximum final difference | 0.231351 nm |
| Final center, coarse minus fine | -0.231351 nm |
| Final edge, coarse minus fine | +0.212589 nm |

The coarse mesh (2.5 nm nominal spacing) predicts a center displacement during
phase evolution of -0.190143 nm. The fine mesh (1.25 nm) predicts +0.031940 nm.
Thus even the sign of net center motion over 100 s is unresolved. The edge moves
down 7.167025 nm versus 7.371618 nm. Differences are strongest near the axis and
outer edge, but the 25–225 nm radial interval still has 0.031178 nm RMS difference.
The difference is not merely a single poorly extracted center point.

Initial contour extraction relative to the analytic initial profile improves
from 0.006932 to 0.000921 nm RMS on refinement. However, the 0.061902 nm difference
in subsequent evolution remains after subtracting initial shapes. Initial-shape
sampling alone therefore does not explain the final discrepancy.

Final phi extrema improve from [-0.0000730, 1.0002167] to approximately
[0.000000395, 1.0]. This is encouraging, but bounded extrema and mass conservation
do not establish spatial accuracy. Two meshes cannot establish a convergence
order or prove the finer answer is the converged answer.

The 0.065416 nm mesh sensitivity exceeds the earlier 0.052040 nm native/ANSYS
surface RMS difference (both on the coarse mesh). Agreement between mechanical
input paths remains a useful consistency result, but it is not evidence of
spatially converged physics. Both paths can share discretization error.

Refinement also changes the P1 representation of frozen energy on the phase
mesh. The present results combine phase-discretization, initial-profile, and
energy-transfer representation effects; they do not isolate an error in the
constitutive law, axisymmetric equation, or wetting term.

Next: first check dt=0.25 s on the 200 × 440 phase mesh, since the previous time
study was only on the coarse mesh. A configuration is provided in
configs/model_mesh_fine_dt025.json. If fine-mesh time sensitivity remains small,
add a third spatial resolution at matched dt, and repeat the ANSYS-driven case
on the selected mesh. Defer quantitative claims about net center evolution and
long-run predictions until that sensitivity is understood.
