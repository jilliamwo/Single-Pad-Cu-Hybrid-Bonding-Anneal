# Corner-energy sampling across existing phase grids

Only saved arrays and the original ANSYS CSV were evaluated. No FEA or phase
simulation was run. One fixed ANSYS linear interpolant was evaluated at each
native run's saved phase nodes, using the archived source-mask convention.

| Phase grid | Upper 10 nm relative RMS mismatch | Outer 25 nm share of squared mismatch |
|---|---:|---:|
| 100 × 220 | 119.62% | 98.55% |
| 200 × 440 | 142.77% | 99.03% |
| 400 × 880 (existing benchmark only) | 141.47% | 99.02% |

At coordinates shared with the routine grid, native energy values agree to
better than 8e-7 J/m³ across these runs. The higher mismatch on finer sampling
therefore reflects the denser evaluation of the existing field differences,
not a changed mechanics problem. These are unweighted nodal statistics;
not quadrature-integrated errors, nor phase evolution with ANSYS on finer meshes.

This rules out the hypothesis that the near-surface discrepancy disappears
simply by sampling the existing mechanical fields more densely. It does not
separate differences in mechanical solutions from DG0/nodal postprocessing.
No basis has been established for changing the material law or rescaling energy.

Next, measure how much energy forcing matters to evolution: use
configs/model_energy_off_dt05.json on the 100 × 220 mesh for 100 s. Compare to
runs/native_dt05_100s. The mechanics calculation and initial deformed profile
remain active and identical; only mechanical_energy_multiplier changes from 1
to 0. Do not use --source none, which also changes the initial shape.

This ablation measures the effect of the whole mechanical-energy term, not the
error due specifically to the corner mismatch. Small influence would reduce its
priority at this duration; large influence would motivate a more targeted field
comparison. It cannot validate either energy field or predict long-time effects.
