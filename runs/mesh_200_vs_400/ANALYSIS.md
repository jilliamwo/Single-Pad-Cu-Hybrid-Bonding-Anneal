# Three-mesh analysis at matched dt=0.25 s

| Refinement | Final surface RMS change | Maximum change |
|---|---:|---:|
| 2.5 to 1.25 nm | 0.065438 nm | 0.231366 nm |
| 1.25 to 0.625 nm | 0.016414 nm | 0.059390 nm |

RMS differences shrink by 3.9867, giving an empirical spatial order
1.9952. This is consistent with approximately second-order convergence of
the extracted surface in this study. It is not proof of convergence of all
fields, of the mechanics solution, or of the continuum physical model.
The study combines phase discretization and transfer representation changes.

The 400 × 880 run completed 400 steps at dt=0.25 s without rejection or solver
rebuild. Maximum Newton iterations: 3. Final mass drift was
3.55e-14; all surface crossings were found, and nodal phi extrema were within
[0,1]. Its center rises by 0.093291 nm over 100 s, compared with 0.033132 nm on
the 200 × 440 mesh. Final center height is 2.019617 nm below the deformed
dielectric edge. The finer two meshes agree on the upward direction of net
center movement, but its magnitude remains resolution-sensitive.

The 200 × 440 mesh differs from this benchmark by 0.016414 nm RMS and up to
0.059390 nm at the center. This is a measured pairwise difference, not a strict
error bound. It is a reasonable working compromise for selected quantitative
comparisons, while 100 × 220 remains appropriate for inexpensive exploratory
runs with explicit spatial sensitivity. Do not equate the coarse mesh's close
native/ANSYS agreement with converged physical accuracy.

## User runtime constraint

The user reports this 400 × 880 run took multiple hours and explicitly asks to
avoid meshes this fine in future. Preserve it as an existing benchmark; do not
request or launch another 400 × 880 or finer uniform run. Use 100 × 220 for
routine screening and reserve 200 × 440 for selected checks. This is a ceiling,
not a promise of short runtime. Reuse existing results and favor postprocessing
or short controlled tests before asking for another simulation. Any future
higher-resolution approach requires a separate user decision; do not infer
authorization from this completed benchmark.

Next priority is diagnosing the existing near-surface energy-transfer mismatch
from saved fields. No new expensive simulation is needed for that initial
investigation. A later matched native/ANSYS check can use the 200 × 440 mesh,
if warranted, rather than rerunning this benchmark.
