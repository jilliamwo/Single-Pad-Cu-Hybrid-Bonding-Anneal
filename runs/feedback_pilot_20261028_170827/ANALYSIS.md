# Two-second feedback pilot review

Both cases completed eight fixed 0.25 s phase steps, with no rejected steps or
solver rebuilds. Case 1 performed eight mechanics updates; case 2 performed four.
Initial sampled surfaces and initial energy arrays are exactly identical.
The phase mesh is 100x220 and mechanics is 80x80 per layer; no fine mesh was used.

| Metric | Every0.25 s versus every 0.5 s |
|---|---:|
| Final surface RMS difference | 0.00002636 nm |
| Maximum surface difference (at edge) | 0.00019009 nm |
| Sampled energy relative L2 difference | 0.004845% |
| Maximum equivalent-plastic-strain difference | 3.603e-6 |

Both surfaces evolve materially despite the small interval difference: case 1
center drops 0.14493 nm and edge drops 2.50921 nm from the initial contour.
The final center remains 2.26634 nm below the initial deformed dielectric edge.
These are contour heights; the mechanics mass-equivalent graph is a different
surface representation, whose maximum discrepancy from the contour during this
pilot reaches 0.24056 nm. That discrepancy is much larger than interval sensitivity
and remains a modeling issue to examine rather than an error removed by faster
mechanics updates.

Relative phase-mass drift is 5.12e-15 and 6.90e-15. Maximum relative mechanics-volume
drift is 7.83e-16 and 2.24e-16. Final intrinsic yield violations are below 5.1e-7 Pa,
consistent with roundoff. These support numerical consistency of this short test,
not validation of the material/history mixture or physical two-way model.

Total wall time is 82.76 s: case 1 takes 29.52 s and case 2 takes 22.90 s, with shared
initialization accounting for the remainder. The less frequent update saves
about 22% of case runtime here, but long-run timing and accuracy remain untested.

Next: extend the same paired interval study to 10 s with the same mesh and phase
timestep. Retain both update intervals for this longer comparison; do not infer
long-time convergence or choose a universal update interval from the 2 s test.
The stress-free 22 C reference is still unverified against the ANSYS setup.
