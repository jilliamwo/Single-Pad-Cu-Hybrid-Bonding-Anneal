# Analysis of the paired 100 s runs

Both runs completed 100 steps at dt=1 s with no rejected steps or solver rebuilds.
Maximum Newton iterations per step: 3 in each run. Mass drift magnitude remained
below 1e-14. All sampled radial positions had an extracted surface. Final phi
overshoot above one was about 2.2e-4; undershoot below zero was less than 8.4e-5.
These are favorable diagnostics, not a mesh/time convergence demonstration.

| Metric | Value |
|---|---:|
| Initial surface RMSE, native versus ANSYS input | 0.052894 nm |
| Final surface RMSE | 0.052040 nm |
| Maximum final absolute surface difference | 0.067832 nm |
| RMS difference in surface change from each initial state | 0.015677 nm |
| Final center difference, native minus ANSYS | -0.067624 nm |
| Difference in center change during evolution | -0.004234 nm |
| Final edge difference, native minus ANSYS | +0.004889 nm |

The native center moved down 0.192454 nm; the ANSYS-driven center moved down
0.188220 nm. Thus most of their center offset predates phase evolution. The Cu
edge moved down 7.163018 nm and 7.050786 nm, respectively. Edge evolution is much
larger than center evolution in both runs; matching it does not independently
validate the shared wall energy or kinetics.

Relative to each run's own deformed dielectric edge, final center heights are
-2.313861 nm (native) and -2.128421 nm (ANSYS), a difference of -0.185440 nm.
This is not the same as the absolute surface difference: the native dielectric
edge displacement is 1.364816 nm versus 1.247000 nm in the imported ANSYS data.
The 0.117816 nm reference-plane difference adds to the 0.067624 nm Cu-center
height difference. These two metrics must not be conflated.

Transferred energy remains a separate concern. On common P1 sample nodes strictly
inside reference Cu, the unweighted relative RMS energy difference is 10.07%.
For nodes within the upper 10 nm of reference Cu it is 119.61%, with mean values
32810 J/m3 (native) and 24343 J/m3 (ANSYS). These are nodal statistics, not
axisymmetric volume-weighted integrals or an isolated constitutive-law error.
The upper layer includes corner effects and different transfer/interpolation
conventions. A small near-surface energy denominator also affects relative errors.
Do not infer energy-field equivalence from the close evolved shapes. Controlled
energy-off and transfer-sensitivity comparisons would be needed to isolate causes.

Recommended order: (1) halve the phase time step at fixed mesh; (2) refine the
phase mesh at fixed time step and interface width; (3) investigate near-surface
energy transfer if claiming quantitative stress-coupled predictions; (4) extend
the duration. Verify the assumed 22 C stress-free temperature against the ANSYS
setup before treating calibration as fully specified. The comparison validates
short-time agreement between two mechanical-input paths for this SiO2 case,
not transient ANSYS equivalence or feedback coupling.
