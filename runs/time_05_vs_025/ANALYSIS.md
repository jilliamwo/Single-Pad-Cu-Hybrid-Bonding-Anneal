# Time-step sensitivity analysis at 100 s

All three runs completed, at actual fixed steps of 1, 0.5, and 0.25 s, with
100, 200, and 400 accepted steps. No rejected steps or solver rebuilds occurred.
Mass drift was below 1e-14 and all diagnostic surface crossings were found.
Initial profiles agreed within 4e-12 nm, and transferred energies agreed within
floating-point solve variation. Thus this study isolates time discretization.

| Refinement | Final surface RMS difference (nm) | Maximum difference (nm) |
|---|---:|---:|
| 1 to 0.5 s | 0.002099 | 0.004062 |
| 0.5 to 0.25 s | 0.001060 | 0.002058 |

The RMS reduction factor is 1.9809, corresponding to an empirical order
log2(ratio) = 0.9861. Successive difference vectors have alignment
0.999963. This is evidence of
approximately first-order temporal convergence for the extracted final surface
on this mesh and at this duration, not a universal accuracy guarantee.

The center relative to the deformed dielectric edge approaches -2.310374 nm
at dt=0.25 s, versus -2.313861 nm at dt=1 s. The last refinement changes the
center by 0.001176 nm and the edge by 0.002030 nm. Small phi overshoot persists
(maximum approximately 1.000218); reducing dt does not eliminate it.

The full 1-to-0.25 s surface RMS change is 0.003159 nm,
well below the earlier 0.052040 nm native/ANSYS-input RMS discrepancy. The native
trajectory is weakly sensitive to dt at this scale; an ANSYS-input time study
has not been performed, so this does not rigorously bound both trajectories.

Proceed with the prepared mesh study at dt=0.5 s, using the existing dt=0.5 s run
as reference. Its remaining time sensitivity is small enough to make this a useful
next diagnostic. If mesh differences are themselves around 0.001-0.002 nm, refine
time on the finer mesh before claiming spatial convergence. Fixed 5 nm interface
width and unchanged mechanics isolate phase discretization plus energy-transfer
representation. The 22 C reference and near-surface energy mismatch remain open.
