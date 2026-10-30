# Ten-second feedback-effect review

Both cases completed 40 fixed 0.25 s phase steps without rejected steps, rebuilds,
or missing contour points. Initial surface and sampled energy arrays were
independently checked for exact equality. Updated mechanics ran every 0.5 s;
the control retained the identical initial mapped energy throughout.

| Updated minus frozen at 10 s | Result |
|---|---:|
| Contour RMS difference | 0.001763 nm |
| Maximum difference, at r=250 nm | 0.010560 nm |
| Center difference | -0.000000957 nm |
| Sampled energy relative L2 difference | 1.9656% |
| Maximum equivalent plastic strain difference | 0.004599 |

The updated contour edge moves -3.88200 nm, versus -3.89256 nm for the frozen
case: feedback slightly reduces edge recession by about 0.27% of the frozen
edge movement. The effect develops gradually during the saved trajectory and
is concentrated near the pad edge. The center effect is negligible here.
The updated final surface reproduces the previous interval-2 run within
3.5e-12 nm, supporting reproducibility.

The 10 s interval study gave 0.00001497 nm RMS and 0.00009031 nm maximum
surface differences. The feedback effect exceeds those by factors of about
118 and 117 respectively. This separates the feedback signal from measured
update-interval sensitivity; it does not establish physical accuracy, timestep
or mesh independence, or an acceptable engineering tolerance.

Maximum phase-mass drift over both trajectories is about 1.22e-14. Updated
mechanics-volume drift stays below 5.60e-16; final intrinsic yield violation
is below 3.8e-7 Pa. The frozen case has no geometry updates, so its reported
zero mechanics-volume drift is not an independent phase-to-geometry check.
The accumulated-plastic-strain maximum increases from 0.03912 initially to
0.04372 in the updated case; pointwise changes include the explicit history
mixing approximation and are not a measure of plastic dissipation.

The maximum mass-equivalent versus contour discrepancy is 0.24013 nm, much
larger than the 0.01056 nm feedback signal. This is a distinction between surface
representations, not a quantified error bound, but it limits physical claims.
The unverified 22 C stress-free reference, fixed column coordinates, and
homogeneous new-Cu history mixing remain unresolved assumptions.

Total runtime was 141.81 s, including shared initial mechanics. Updated and
frozen phase cases took 76.14 and 31.57 s, respectively. Next: extend this same
comparison to 100 s on the existing practical mesh, aligning with the previously
studied time horizon. Keep dt=0.25 s and mechanics updates every 0.5 s. Do not
refine the mesh for this check. Rough linear runtime scaling suggests about
18-20 minutes for the pair; nonlinear convergence can change that estimate.
Review conservation, geometry discrepancy, plastic evolution, and effect growth
before extending further or claiming the two-way model is physically validated.
