# Ten-second update-interval review

Both cases completed 40 fixed 0.25 s steps without rejections or solver rebuilds.
Initial surface and energy arrays are exactly identical. Mechanics updated 40
versus 20 times, using the existing 100x220 phase and 80x80 mechanics meshes.

| Difference: 0.25 s updates minus 0.5 s updates | Result |
|---|---:|
| Surface RMS | 0.00001497 nm |
| Maximum surface difference, at edge | 0.00009031 nm |
| Sampled energy relative L2 | 0.003606% |
| Maximum equivalent plastic strain difference | 1.935e-6 |

These endpoint differences are smaller than in the 2 s pilot; the two endpoints
alone do not establish a monotonic error trend or a convergence order. At 10 s,
case 1's center has moved -0.19071 nm and edge -3.88191 nm relative to its initial
contour, so the small interval difference is not due to an unchanged surface.
The final center is 2.31212 nm below the initial deformed dielectric edge.

Relative phase-mass drift is -8.46e-15 and -5.34e-15. Mechanics-volume drift stays
below 7.83e-16 and 5.60e-16. Final intrinsic yield violations are below 5.4e-7 Pa.
No extracted surface points are missing. The maximum mass-equivalent versus
contour discrepancy remains about 0.241 nm over the trajectory, much larger than
the update-interval effect and still a representation assumption to assess.

Total runtime: 212.95 s (3.55 minutes). Case runtimes: 108.93 s and 72.38 s.
The 0.5 s update interval saves about 33.5% of case runtime here. This supports
using that interval provisionally for the next short study, not for arbitrary
long-time conditions or other materials.

Next: compare 10 s updated mechanics every 0.5 s against frozen mechanics energy,
with identical initial energy transfer. Do not compare directly against an older
static-run energy mapping: that would confound feedback with coordinate-transfer
differences. The new feedback_strength study isolates subsequent updates while
holding initial surface, initial sampled energy, mesh, and phase timestep fixed.
This remains experimental; the 22 C stress-free reference and homogeneous
fractional-history mixing assumptions have not received physical validation.
