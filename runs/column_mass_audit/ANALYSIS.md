# Routine-mesh mass-transfer audit

Reused the saved initial/final P1 phase fields from `runs/native_dt05_100s`.
No phase or mechanics simulation was run. The original phase mesh has 100 by
220 rectangles (44,000 triangles), and the trajectory uses static one-way energy.

| Radial projection nodes | RMS versus contour (nm) | Maximum versus contour (nm) | RMS versus 401 nodes (nm) |
|---|---:|---:|---:|
| 51 | 0.06357 | 0.18615 | 0.00317 |
| 101 | 0.06356 | 0.18848 | 5.63e-11 |
| 201 | 0.06356 | 0.18848 | 1.16e-11 |
| 401 | 0.06356 | 0.18848 | 0 |

RMS values use uniformly spaced radii, not area weighting. These compare changes
from the initial surface, retaining the same initial mechanics geometry.
The maximum discrepancy is near the pad edge. At 201 projection nodes, the
mass-equivalent center changes by -0.10609 nm and the edge by -7.07555 nm.
All resulting reference heights stay inside the supported mechanical domain.

The contour-derived reference volume changes by -0.009270%, while independent
integration of the saved phase change is zero to roundoff. Every projection
preserves that integrated phase change to roundoff. The 101-, 201-, and 401-node
profiles agree numerically for this field; this is radial projection insensitivity
on the saved mesh, not phase-field spatial convergence. The existing 201-node
sampling can remain unchanged. Refining the phase mesh is not required for this
transfer-resolution check.

The new map is mass conserving, but it defines a different surface from the
phi=0.5 contour. That difference is modest in absolute size on this routine mesh,
but can matter for subnanometer conclusions and for energy near the Cu edge.
It must not be described as physical validation or a correction to saved results.

Only t=0 and t=100 s phase fields were saved. These endpoints cannot determine
plastic-history evolution or feedback update-interval sensitivity. Next is a
short controlled two-way pilot on the routine mesh, comparing mechanics update
intervals before any long-duration run. No finer phase mesh is warranted here.

Files: `audit.json`, `surface_changes.csv`, `surface_comparison.png`.
