# Saved-field energy-transfer audit

No mechanics or phase-field simulation was run. This analysis uses the existing
native and ANSYS-driven 100 s runs on their common phase grid. The ANSYS source
interpolation was reconstructed and reproduces saved interior energy values
exactly. A small boundary tolerance excludes nominal-boundary nodes perturbed
by floating-point roundoff; this makes the interior metric 10.11%, rather than
the earlier 10.07% calculated with strict coordinate inequalities.

| Region | Nodal relative RMS energy mismatch | Absolute RMS (J/m³) |
|---|---:|---:|
| Interior Cu | 10.11% | 37059 |
| Upper 10 nm of reference Cu | 119.62% | 61662 |
| Upper layer, r < 25 nm | 39.22% | 2142 |
| Upper layer, r = 25–225 nm | 34.07% | 8216 |
| Upper layer, r > 225 nm | 130.97% | 192607 |

The outermost 25 nm accounts for 98.55% of the upper-layer sum of squared
energy differences. Thus the 120% headline mismatch is dominated by the high-energy
outer corner, rather than a uniform 120% error across Cu. Nevertheless, the middle
region still differs by 34% relatively, despite much smaller absolute energies.
These are unweighted nodal statistics, not axisymmetric volume integrals.

No interior or upper-layer sample in either input is spuriously zero from being
outside the interpolator's convex hull. Switching only the ANSYS source mask
from the displaced initial top to the reference Cu top removes three source
points, but changes the upper-layer interpolated energy by less than 3e-11 J/m³
at any sampled node. That mask choice does not explain this sampled discrepancy.
This test does not cover every possible quadrature location or the field above
the reference surface.

Remaining plausible contributors are corner mechanical resolution, nodal ANSYS
postprocessing versus cellwise native DG0 energy, and subsequent P1 interpolation.
Saved fields do not separate those effects or prove which is dominant. The shared
bulk strain-energy definition was checked earlier; this audit does not justify
an empirical energy rescaling or changes to the material law.

The next work should inspect corner profiles and interpolation/averaging at fixed
existing resolution. If evolution sensitivity must be isolated, a short energy-off
comparison on the routine mesh is cheaper than further uniform refinement. Such
an ablation measures sensitivity, not which mechanical field is correct. No new
simulation is required from the user at this point.

Reproduce with:
`python new/analyze_energy_transfer.py`
