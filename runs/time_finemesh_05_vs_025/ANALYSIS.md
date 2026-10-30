# Fine-mesh time sensitivity at 100 s

The 200 × 440 run at dt=0.25 s completed 400 accepted steps, no rejected steps,
and no solver rebuilds. Actual steps remained 0.25 s; maximum Newton iterations
were 3. Final relative mass drift was 1.313e-14. All surface
crossings were found, and sampled nodal phi extrema stayed inside [0,1].

Refining time from 0.5 to 0.25 s changes the final surface by 0.001092 nm RMS,
with maximum difference 0.002115 nm. The center changes by 0.001191 nm.
This is about 60 times smaller than the previous 0.065416 nm mesh difference.

At matched dt=0.25 s, comparing the existing 100 × 220 and 200 × 440 runs gives:
- Final surface RMS difference: 0.065438 nm.
- Maximum difference: 0.231366 nm.
- Center, coarse minus fine: -0.231366 nm.
- Edge, coarse minus fine: 0.212675 nm.

Thus time refinement does not remove the spatial sensitivity. The fine-mesh
center rises by 0.033132 nm from its
initial height; its final height remains 2.079007 nm below the deformed dielectric
edge. The small upward movement is robust to this time-step change but is not
yet established under spatial refinement.

Next: a third phase mesh (400 × 880; nominal 0.625 nm spacing) at dt=0.25 s,
using configs/model_mesh_finer_dt025.json. It has four times the phase cells of
the 200 × 440 run and will require more memory and solve time; no production run
was launched during this analysis. Compare with the existing 200 × 440 result.
The mechanics mesh, 5 nm interface-width parameter, and physical inputs remain
fixed. These comparisons include phase-discretization and transferred-energy
representation sensitivity, not a refinement of the mechanics solution.

Three meshes will let us assess whether spatial differences decrease. They do
not automatically establish an asymptotic convergence regime, and time sensitivity
may need reevaluation if remaining spatial differences approach 0.001–0.002 nm.
