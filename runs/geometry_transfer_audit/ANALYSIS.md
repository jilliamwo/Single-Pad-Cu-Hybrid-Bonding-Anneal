# Geometry transfer audit

The existing native 100 s run was postprocessed on its existing mechanics grid.
No mechanical or phase-field solve was performed.

- The new axisymmetric cell fractions reproduce the supplied reference surface
  graph volume to approximately 1.1e-16 relative error.
- Fractions change in 97 cells, including 86 partly filled cells. Centroid binary
  occupancy changes in only 26 cells.
- Binary final volume is 0.009157% below the reconstructed graph volume.
- The reconstructed graph itself loses 0.009270% volume relative to its initial
  reference geometry, despite diffuse phase mass being conserved to ~3.7e-15.

The cell-fraction method removes the additional graph-to-grid volume error and
captures subcell movements without mesh refinement. It does not remove the
separate diffuse-field-to-graph discrepancy. No shape correction was imposed.

The fractions currently serve as a candidate geometry map and diagnostics in
experimental mechanical trials. Binary mechanics remains active until fractional
Cu/history constitutive treatment is implemented and tested. No new user run is
needed yet; the audit is reproducible with `python new/audit_geometry_transfer.py`.
