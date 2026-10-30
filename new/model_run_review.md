# Review of the saved native 100-second run

Reviewed 2026-10-14. Source: `runs/native_100s` (preserved unchanged).

- 44 accepted steps of 1 s; last accepted time 44 s.
- The six failed attempts start at the same time, with dt reduced from 1 s to
  0.03125 s. All were recorded as -999 by the old catch-all PETSc handler.
- The recorded exception is PETSc error 73: a read-locked vector in
  `SNESLineSearchApply_BT` / `MatMult`.
- Accepted steps immediately before failure required only two Newton iterations.
  This does not establish convergence of the next step, but the exception alone
  is not evidence of a physical instability or poor material calibration.
- `time_history.csv` and `final_surface.csv` only describe t=0. No 44 s profile
  was saved because the failure preceded the first 50 s measurement.

PETSc describes this error as an attempt to write a vector while it is locked:
[VecSetErrorIfLocked source documentation](https://petsc.org/release/include/petscvec.h.html).
The exact source of the retained lock is unresolved. A coarse 60-step test did
not reproduce it. No dependency upgrade, material change, or phase-field equation
change was made to try to conceal the failure.

The runner now restores the last accepted solution, saves diagnostics, and
attempts one fresh solver construction after a PETSc exception. A repeated
exception is fatal, with full messages saved in `solver_errors.json`. Normal
SNES non-convergence retains adaptive time-step reduction.

Fault-injection verification passes: corrupting the tentative solution and
raising error 73 on the second solve produces the same final surface as the
uninterrupted run after recovery. Repeated errors stop after one rebuild and
preserve the accepted state. This validates recovery logic, not resolution of
the original full-mesh error. The user should run the documented 100 s retry.

The thermal reference remains 22 C as an explicitly unverified input assumption.
It has not been inferred from the annealing temperatures or newly confirmed.
