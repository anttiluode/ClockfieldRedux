# Vortex assay: declared plan

Written before running the vortex sweep. This is a deterministic numerical
experiment, not a blinded discovery trial.

Test a winding-one vortex of the undamped, full complex legacy equation at
tau = 0, 0.1, 0.4, with mu2 = 1.4 and lambda = 0.55. Tau is a global coupling
change, not a sourced gravitational potential. Fixed outer winding supports the
topological sector; an isolated global vortex has long-range phase structure.

1. Solve the static nonlinear equation on a square with fixed outer field.
   Measure radii at 50% and 90% vacuum amplitude. Require maximum dimensionless
   static residual < 1e-8.
2. Check the lowest eigenvalue of the full two-component static Hessian, whose
   inertia gives the sign of the linear dynamical spectrum. Keep negative signs.
   A lowest eigenvalue >= -1e-7 is a finite-box linear stability pass, not a proof
   of infinite-domain or nonlinear stability. Check eigenvector residual < 1e-7.
3. Independently solve the radial vortex and its amplitude-mode spectrum. A
   candidate clock must have positive squared frequency below 0.99 times the
   vacuum squared mass threshold; its frequency must change by < 2% on doubling
   radial domain and halving radial spacing. Report its spatial localization.
4. Excite the candidate radial mode in the full two-dimensional dynamics and
   measure its period from at least four upward crossings. Require measured
   frequency within 2% of the radial spectral prediction. Repeat with half
   timestep and half amplitude; require < 1% frequency changes.
5. Independently perturb both field components with seeded localized noise.
   Evolve for eight vacuum radial periods. Require winding retained and
   perturbation RMS growth < 5. This is a finite-duration survival test only.
6. Repeat square-grid sizes at dx = 0.25 and 0.125, half-width 8; repeat the
   coarse grid at half-width 16. Require radii and clock frequency changes < 2%.
   Retain failures; do not silently replace a failed gate with a finer result.

Report dimensionless core radius / xi, tail-region radius / xi, and candidate
clock frequency / vacuum frequency. A varying ratio can be observed internally;
only universal changes of units are unobservable. No gravity claim follows.

Implementation plan: add the field and spectrum solver, independent analytic
and finite-difference tests, a reproducible runner, raw results and offline visual
report; update the README interpretation without removing earlier results.
