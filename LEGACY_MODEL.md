# Direct test of the original Clockfield force

## Result

The nonlinear simulation confirms the predicted mismatch between a radial oscillator and phase travel with a fixed lattice ruler. All five nonzero tested couplings fail the 1% common-rate criterion. At tau=0.4, normalized radial and travel rates are 0.245516 and 0.172820. At the source default tau=5, they are 0.005307 and 0.001432.

This is evidence about this equation in a homogeneous, weak-probe regime. It is not a general disproof of emergent spacetime, and it does not test inertial versus gravitational mass.

**Interpretation update:** with the characteristic wave length
xi=c_eff,0/sqrt(2mu²), the dimensionless ratio xi omega_radial / C equals one,
where C=Gamma0 c_eff,0. Thus the fixed-ruler mismatch does not rule out a common
linear wave geometry. A characteristic length is not automatically a constructed
material rod. The subsequent [vortex assay](VORTEX_MODEL.md) measures an actual
defect's sizes and candidate breathing clock. Changing the global dimensionless
coupling tau v² can be internally observable; it is not a gravitational redshift
experiment or necessarily a mere change of units.

## Source and fidelity

Source: [ClockfieldBornRule/clockfield_core.py](https://github.com/anttiluode/ClockfieldBornRule/blob/main/clockfield_core.py), inspected source blob `f4cc9fe9e72bc138c874337bf9e8c10f4a521f13`.

We retain the full state-dependent force, with mu²=1.4, lambda=0.55 and c0²=1:

$$\partial_t^2\phi=\Gamma(\phi)^2\left[c_{eff}(\phi)^2\nabla^2\phi+\mu^2\phi-\lambda|\phi|^2\phi\right],$$

$$\Gamma=(1+\tau|\phi|^2)^{-2},\qquad c_{eff}^2=(1+\tau|\phi|^2)^{-1}.$$

The implemented one-dimensional periodic slice is the original two-dimensional equation restricted to fields independent of y. A test compares its force directly against a two-dimensional five-point stencil. Unlike the source's vortex initialization, this assay starts near the homogeneous vacuum. The main assay removes damping explicitly; a separate function reproduces the source's damped scalar recurrence. It is not a reproduction of the source's vortex run.

Here phi is the original physical complex order parameter, with real and imaginary components. The two rows of the simulation array are independent probe experiments, not coupled fields. The code recomputes beta, Gamma and c_eff at every integration step. No linearization is substituted into its force.

## Prediction and measurement

Write phi = v + a + i b, with v²=mu²/lambda. To first order about equilibrium, neglecting damping:

$$\ddot a=\Gamma_0^2(c_{eff,0}^2\nabla^2a-2\mu^2a),\qquad\ddot b=\Gamma_0^2c_{eff,0}^2\nabla^2b.$$

The radial zero-wave-number mode oscillates at Gamma0 sqrt(2 mu²). Long-wavelength phase propagation has speed Gamma0 c_eff,0. Relative to tau=0:

$$r_a=\Gamma_0,\qquad r_b=\Gamma_0 c_{eff,0},\qquad r_b/r_a=(1+\tau v^2)^{-1/2}.$$

A finite-grid packet has a lattice group-velocity correction. Using the same carrier, width and grid in each tau run and its reference cancels that correction in the linear limit; raw speeds still need spatial refinement.

The oscillator is measured from interpolated upward zero crossings. Travel is measured from the circular centroid of the phase packet's spatial Hilbert envelope. The pulse velocity is initialized using the right-moving lattice phase branch, but thereafter both components evolve under the full nonlinear force. A slope fit supplies the measured speed. The fixed-circumference pulse-clock frequency has the same normalized rate as that speed.

To avoid an unnecessarily long computation when Gamma0 is small, integrate in s=Gamma0*t. This is a constant change of numerical units for each run; the entire force is divided by Gamma0² and all observation times are converted back to t. It does not impose equal measured rates on the two probes.

## Numerical results and the failed coarse-grid check

The six-case sweep agrees with predicted normalized rates to a maximum relative error of 2.54e-5 (about 0.00254%). At tau=0.4 and 5:

- Halving the integration step changes rates by at most about 0.0055%.
- Halving the probe amplitude changes rates by less than 0.000006%.
- Refining 256 to 512 sites changes raw pulse speeds by about 0.839%, **failing** the declared 0.5% criterion.
- A supplemental refinement from 512 to 1024 sites changes raw pulse speeds by about 0.2085%, below that tolerance.

The original coarse-grid gate therefore did not pass in full. The follow-up refinement supports convergence while preserving that failure in the JSON and summary. Main displayed normalized rates are the 256-site values, not continuum-extrapolated results. No thresholds were relaxed. The script reports `first_grid_refinement: false`; exit success means the prediction, step, amplitude and supplemental finer-grid checks passed, not that every original gate passed.

No conserved energy is claimed for this state-dependent, non-divergence-form force. The Hamiltonian conservation tests in the original Redux suite apply to that different model only.

## Original damping matters

The source recurrence includes `-0.005*(q-q_previous)` at dt=0.02. Its continuous-time friction approximation is 0.005/0.02=0.25. Treating the coefficient as friction 0.005 would change the model.

The separate audit uses precisely that scalar recurrence for the source's 3000 steps (60 time units). At tau=5, the linearized radial mode is overdamped; the measured probe has no zero crossings. Tau=1 is underdamped in the linearized equation but gives no upward crossings in that short window. These are different situations, and neither supports an extracted clock frequency from that run.

## What remains open

Two fixed-ruler clocks fail to share a single multiplicative slowdown in the undamped homogeneous assay. That directly challenges identifying Gamma alone with every clock rate in this construction. It does not exclude a richer effective geometry with changing spatial rulers: material lengths can change, and bound-object clocks must be constructed and measured rather than assumed. Vortices, spatially varying backgrounds, nonlinear bound clocks and moving rulers are not tested here.

## Run

```bash
python run_legacy_experiment.py
python -m unittest discover -s tests -v
```

Open `legacy_results/report.html` locally, or read [the result table](legacy_results/SUMMARY.md). The full numerical checks, parameters, measured rates and damping audit are in [results.json](legacy_results/results.json).
