# Original Clockfield PDE: homogeneous-clock protocol

Declared before running this assay. This is a deterministic numerical test of an analytically predicted disagreement, not a blinded discovery trial.

Source: https://github.com/anttiluode/ClockfieldBornRule/blob/main/clockfield_core.py

Use the full nonlinear complex-field force from that source with mu²=1.4, lambda=0.55, c0²=1. A one-dimensional periodic slice is equivalent to its two-dimensional equation for fields uniform in the other direction. The equilibrium is v=sqrt(mu²/lambda). This test does not insert a vortex.

Main assay: omit damping to isolate the clock-rate claim. Evolve two independent weak perturbations: uniform radial displacement and a travelling phase packet. Recompute beta, Gamma and c_eff from the evolving field; do not freeze coefficients or substitute analytic frequencies in the measurements.

Tau sweep: 0, 0.1, 0.2, 0.4, 1, 5. Default epsilon=1e-4 times equilibrium amplitude. Domain length 128, N=256, packet width 8, carrier 2*pi*12/128. Integrate using velocity Verlet in scaled numerical time s=Gamma_equilibrium*t; this constant change of units is undone for all reported physical rates. Default ds=0.025. Sample every approximately 0.1 units of s. Pulse travel measured by Hilbert-envelope circular centroid using a fixed lattice ruler. Radial frequency measured by upward zero crossings. Require >=4 crossings, finite bounded fields, and pulse centroid coherence >0.5.

Prediction relative to tau=0: radial rate Gamma0; pulse rate Gamma0/sqrt(1+tau*beta0). Clock ratio pulse/radial = 1/sqrt(1+tau*beta0). Common-rate tolerance 1%; prediction relative error tolerance 0.5%.

Controls: at tau=0.4 and 5 halve ds, halve epsilon, and double N. Each rate must change by <0.5%. Check one-dimensional full-field acceleration against the two-dimensional stencil for a y-uniform field, equilibrium, and the finite-difference linearization.

Original damping audit: separately run the source's actual uniform scalar recurrence with dt=0.02, damping=0.005, 3000 steps, epsilon=1e-4. This coefficient corresponds approximately to a continuous friction rate 0.005/0.02=0.25, not 0.005. Report observed cycles and the linearized damping regime. Do not infer a measured ticking frequency where the window lacks enough cycles or the oscillator is overdamped.

Scope: fixed potential and equilibrium, varying model coupling tau. No material rulers, self-supported defects, localized bound clocks, mass measurements or equivalence-principle test. This homogeneous assay must be interpreted before a defect experiment is designed.

## Supplemental refinement (reported separately)

The 256→512 raw pulse-speed check fails the frozen 0.5% tolerance. A 512→1024 refinement is also reported, without deleting or relabelling the failed coarse-grid check. The published summary distinguishes the original failure from the supplemental finer-grid pass. No continuum-extrapolated rate is claimed.
