# A vortex as an object, ruler and clock

The full complex field supports a winding-one vortex that is linearly stable in
the tested fixed-boundary boxes and survives the tested disturbances. Its
breathing clock is a candidate at tau=0 and 0.1; no bound radial breathing mode
is resolved at tau=0.4. Its characteristic sizes and clock ratios change when
the global coupling changes. This is progress toward a toy world containing
objects, not a demonstration of universal matter geometry or gravity.

[Measurements](vortex_results/SUMMARY.md) · [Declared protocol](VORTEX_PROTOCOL.md)
· [Raw results](vortex_results/results.json)

## Why a vortex instead of a real kink?

The physical field is complex. A real kink can unwind through its imaginary
component; testing only real perturbations misses that instability. The test
suite checks for a negative transverse kink eigenvalue at tau=0, where the
infinite-line result is omega²=-mu²/2. Negative squared frequencies must remain
negative: taking sqrt(abs(omega²)) disguises growth as an oscillation.

A two-dimensional vortex instead winds once around the complex vacuum circle.
Here its outer boundary fixes that winding. The centre is solved from the
nonlinear equation, not inserted as a permanent Gaussian hole or frozen mask.
Both field components are free to move at every interior point.

This is a global vortex, with long-range phase gradients. The fixed boundary
supplies its topological sector. It is not a finite-energy, isolated particle
in an infinite plane, and winding alone is not a stability proof. Its amplitude
approaches the vacuum algebraically; r90 is a 90%-amplitude radius, not an
exponential tail length or the kink's tail scale.

## Same legacy force, new boundary conditions

The source is [ClockfieldBornRule/clockfield_core.py](https://github.com/anttiluode/ClockfieldBornRule/blob/main/clockfield_core.py),
inspected blob `f4cc9fe9e72bc138c874337bf9e8c10f4a521f13`. As in the preceding
homogeneous assay, damping is explicitly removed. The original damped dynamics
could erase an oscillator; these results do not establish an operational clock
in that version. The source's periodic boundaries are replaced by a square with
fixed phi=v exp(i theta) on its edge, v²=mu²/lambda.

Let beta=|phi|². The exact undamped force can be rewritten as

$$\ddot\phi=a(\beta)[\nabla^2\phi+F(\beta)\phi],$$

$$a=(1+\tau\beta)^{-5}>0,\qquad F=(1+\tau\beta)(\mu^2-\lambda\beta).$$

We keep mu²=1.4 and lambda=0.55, and vary tau=0, 0.1, 0.4. These are different
global couplings, not different positions in one gravitational solution. The
source default tau=5 is outside this sweep; no vortex claim is made there.

## Full complex stability

Static equilibria solve R=lap(phi)+F phi=0 by sparse Newton iteration. The
Jacobian of -R is a real symmetric matrix K acting on both components. Writing
phi=x+i y, its local block is

$$K=-\nabla^2 I-FI-2F'\begin{pmatrix}x^2&xy\\xy&y^2\end{pmatrix}.$$

At equilibrium the derivative of a multiplies R=0, so the dynamical eigenproblem
is omega² delta=a K delta. This operator is similar to the symmetric matrix
sqrt(a) K sqrt(a). Congruence preserves the number of negative directions:
the lowest eigenvalues of K therefore test the **sign** of the full dynamical
spectrum. K's numerical eigenvalues are not themselves clock frequencies.

The full square grid admits angular and nonradial perturbations. No rotational
symmetry restriction is imposed on this stability test. The eigensolver finds
the three smallest algebraic eigenvalues, keeps their signs and checks residuals.
The lowest is positive in all nine cases. It decreases when the box is enlarged,
consistent with weaker boundary confinement; no claim of a nonzero infinite-box
translation gap follows. A separate force-Jacobian test checks both components
against finite differences.

There is a static functional whose Hessian is K. This does not establish a
conserved energy for the original nonlinear dynamics: a configuration-dependent
kinetic metric would generally add velocity terms absent from this equation.

## Size and a candidate breathing clock

The radius at 50% vacuum amplitude is r50; the radius at 90% is r90. Both are
measured on the positive horizontal axis of the solved square-grid vortex.
They are geometric sizes, not yet a demonstrated rigid rod or force-balanced
distance between two objects.

For an independent spectral calculation, solve phi=f(r) exp(i theta):

$$f''+f'/r-f/r^2+F(f^2)f=0.$$

The radial amplitude perturbation obeys

$$\omega^2u=a[-u''-u'/r+u/r^2-(F+2f^2F')u].$$

A cell-centred radial finite-volume discretization handles the origin without
an artificial hole. Transforming z=sqrt(r/a)u gives a symmetric tridiagonal
eigenproblem. The calculation tests radii 16 and 32 and spacings 0.025 and
0.0125. This radial-sector eigenproblem supplies a candidate clock; the separate
full 2D Hessian supplies the stability test.

The squared continuum radial threshold is omega_vac²=Gamma0² 2mu². A bound
candidate must lie at least 1% below that squared threshold, remain positive,
and change by less than 2% under radial box/grid refinement. Its field-amplitude
norm fraction inside 4 r50 is also reported. It is about 89% at tau=0, 56% at
tau=0.1, and only 1.6% at tau=0.4 in the radius-16 radial box. The latter drops
to about 0.11% when that box doubles: the mode is extended.

At tau=0 and 0.1 the spectral frequency / vacuum frequency is respectively
0.9019 and 0.9841. In full 2D evolution, the corresponding measured ratios in
the larger square box are 0.90145 and 0.98405. These are different from the real
kink's numbers: a vortex is a different object. At tau=0.4 the lowest radial
mode lies above threshold. Its box oscillations are recorded but explicitly
not counted as an internal clock. This finite search does not rigorously exclude
arbitrarily shallow modes or other angular-sector resonances.

The radial amplitude sector decouples from the gapless phase sector at linear
order in a rotationally symmetric continuum. That is a limitation of calling
this a clock: finite-amplitude radiation, symmetry-breaking disturbances and
long-term persistence need further tests.

## Time evolution and controls

The full nonlinear, two-component square-grid force is integrated by velocity
Verlet in original time units. A small radial-mode disturbance is prepared with
zero initial velocity. Its weighted projection is sampled; at least four upward
crossings are required to report a period. Initial preparation uses a computed
mode, but the reported time-domain frequency comes from evolved measurements.
An analytic vacuum lattice-mode test independently verifies the time integrator
and period estimator.

A separate test adds seeded complex noise of peak amplitude 0.001 v, localized
near the core. Both experiments run for eight vacuum radial periods. The noise
RMS stays below its starting value and the interior contour retains winding one
in every run. This tests one seed, small disturbances, zero initial velocities,
and a finite duration; it is not arbitrary-perturbation nonlinear stability.

All nine square cases pass the recorded numerical checks. Halving dx changes
the reported radii by at most 0.30%; doubling the box changes the measured
frequency by at most 1.20% (the largest change belongs to the extended tau=0.4
mode). Halving the timestep or probe amplitude changes frequencies by less than
0.005%. Raw traces and eigenvalue residuals are retained. Tightening these gates
or probing a wider parameter range would be a new experiment.

## What an inhabitant could notice

The vacuum wave scales obey C=Gamma0 c_eff,0,
omega_vac=Gamma0 sqrt(2mu²), and xi=c_eff,0/sqrt(2mu²), so
xi omega_vac / C=1. That cancellation calibrates the linear wave sector.

However, tau v² is dimensionless. Changing it can change internal structure;
it is not generally an invisible change of units. Here r50/xi rises from 1.347
to 1.598, r90/xi rises from 3.919 to 4.266, and the breathing/vacuum clock ratio
changes before the localized mode is no longer resolved. Even r90/r50 changes.
One common rescaling does not describe every measured property of this object.

This does not by itself violate the equivalence principle: globally changing a
coupling is not placing identical matter in different gravitational potentials.
No sourced background, free-fall comparison, light bending, Einstein dynamics,
or universal material metric has been demonstrated.

The result is narrower and useful: the model now has an explicitly constructed
defect with tested stability, measured sizes, and a candidate internal clock in
part of its parameter range. The next questions concern its lifetime, motion,
interactions, and clocks in a dynamically supported nonuniform background.

## Reproduce

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python run_vortex_experiment.py --out my_vortex_run
```

Open `my_vortex_run/report.html`. A complete run took about 31 seconds on the
recording machine with one BLAS thread; hardware and numerical libraries affect
runtime. No GPU is required. The script returns nonzero if a declared numerical
gate fails and still writes the measurements for inspection.
