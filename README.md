# ClockfieldRedux

Research proposal, 24 October 2026. No simulator or experimental results yet.

## The question

Can a persistent structure in a dynamical medium change independently constructed clocks by a common factor, so that an effective local time description becomes useful?

Agreement must be measured from dynamics. Assigning every process the same local timestep demonstrates an imposed common clock, not its emergence. A shared rate is a necessary diagnostic for the proposed scalar clockfield description, not sufficient evidence for general relativity.

## What the cold-atom experiment contributes

Giovanni Barontini's *Testing the problem of time with cold atoms*, Physical Review Research 8, L022047 (2026), constructs an operational internal time using coarse-grained entropy exchange in a trapped condensate. It motivates comparing internal clocks. It does not establish that entropy exchange determines gravitational time dilation or wave speed.

- Paper: https://arxiv.org/html/2509.07745v3
- Published article: https://doi.org/10.1103/1h9j-df4k
- Analogue-gravity background: https://arxiv.org/abs/gr-qc/0505065

A clock observable can stop changing while other observables continue evolving. In particular, constant von Neumann entropy does not imply a stationary quantum state. Entropy exchange, entropy production, coherent particle current, and state distinguishability must remain separate measured quantities.

## Three possible starting points

1. **Repair the supplied Gemini animation.** Useful for demonstrating a chosen refractive medium, but its assigned speed and particle force cannot establish emergent gravity.
2. **Test clocks in a specified conservative substrate.** Recommended first scientific direction. Derive evolution from one declared Hamiltonian or action, construct different probes, and measure their responses independently.
3. **Construct a constrained relational-time quantum model.** Addresses a different question: conditional dynamics relative to an internal clock. Renaming numerical time lambda is insufficient, and tracing out a stationary clock is not equivalent to conditioning on its reading.

These approaches answer different questions and should not be combined into one claimed result.

## Proposed first experiment: clock agreement

The substrate and its parameters remain to be selected before implementation. A two-component condensate-like Hamiltonian is one candidate, not an accepted physical theory of spacetime.

Use a homogeneous reference state and a controlled altered background. Initially a prescribed static background is acceptable for diagnosing probe response, provided it is explicitly labelled prescribed. A self-supported defect is a later, separate requirement.

Compare:

- frequency of a localized oscillatory mode;
- round-trip time of a propagating pulse over a specified distance;
- frequency of coherent exchange between two coupled modes.

For each clock, normalize its measured frequency to its own homogeneous reference. Compare the normalized rates at matched locations and background strengths. Use inverse round-trip time for the pulse clock. Measure probe amplitudes, dispersion, damping, and sensitivity to the probe itself.

Keep ruler conventions explicit. A cavity fixed in lattice coordinates and a cavity whose boundaries are material structures need not represent the same physical length. Report both where feasible; do not interpret a changing distance as clock slowing.

### Outcomes and what they mean

| Observation | Supported interpretation |
|---|---|
| Clock rates change differently | Process-dependent response; the proposed common scalar clock description fails in that regime. |
| Agreement only after multiplying all updates by a shared factor | Implemented common-rate control works; emergence not demonstrated. |
| Agreement for one fitted pair, but not a third clock | Fit does not transfer across mechanisms. |
| Agreement persists across independent probes and held-out backgrounds | Evidence for an approximate shared effective time in the tested sector. |
| Clocks agree but trajectories require unrelated force rules | Clock agreement alone has not unified motion and propagation. |

No numerical pass thresholds are registered yet. Select tolerances using numerical convergence and measurement uncertainty, freeze them before held-out evaluation, and report failures.

## Controls required for interpretation

- Homogeneous baseline and timestep/grid convergence.
- An explicitly imposed common-rate positive control, labelled as such.
- A selective perturbation that changes one clock mechanism, to verify the analysis can detect disagreement.
- Background strengths, probe frequencies, and locations excluded from fitting.
- Local references: no grid-wide percentile normalization acting as an instantaneous physical coupling.
- Energy and particle-number accounting appropriate to the declared equations.
- Separate coherent exchange from entropy measures; do not infer one from the other without deriving and validating the relationship.

Even a common slowdown can be built into a Hamiltonian by scaling every term identically. That is a valid model assumption, but not an explanation of universality. Report exactly which couplings produce any agreement and whether it survives independent changes to them.

## Audit of the supplied Gemini prototype

This audit refers to the pasted prototype discussed before this repository was initialized, not to code currently in this repository.

- Its supposed topological solitons are Gaussian masks with separately integrated coordinates. No topological stability is demonstrated.
- Attraction is explicitly supplied as acceleration proportional to minus the gradient of wave speed.
- Its bridge defines Psi = -log(q)/alpha and then c = c0 exp(-alpha Psi), so c = c0 q. Alpha cancels from the dynamics through this path.
- The dynamic version's explicit two-dimensional phase-diffusion ratio is 0.4116; its stability limit is 0.25. Executing its core produced a maximum absolute phase of approximately 5.96e10 after the 35-step warm-up. Bounded sine outputs do not establish numerical stability.
- Its proposed metric ds^2 = -c0^2 q^2 dt^2 + q^-2 dx^2 gives null coordinate speed c0 q^2, whereas the code uses c0 q.
- Diffusive wakes and shrinking orbits do not by themselves demonstrate gravitational radiation or frame dragging.

Repairing these issues would make the prototype more interpretable; it would not turn its assumed coupling laws into derived physics.

## The longer-term requirement

If clock agreement survives, next ask whether propagation, material rulers, and localized excitations' motion can be described by a consistent effective geometry derived from the same substrate. Self-supported structures, backreaction, conservation laws, and the range of validity must then be measured.

The target is a falsifiable model of shared internal time. Emergent gravity remains a hypothesis, not a result.
