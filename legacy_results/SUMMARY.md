# Original Clockfield PDE: measured clocks

Full nonlinear force; homogeneous equilibrium; weak probes; main assay damping disabled.

| tau | Radial clock | Phase travel | Travel / radial | Predicted ratio | Common rate? |
|---|---:|---:|---:|---:|---|
| 0 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | yes |
| 0.1 | 0.635371 | 0.567259 | 0.892800 | 0.892805 | no |
| 0.2 | 0.439106 | 0.357444 | 0.814026 | 0.814034 | no |
| 0.4 | 0.245516 | 0.172820 | 0.703904 | 0.703914 | no |
| 1 | 0.079553 | 0.042249 | 0.531075 | 0.531085 | no |
| 5 | 0.005307 | 0.001432 | 0.269896 | 0.269903 | no |

Validation: {"prediction": true, "half_step": true, "half_amplitude": true, "first_grid_refinement": false, "second_grid_refinement": true}.

The first spatial refinement is retained even if it fails the declared 0.5% raw-rate tolerance. The second refinement tests whether that discrepancy shrinks. The normalized clock disagreement is far larger than the reported numerical errors.

## Original damping audit

| tau | Linearized regime | Upward crossings in original 60 time units | Measured frequency |
|---|---|---:|---:|
| 0 | underdamped | 16 | 0.2659167538969575 |
| 0.1 | underdamped | 10 | 0.16824737973492934 |
| 0.2 | underdamped | 7 | 0.1153778380604842 |
| 0.4 | underdamped | 3 | not resolved |
| 1 | underdamped | 0 | not resolved |
| 5 | overdamped | 0 | not resolved |

The original coefficient 0.005 at dt=0.02 approximates friction 0.25. An unresolved or overdamped amplitude probe is not a zero-frequency universal clock. This is not a vortex experiment or a material-ruler comparison.
