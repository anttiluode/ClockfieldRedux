# Vortex: object, ruler, clock

Full complex legacy field, undamped. Fixed winding-one outer boundary. Finite-box stability and finite-duration survival are not infinite-domain stability proofs.

| tau | r50 / xi | r90 / xi | Radial mode / vacuum | Bound clock candidate? | Full complex linear stability |
|---|---:|---:|---:|---|---|
| 0 | 1.3474 | 3.9190 | 0.9019 | yes | pass |
| 0.1 | 1.4269 | 4.0475 | 0.9841 | yes | pass |
| 0.4 | 1.5977 | 4.2659 | 1.0050 | not resolved | pass |

The spectral ratio in the table is a radial-sector calculation. At tau=0.4 the lowest reported radial mode is extended and above the vacuum threshold; it is not an internal clock. Full 2D time measurements and all convergence checks follow.

| tau | dx | half-width | Measured omega / vacuum | Prediction error | Noise RMS growth |
|---|---:|---:|---:|---:|---:|
| 0 | 0.25 | 8 | 0.90148 | 0.047% | 1.000 |
| 0 | 0.125 | 8 | 0.90187 | 0.004% | 1.000 |
| 0 | 0.25 | 16 | 0.90145 | 0.051% | 1.000 |
| 0.1 | 0.25 | 8 | 0.98519 | 0.112% | 1.000 |
| 0.1 | 0.125 | 8 | 0.98525 | 0.118% | 1.000 |
| 0.1 | 0.25 | 16 | 0.98405 | 0.004% | 1.000 |
| 0.4 | 0.25 | 8 | 1.01680 | 1.178% | 1.000 |
| 0.4 | 0.125 | 8 | 1.01614 | 1.112% | 1.000 |
| 0.4 | 0.25 | 16 | 1.00467 | 0.028% | 1.000 |

## Checks

- tau=0, half grid spacing: r50 0.011%, r90 0.212%, oscillation 0.043% — pass (2% limit).
- tau=0, double box: r50 0.000%, r90 0.001%, oscillation 0.003% — pass (2% limit).
- tau=0.1, half grid spacing: r50 0.275%, r90 0.291%, oscillation 0.006% — pass (2% limit).
- tau=0.1, double box: r50 0.000%, r90 0.000%, oscillation 0.116% — pass (2% limit).
- tau=0.4, half grid spacing: r50 0.151%, r90 0.092%, oscillation 0.065% — pass (2% limit).
- tau=0.4, double box: r50 0.000%, r90 0.000%, oscillation 1.193% — pass (2% limit).
- tau=0: half-step 0.0044%, half-amplitude 0.0000% (1% limit).
- tau=0.1: half-step 0.0021%, half-amplitude 0.0000% (1% limit).
- tau=0.4: half-step 0.0003%, half-amplitude 0.0000% (1% limit).

All declared numerical checks: **PASS**.
A successful stability check is not a successful universal-clock test. Changing tau changes dimensionless internal structure. A characteristic wave length is not automatically a material ruler.

An above-threshold radial mode is a box oscillation, not a bound-clock candidate; its frequency agreement gate is still recorded, never silently waived.

See [VORTEX_MODEL.md](../VORTEX_MODEL.md), [protocol](../VORTEX_PROTOCOL.md), and [raw results](results.json).
