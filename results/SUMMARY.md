# Measured clock comparison

Deterministic linear toy model; no emergent-gravity claim.

| Case | Oscillation | Pulse | Exchange | Agreement within 1% |
|---|---:|---:|---:|---|
| reference | 1.0000 | 1.0000 | 1.0000 | yes |
| imposed h=0.85 | 0.8500 | 0.8500 | 0.8500 | yes |
| imposed h=0.65 | 0.6500 | 0.6500 | 0.6500 | yes |
| imposed h=0.45 | 0.4500 | 0.4500 | 0.4500 | yes |
| stiffness c=0.85 | 1.0000 | 0.8049 | 1.0000 | no |
| stiffness c=0.65 | 1.0000 | 0.5474 | 1.0000 | no |
| stiffness c=0.45 | 1.0000 | 0.3048 | 1.0000 | no |
| candidate s=0.15 | 1.1500 | 0.7606 | 0.9136 | no |
| candidate s=0.3 | 1.3000 | 0.5316 | 0.8390 | no |
| candidate s=0.45 | 1.4500 | 0.3298 | 0.7742 | no |

Validation: {"energy": true, "sampling": true, "spatial": true, "imposed_control": true, "selective_control": true}.

The imposed Hamiltonian multiplier makes clocks agree by construction. All tested candidate backgrounds fail the common-rate criterion. Pulse rates use a fixed lattice ruler. See results.json for raw measurements and refinement checks.
