# Version 1 protocol

Frozen before the parameter sweep was executed. Engineering validation thresholds, not statistical significance or a blinded discovery trial. Model choices and controls were already known analytically.

- Reference: m=1, c=1, coupling=0.8, h=1, L=128, N=256.
- Imposed common-rate controls: h=0.85, 0.65, 0.45, all other parameters fixed.
- Selective stiffness controls: c=0.85, 0.65, 0.45, all other parameters fixed.
- Candidate backgrounds: m=1+s, c=1-s, coupling=0.8, h=1, s=0.15,0.30,0.45. These are prescribed parameter changes, not self-supported defects or density laws.
- Normalize each clock to its own reference rate. Agreement requires max(rate)/min(rate)-1 < 0.01. Positive controls must also match h within absolute 0.002. Selective controls must fail agreement.
- Candidate agreement is exploratory: no fit is performed and no held-out learning claim is made.
- Relative energy drift must be < 1e-10. Signals must have >=4 upward crossings and circular packet coherence >0.5.
- Sampling check: 2048 versus 4096 samples, every case, rate changes <0.002 relative.
- Spatial check: reference and strongest candidate, N=256,512,1024,2048 at fixed physical length, width, carrier. Relative pulse-rate difference between the final two grids <0.01. Exact-in-time evolution removes integration timestep error; sampling and spatial errors remain.
- Pulse clock is measured centroid speed divided by the fixed ring circumference, an inferred circulation rate. It is not a measured cavity round trip or a material-ruler clock.
