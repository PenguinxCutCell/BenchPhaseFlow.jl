# TaylorGreen 2D Mono

Taylor–Green vortex convergence for the Navier–Stokes prototype (no plotting).

Solves the manufactured Taylor–Green flow on [0, 2π]² with exact velocity
Dirichlet boundaries, computes weighted L2 errors for u and v, estimates
pairwise convergence rates, and writes a CSV to `results/NavierStokes`.

**CSV source:** `results/NavierStokes/TaylorGreen_2D_Mono.csv`

| h | error_u | error_v | rate_u | rate_v |
| --- | --- | --- | --- | --- |
| 0.785398 | 0.0125523 | 0.0125523 | NaN | NaN |
| 0.392699 | 0.00372066 | 0.00372068 | 1.75432 | 1.75431 |
| 0.19635 | 0.00101272 | 0.00101261 | 1.87733 | 1.87749 |
| 0.0981748 | 0.000261897 | 0.000261857 | 1.95116 | 1.95123 |
| 0.0490874 | 6.64827e-05 | 6.64723e-05 | 1.97795 | 1.97795 |
| 0.0245437 | 1.67361e-05 | 1.67335e-05 | 1.99002 | 1.99002 |
