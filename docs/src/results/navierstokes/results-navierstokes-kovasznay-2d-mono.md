# Kovasznay 2D Mono

Kovasznay flow convergence study (steady Navier–Stokes, Re ≈ 40) — CSV only.

Solves the manufactured Kovasznay solution on [0,1]² for a set of uniform grids,
enforcing exact Dirichlet data on all boundaries. Reports L2/Linf errors for
u, v, p and estimated pairwise convergence rates. Outputs are written to
`results/NavierStokes/Kovasznay_2D_Mono.csv` (no plots).

**CSV source:** `results/NavierStokes/Kovasznay_2D_Mono.csv`

| N | h | err_u_L2 | err_u_Linf | err_v_L2 | err_v_Linf | err_p_L2 | err_p_Linf | rate_u_L2 | rate_u_Linf | rate_v_L2 | rate_v_Linf | rate_p_L2 | rate_p_Linf |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 32 | 0.03125 | 0.00408245 | 0.0416132 | 0.00256848 | 0.0325682 | 1.55727e+10 | 6.03128e+10 | NaN | NaN | NaN | NaN | NaN | NaN |
| 64 | 0.015625 | 0.000148288 | 0.000423464 | 6.9532e-05 | 0.000432081 | 1.69798e+08 | 9.45394e+08 | 4.78296 | 6.61866 | 5.2071 | 6.23602 | 6.51906 | 5.99541 |
| 128 | 0.0078125 | 4.3145e-05 | 0.000138386 | 1.91839e-05 | 7.92033e-05 | 1.36733e+07 | 1.08528e+08 | 1.78113 | 1.61354 | 1.85778 | 2.44767 | 3.63439 | 3.12284 |
| 256 | 0.00390625 | 8.75004e-06 | 5.83649e-05 | 5.0024e-06 | 9.55056e-05 | 5.72585e+06 | 6.4527e+07 | 2.30183 | 1.24553 | 1.9392 | -0.270025 | 1.2558 | 0.750096 |
