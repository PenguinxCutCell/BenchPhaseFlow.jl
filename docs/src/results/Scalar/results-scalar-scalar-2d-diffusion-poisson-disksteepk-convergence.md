# Scalar 2D Diffusion Poisson DiskSteepK Convergence

2D variable-coefficient Poisson problem in a disk with a smooth steep k(r).
Solve -div(k grad u) = f on r <= R, u = 0 on r = R.
Exact solution (non-radial, m=2 by default):
    u(x,y) = (1 - r^2/R^2) * P_m(x,y), with P_m = Re((x + i y)^m).

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Poisson_DiskSteepK_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_2D_Diffusion_Poisson_DiskSteepK | 0.025 | L^2 | 1185 | [40, 40] | 0.00460489 | 0.00450753 | 0.000941922 | 0 | NaN | NaN | NaN |
| Scalar_2D_Diffusion_Poisson_DiskSteepK | 0.0125 | L^2 | 4865 | [80, 80] | 0.00184726 | 0.00182797 | 0.000266283 | 0 | 1.31778 | 1.3021 | 1.82265 |
| Scalar_2D_Diffusion_Poisson_DiskSteepK | 0.00625 | L^2 | 19789 | [160, 160] | 0.000714572 | 0.000710932 | 7.20366e-05 | 0 | 1.37024 | 1.36246 | 1.88616 |
| Scalar_2D_Diffusion_Poisson_DiskSteepK | 0.003125 | L^2 | 79813 | [320, 320] | 0.000270994 | 0.000270221 | 2.04519e-05 | 0 | 1.39882 | 1.39557 | 1.8165 |
