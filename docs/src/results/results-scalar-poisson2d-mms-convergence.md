# Poisson2D MMS Convergence

2D MMS Poisson, homogeneous Dirichlet on all borders.
u_exact = sin(pi*(x-dx)/(Lx-dx)) * sin(pi*(y-dy)/(Ly-dy))
-Laplace(u) = f = ( (pi/(Lx-dx))^2 + (pi/(Ly-dy))^2 ) * u_exact

**CSV source:** `results/scalar/Poisson2D_MMS_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Poisson2D_MMS_Convergence | 0.125 | L^2 | 64 | [8, 8] | 0.00741801 | 0.00741801 | 0 | 0 | NaN | NaN | NaN |
| Poisson2D_MMS_Convergence | 0.0625 | L^2 | 256 | [16, 16] | 0.00171724 | 0.00171724 | 0 | 0 | 2.11094 | 2.11094 | NaN |
| Poisson2D_MMS_Convergence | 0.03125 | L^2 | 1024 | [32, 32] | 0.000414763 | 0.000414763 | 0 | 0 | 2.04973 | 2.04973 | NaN |
| Poisson2D_MMS_Convergence | 0.015625 | L^2 | 4096 | [64, 64] | 0.000102005 | 0.000102005 | 0 | 0 | 2.02365 | 2.02365 | NaN |
| Poisson2D_MMS_Convergence | 0.0078125 | L^2 | 16384 | [128, 128] | 2.52981e-05 | 2.52981e-05 | 0 | 0 | 2.01154 | 2.01154 | NaN |
