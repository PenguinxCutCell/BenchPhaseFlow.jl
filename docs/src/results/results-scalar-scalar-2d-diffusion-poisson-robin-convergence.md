# Scalar 2D Diffusion Poisson Robin Convergence

Steady 2D Poisson problem with an embedded circular Robin boundary.

The domain is [0, 4]^2 with a circular obstacle of radius `ly / 4` centered at
`(lx/2, ly/2) + (0.01, 0.01)`. Dirichlet zero data is imposed on the outer box
and Robin(1, 1, 1) data on the embedded boundary. The manufactured solution
`u(x, y) = 7/4 - ((x - cx)^2 + (y - cy)^2)/4` is used to measure both L2 and H1
(gradient) convergence.

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Poisson_Robin_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut | h1_all_err | h1_full_err | h1_cut_err | h1_empty_err | pair_order_h1_all | pair_order_h1_full | pair_order_h1_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_2D_Diffusion_Poisson_Robin | 0.25 | L^2 | 37 | [8, 8] | 0.00541294 | 0.0038831 | 0.00377113 | 0 | NaN | NaN | NaN | 0.0888754 | 0.0761775 | 0.0457802 | 0 | NaN | NaN | NaN |
| Scalar_2D_Diffusion_Poisson_Robin | 0.125 | L^2 | 173 | [16, 16] | 0.00121629 | 0.000881325 | 0.000838236 | 0 | 2.15392 | 2.13946 | 2.16957 | 0.0445567 | 0.0410268 | 0.017381 | 0 | 0.996142 | 0.892798 | 1.39721 |
| Scalar_2D_Diffusion_Poisson_Robin | 0.0625 | L^2 | 741 | [32, 32] | 0.000419804 | 0.000261555 | 0.000328366 | 0 | 1.53471 | 1.75256 | 1.35205 | 0.0229004 | 0.0212967 | 0.00841897 | 0 | 0.96027 | 0.945938 | 1.0458 |
| Scalar_2D_Diffusion_Poisson_Robin | 0.03125 | L^2 | 3090 | [64, 64] | 0.000120944 | 7.1036e-05 | 9.78838e-05 | 0 | 1.79538 | 1.88049 | 1.74616 | 0.0118053 | 0.0108932 | 0.00454992 | 0 | 0.955941 | 0.967195 | 0.887803 |
| Scalar_2D_Diffusion_Poisson_Robin | 0.015625 | L^2 | 12608 | [128, 128] | 3.47577e-05 | 1.69713e-05 | 3.03328e-05 | 0 | 1.79893 | 2.06545 | 1.69019 | 0.00579711 | 0.00548693 | 0.00187086 | 0 | 1.02603 | 0.989362 | 1.28214 |
| Scalar_2D_Diffusion_Poisson_Robin | 0.0078125 | L^2 | 50954 | [256, 256] | 1.00016e-05 | 4.0563e-06 | 9.14213e-06 | 0 | 1.7971 | 2.06486 | 1.73027 | 0.00290877 | 0.00276324 | 0.000908556 | 0 | 0.994924 | 0.98964 | 1.04205 |
