# Scalar 1D Diffusion Heat Moving ConstantBC Convergence

1D heat equation with a moving interval and constant Dirichlet boundary data.
The solution is manufactured to remain identically 1 when both the immersed
boundary and the domain edges enforce `u = 1`. This script verifies that the
numerical solution preserves this constant state while the cut location moves.

**CSV source:** `results/scalar/Scalar_1D_Diffusion_Heat_Moving_ConstantBC_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_1D_Diffusion_Heat_Moving_ConstantBC | 0.0625 | L^2 | 7 | [8] | 1.33804e-06 | 1.33648e-06 | 6.45012e-08 | 0 | NaN | NaN | NaN |
| Scalar_1D_Diffusion_Heat_Moving_ConstantBC | 0.03125 | L^2 | 15 | [16] | 3.85023e-08 | 3.8496e-08 | 6.93673e-10 | 0 | 5.11903 | 5.11759 | 6.53893 |
| Scalar_1D_Diffusion_Heat_Moving_ConstantBC | 0.015625 | L^2 | 31 | [32] | 2.68751e-09 | 2.68745e-09 | 1.7485e-11 | 0 | 3.8406 | 3.8404 | 5.31006 |
| Scalar_1D_Diffusion_Heat_Moving_ConstantBC | 0.0078125 | L^2 | 63 | [64] | 1.12087e-09 | 1.12087e-09 | 2.61509e-12 | 0 | 1.26164 | 1.26162 | 2.74119 |
