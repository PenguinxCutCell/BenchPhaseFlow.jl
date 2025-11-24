# Scalar 2D Diffusion Heat Moving Convergence

Scalar 2D heat equation with an oscillating circular interface. The radius
oscillates periodically, and the exact solution is imposed on the moving
boundary. This script performs a mesh convergence study and writes a CSV
summary without producing plots or timestamped folders.

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Moving_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_2D_Diffusion_Heat_Moving | 1 | L^2 | 1 | [3, 3] | 1.61454 | 0.231511 | 1.59786 | 0 | NaN | NaN | NaN |
| Scalar_2D_Diffusion_Heat_Moving | 0.5 | L^2 | 13 | [6, 6] | 0.219066 | 0.149636 | 0.159997 | 0 | 2.88169 | 0.629628 | 3.32002 |
| Scalar_2D_Diffusion_Heat_Moving | 0.25 | L^2 | 69 | [11, 11] | 0.0789012 | 0.0705561 | 0.0353161 | 0 | 1.47325 | 1.08461 | 2.17964 |
| Scalar_2D_Diffusion_Heat_Moving | 0.125 | L^2 | 293 | [21, 21] | 0.0289273 | 0.0277954 | 0.00801294 | 0 | 1.44761 | 1.34393 | 2.13992 |
| Scalar_2D_Diffusion_Heat_Moving | 0.0625 | L^2 | 1273 | [42, 42] | 0.0123655 | 0.0120966 | 0.00256476 | 0 | 1.22611 | 1.20024 | 1.64351 |
| Scalar_2D_Diffusion_Heat_Moving | 0.03125 | L^2 | 5217 | [83, 83] | 0.00541135 | 0.00536801 | 0.000683494 | 0 | 1.19226 | 1.17214 | 1.90782 |
