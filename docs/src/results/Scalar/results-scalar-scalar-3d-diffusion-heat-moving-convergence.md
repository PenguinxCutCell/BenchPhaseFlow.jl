# Scalar 3D Diffusion Heat Moving Convergence

Scalar 3D heat equation with an oscillating spherical interface. The radius
oscillates periodically, and the exact solution is imposed on the moving
boundary. This script performs a mesh convergence study and writes a CSV
summary without producing plots or timestamped folders.

**CSV source:** `results/scalar/Scalar_3D_Diffusion_Heat_Moving_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_3D_Diffusion_Heat_Moving | 1 | L^2 | 1 | [3, 3, 3] | 2.50159 | 0.927856 | 2.32315 | 0 | NaN | NaN | NaN |
| Scalar_3D_Diffusion_Heat_Moving | 0.5 | L^2 | 19 | [6, 6, 6] | 0.369591 | 0.116326 | 0.350807 | 0 | 2.75884 | 2.99573 | 2.72733 |
| Scalar_3D_Diffusion_Heat_Moving | 0.333333 | L^2 | 147 | [8, 8, 8] | 0.224538 | 0.195606 | 0.110254 | 0 | 1.22908 | -1.28174 | 2.85462 |
| Scalar_3D_Diffusion_Heat_Moving | 0.25 | L^2 | 389 | [11, 11, 11] | 0.180547 | 0.162087 | 0.0795309 | 0 | 0.757965 | 0.653381 | 1.13543 |
