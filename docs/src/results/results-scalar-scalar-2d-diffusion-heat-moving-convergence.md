# Scalar 2D Diffusion Heat Moving Convergence

Scalar 2D heat equation with an oscillating circular interface. The radius
oscillates periodically, and the exact solution is imposed on the moving
boundary. This script performs a mesh convergence study and writes a CSV
summary without producing plots or timestamped folders.
# Might need to adjust interface centroid computation for moving bodies : bary_interface vs compute_interface_centroid()

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Moving_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_2D_Diffusion_Heat_Moving | 1 | L^2 | 5 | [4, 4] | 1.71959 | 1.5846 | 0.667873 | 0 | NaN | NaN | NaN |
| Scalar_2D_Diffusion_Heat_Moving | 0.5 | L^2 | 21 | [7, 7] | 0.328935 | 0.224704 | 0.240221 | 0 | 2.38619 | 2.81802 | 1.47521 |
| Scalar_2D_Diffusion_Heat_Moving | 0.25 | L^2 | 97 | [13, 13] | 0.157871 | 0.147073 | 0.057382 | 0 | 1.05905 | 0.611492 | 2.06569 |
| Scalar_2D_Diffusion_Heat_Moving | 0.125 | L^2 | 449 | [26, 26] | 0.043205 | 0.0404496 | 0.0151824 | 0 | 1.86948 | 1.86234 | 1.9182 |
| Scalar_2D_Diffusion_Heat_Moving | 0.0625 | L^2 | 1925 | [51, 51] | 0.0166473 | 0.0164591 | 0.00249611 | 0 | 1.37591 | 1.29724 | 2.60464 |
| Scalar_2D_Diffusion_Heat_Moving | 0.03125 | L^2 | 7909 | [102, 102] | 0.00570675 | 0.00568447 | 0.000503859 | 0 | 1.54455 | 1.53379 | 2.30859 |
