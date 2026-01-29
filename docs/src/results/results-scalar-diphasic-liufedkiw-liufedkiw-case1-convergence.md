# LiuFedkiw Case1 Convergence

Liu-Fedkiw diphasic diffusion benchmark (Case 1).

Problem: u_xx = 0 on [0, 1] with u(0) = 0, u(1) = 2, interface at x = 0.5, jump
conditions [u] = 1 and [u_x] = 0. For convergence, the analytical solution is
evaluated on the effective domain using L_eff = L - dx:
u_left(x) = (x - dx) / L_eff and u_right(x) = u_left(x) + 1.

**CSV source:** `results/scalar/diphasic/LiuFedkiw/LiuFedkiw_Case1_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_phase1 | inside_cells_phase2 | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | phase1_all_err | phase1_full_err | phase1_cut_err | phase1_empty_err | phase2_all_err | phase2_full_err | phase2_cut_err | phase2_empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LiuFedkiw_Case1 | 0.125 | L^2 | 7 | 3 | 4 | [3, 4] | 3.62597e-16 | 3.31005e-16 | 2.09812e-16 | 0 | 3.12965e-16 | 2.3222e-16 | 2.09812e-16 | 0 | 3.62597e-16 | 3.31005e-16 | 1.4803e-16 | 0 | NaN | NaN | NaN |
| LiuFedkiw_Case1 | 0.0625 | L^2 | 15 | 7 | 8 | [7, 8] | 2.01502e-16 | 1.703e-16 | 1.07707e-16 | 0 | 1.75542e-16 | 1.53034e-16 | 8.59975e-17 | 0 | 2.01502e-16 | 1.703e-16 | 1.07707e-16 | 0 | 0.847573 | 0.958769 | 0.961982 |
| LiuFedkiw_Case1 | 0.03125 | L^2 | 31 | 15 | 16 | [15, 16] | 2.73291e-15 | 2.63238e-15 | 7.77668e-16 | 0 | 2.57566e-15 | 2.45546e-15 | 7.77668e-16 | 0 | 2.73291e-15 | 2.63238e-15 | 7.34407e-16 | 0 | -3.76157 | -3.95021 | -2.85204 |
| LiuFedkiw_Case1 | 0.015625 | L^2 | 63 | 31 | 32 | [31, 32] | 5.69951e-15 | 5.67283e-15 | 5.50825e-16 | 0 | 3.05489e-15 | 3.00578e-15 | 5.45512e-16 | 0 | 5.69951e-15 | 5.67283e-15 | 5.50825e-16 | 0 | -1.0604 | -1.1077 | 0.49756 |
| LiuFedkiw_Case1 | 0.0078125 | L^2 | 127 | 63 | 64 | [63, 64] | 1.16277e-14 | 1.16276e-14 | 5.91098e-17 | 0 | 1.52013e-15 | 1.51898e-15 | 5.91098e-17 | 0 | 1.16277e-14 | 1.16276e-14 | 5.86498e-17 | 0 | -1.02866 | -1.03541 | 3.22012 |
| LiuFedkiw_Case1 | 0.00390625 | L^2 | 255 | 127 | 128 | [127, 128] | 4.16865e-14 | 4.1606e-14 | 2.59009e-15 | 0 | 2.5917e-14 | 2.57876e-14 | 2.58632e-15 | 0 | 4.16865e-14 | 4.1606e-14 | 2.59009e-15 | 0 | -1.84201 | -1.83924 | -5.45346 |
