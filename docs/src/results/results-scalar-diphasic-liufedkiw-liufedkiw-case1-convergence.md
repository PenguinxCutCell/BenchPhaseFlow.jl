# LiuFedkiw Case1 Convergence

Liu-Fedkiw diphasic diffusion benchmark (Case 1).

Problem: u_xx = 0 on [0, 1] with u(0) = 0, u(1) = 2, interface at x = 0.5, jump
conditions [u] = 1 and [u_x] = 0. Exact solution is u = x on the left of the
interface and u = x + 1 on the right.

**CSV source:** `results/scalar/diphasic/LiuFedkiw/LiuFedkiw_Case1_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_phase1 | inside_cells_phase2 | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | phase1_all_err | phase1_full_err | phase1_cut_err | phase1_empty_err | phase2_all_err | phase2_full_err | phase2_cut_err | phase2_empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LiuFedkiw_Case1 | 0.125 | L^2 | 7 | 3 | 4 | [3, 4] | 0.104138 | 0.100109 | 0.0286848 | 0 | 0.104138 | 0.100109 | 0.0286848 | 0 | 0.0386045 | 0.031497 | 0.0223214 | 0 | NaN | NaN | NaN |
| LiuFedkiw_Case1 | 0.0625 | L^2 | 15 | 7 | 8 | [7, 8] | 0.0497687 | 0.0489709 | 0.00887559 | 0 | 0.0497687 | 0.0489709 | 0.00887559 | 0 | 0.0186356 | 0.01691 | 0.00783188 | 0 | 1.06518 | 1.03158 | 1.69237 |
| LiuFedkiw_Case1 | 0.03125 | L^2 | 31 | 15 | 16 | [15, 16] | 0.024361 | 0.0241827 | 0.00294212 | 0 | 0.024361 | 0.0241827 | 0.00294212 | 0 | 0.00916555 | 0.00873891 | 0.00276383 | 0 | 1.03066 | 1.01795 | 1.59298 |
| LiuFedkiw_Case1 | 0.015625 | L^2 | 63 | 31 | 32 | [31, 32] | 0.0120554 | 0.0120132 | 0.00100772 | 0 | 0.0120554 | 0.0120132 | 0.00100772 | 0 | 0.00454621 | 0.00444005 | 0.000976712 | 0 | 1.0149 | 1.00936 | 1.54577 |
| LiuFedkiw_Case1 | 0.0078125 | L^2 | 127 | 63 | 64 | [63, 64] | 0.00599708 | 0.00598681 | 0.000350718 | 0 | 0.00599708 | 0.00598681 | 0.000350718 | 0 | 0.00226414 | 0.00223765 | 0.00034528 | 0 | 1.00735 | 1.00476 | 1.52271 |
| LiuFedkiw_Case1 | 0.00390625 | L^2 | 255 | 127 | 128 | [127, 128] | 0.00299097 | 0.00298843 | 0.000123029 | 0 | 0.00299097 | 0.00298843 | 0.000123029 | 0 | 0.00112985 | 0.00112323 | 0.000122071 | 0 | 1.00365 | 1.0024 | 1.51131 |
