# JohansenColella Fixed Dirichlet Convergence

Fixed-domain Johansen–Colella problem on Ω = [-1.5,1.5]×[-1,1] minus three
elliptical cavities. Dirichlet data equal to the analytic solution are imposed
on both the outer box and the ellipse boundaries.

**CSV source:** `results/scalar/JohansenColella_Fixed_Dirichlet_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| JohansenColella_Fixed_Dirichlet | 0.25 | L^2 | 49 | [2, 2] | 0.0215141 | 0.0194584 | 0.00917758 | 0 | NaN | NaN | NaN |
| JohansenColella_Fixed_Dirichlet | 0.125 | L^2 | 228 | [4, 4] | 0.00957428 | 0.00950943 | 0.00111242 | 0 | 1.16805 | 1.03296 | 3.04442 |
| JohansenColella_Fixed_Dirichlet | 0.0625 | L^2 | 944 | [7, 7] | 0.00468286 | 0.00467965 | 0.000173355 | 0 | 1.03177 | 1.02296 | 2.6819 |
| JohansenColella_Fixed_Dirichlet | 0.03125 | L^2 | 3824 | [13, 13] | 0.00296725 | 0.00296698 | 4.00769e-05 | 0 | 0.658263 | 0.657406 | 2.11288 |
