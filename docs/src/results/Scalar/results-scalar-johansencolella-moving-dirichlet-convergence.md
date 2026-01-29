# JohansenColella Moving Dirichlet Convergence

Johansen–Colella moving-ellipse benchmark with Dirichlet data on both the outer
rectangle and all moving elliptical holes. The interface centers translate with
the prescribed velocities from the problem statement.

**CSV source:** `results/scalar/JohansenColella_Moving_Dirichlet_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| JohansenColella_Moving_Dirichlet | 0.5 | L^2 | 0 | [1, 1] | NaN | NaN | NaN | NaN | NaN | NaN | NaN |
| JohansenColella_Moving_Dirichlet | 0.222222 | L^2 | 0 | [2, 2] | 0.00450895 | 0 | 0.00450895 | 0 | NaN | NaN | NaN |
| JohansenColella_Moving_Dirichlet | 0.125 | L^2 | 2 | [4, 4] | 0.00197377 | 0.000841045 | 0.00178561 | 0 | 1.43582 | NaN | 1.60995 |
| JohansenColella_Moving_Dirichlet | 0.0606061 | L^2 | 33 | [7, 7] | 0.000562586 | 0.00040528 | 0.000390194 | 0 | 1.73383 | 1.00849 | 2.10089 |
| JohansenColella_Moving_Dirichlet | 0.031746 | L^2 | 151 | [12, 12] | 0.000159383 | 0.000114594 | 0.000110775 | 0 | 1.95048 | 1.9535 | 1.94725 |
| JohansenColella_Moving_Dirichlet | 0.015748 | L^2 | 719 | [24, 24] | 4.95269e-05 | 3.62004e-05 | 3.38002e-05 | 0 | 1.6672 | 1.64371 | 1.69322 |
