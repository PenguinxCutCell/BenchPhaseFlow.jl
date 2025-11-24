# GibouFedkiw Poisson2D Convergence

2D Poisson problem from Gibou & Fedkiw on Ω = [-1,1]² with a star-shaped cut.
Exact solution: u(x,y) = x² + y². The RHS satisfies Δu = 4.
Dirichlet data equal the analytic solution are enforced on the interface and outer box.

**CSV source:** `results/scalar/GibouFedkiw/GibouFedkiw_Poisson2D_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GibouFedkiw_Poisson2D | 0.125 | L^2 | 34 | [8, 8] | 0.217817 | 0.122243 | 0.18028 | 0 | NaN | NaN | NaN |
| GibouFedkiw_Poisson2D | 0.0625 | L^2 | 203 | [16, 16] | 0.222428 | 0.171438 | 0.141715 | 0 | -0.0302184 | -0.487929 | 0.347247 |
| GibouFedkiw_Poisson2D | 0.03125 | L^2 | 912 | [32, 32] | 0.223709 | 0.19453 | 0.110471 | 0 | -0.00828994 | -0.182311 | 0.359327 |
