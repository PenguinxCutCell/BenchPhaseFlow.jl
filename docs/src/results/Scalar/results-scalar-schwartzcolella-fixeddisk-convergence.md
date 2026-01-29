# SchwartzColella FixedDisk Convergence

Schwartz–Colella 2D heat equation with a fixed circular domain of radius 0.392.
The exact solution and source term are
    a_exact(x,y,t) = 4 / (5π (t+1)) * exp(-(x² + y²)/(5 (t+1))),
    f(x,y,t) = 4 (x² + y² - 5 (t+1)) / (125 π (t+1)³) * exp(-(x² + y²)/(5 (t+1))).
Dirichlet data on the immersed boundary and initial condition use a_exact.

**CSV source:** `results/scalar/SchwartzColella_FixedDisk_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SchwartzColella_FixedDisk | 0.125 | L^2 | 21 | [7, 7] | 0.113794 | 0.0861327 | 0.074365 | 0 | NaN | NaN | NaN |
| SchwartzColella_FixedDisk | 0.0625 | L^2 | 97 | [13, 13] | 0.0003616 | 0.000304748 | 0.000194635 | 0 | 8.29781 | 8.1428 | 8.57771 |
| SchwartzColella_FixedDisk | 0.03125 | L^2 | 441 | [26, 26] | 3.58712e-05 | 2.78409e-05 | 2.26191e-05 | 0 | 3.3335 | 3.45234 | 3.10516 |
| SchwartzColella_FixedDisk | 0.015625 | L^2 | 1877 | [51, 51] | 1.69919e-05 | 1.4883e-05 | 8.19889e-06 | 0 | 1.07798 | 0.903547 | 1.46404 |
