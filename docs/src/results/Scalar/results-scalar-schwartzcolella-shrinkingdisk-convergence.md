# SchwartzColella ShrinkingDisk Convergence

Schwartz–Colella 2D heat equation on a shrinking disk:
Ω(t) = {(x,y): r < 0.392 - t}, 0 ≤ t < 0.392. The exact solution and source are
those specified for the Schwartz–Colella problems, and Dirichlet conditions
match the analytic solution on the moving boundary.

**CSV source:** `results/scalar/SchwartzColella_ShrinkingDisk_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SchwartzColella_ShrinkingDisk | 0.25 | L^2 | 1 | [3, 3] | 0.165268 | 0.0461676 | 0.158689 | 0 | NaN | NaN | NaN |
| SchwartzColella_ShrinkingDisk | 0.125 | L^2 | 13 | [6, 6] | 0.0563622 | 0.0367367 | 0.0427447 | 0 | 1.55201 | 0.329657 | 1.89238 |
| SchwartzColella_ShrinkingDisk | 0.0625 | L^2 | 69 | [11, 11] | 0.000299717 | 0.000188984 | 0.000232628 | 0 | 7.55498 | 7.60282 | 7.52158 |
| SchwartzColella_ShrinkingDisk | 0.03125 | L^2 | 333 | [22, 22] | 0.000137889 | 0.000115219 | 7.57497e-05 | 0 | 1.12009 | 0.713881 | 1.61871 |
| SchwartzColella_ShrinkingDisk | 0.015625 | L^2 | 1421 | [44, 44] | 5.22006e-05 | 5.1448e-05 | 8.83194e-06 | 0 | 1.40137 | 1.16319 | 3.10044 |
| SchwartzColella_ShrinkingDisk | 0.0078125 | L^2 | 5853 | [88, 88] | 2.69497e-05 | 2.67688e-05 | 3.1171e-06 | 0 | 0.953796 | 0.942561 | 1.50253 |
