# SchwartzColella ExpandingDisk Convergence

Schwartz–Colella 2D heat equation on an expanding disk:
Ω(t) = {(x,y): r < 0.392 + t}, t ≥ 0. Dirichlet data and source follow the
analytic solution specified for the moving-domain test.

**CSV source:** `results/scalar/SchwartzColella_ExpandingDisk_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SchwartzColella_ExpandingDisk | 0.2 | L^2 | 12 | [5, 5] | 0.189453 | 0.15755 | 0.105216 | 0 | NaN | NaN | NaN |
| SchwartzColella_ExpandingDisk | 0.125 | L^2 | 37 | [8, 8] | 0.166914 | 0.147429 | 0.0782631 | 0 | 0.269485 | 0.141268 | 0.629648 |
| SchwartzColella_ExpandingDisk | 0.0625 | L^2 | 161 | [16, 16] | 0.108855 | 0.0944956 | 0.0540375 | 0 | 0.616697 | 0.641701 | 0.534371 |
| SchwartzColella_ExpandingDisk | 0.03125 | L^2 | 725 | [32, 32] | 0.0675463 | 0.0623182 | 0.0260566 | 0 | 0.688462 | 0.600593 | 1.05231 |
| SchwartzColella_ExpandingDisk | 0.015625 | L^2 | 2973 | [63, 63] | 0.0390587 | 0.0345704 | 0.0181789 | 0 | 0.790233 | 0.850118 | 0.519383 |
| SchwartzColella_ExpandingDisk | 0.0078125 | L^2 | 12225 | [126, 126] | 0.0187516 | 0.0176582 | 0.0063095 | 0 | 1.05863 | 0.969196 | 1.52667 |
