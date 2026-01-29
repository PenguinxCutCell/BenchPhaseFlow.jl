# GibouFedkiw Poisson3D Convergence

3D Poisson problem on Ω = [0,1]^3 with interior sphere φ = √((x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2) - 0.3.
Analytical solution: u = exp(-x^2 - y^2 - z^2).

**CSV source:** `results/scalar/GibouFedkiw/GibouFedkiw_Poisson3D_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GibouFedkiw_Poisson3D | 0.125 | L^2 | 19 | [5, 5, 5] | 0.00163228 | 0.000648948 | 0.00149773 | 0 | NaN | NaN | NaN |
| GibouFedkiw_Poisson3D | 0.0625 | L^2 | 281 | [10, 10, 10] | 0.000667492 | 0.000221883 | 0.000629534 | 0 | 1.29007 | 1.5483 | 1.25042 |
| GibouFedkiw_Poisson3D | 0.03125 | L^2 | 2855 | [20, 20, 20] | 0.00025447 | 7.25174e-05 | 0.000243919 | 0 | 1.39125 | 1.6134 | 1.36788 |
| GibouFedkiw_Poisson3D | 0.015625 | L^2 | 26097 | [39, 39, 39] | 9.11428e-05 | 2.40396e-05 | 8.79153e-05 | 0 | 1.4813 | 1.59291 | 1.47221 |
