# GibouFedkiw Heat3D Convergence

3D heat equation on Ω = [0,0.5]^3 with interior sphere φ = √((x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2) - 0.15.
Exact solution: T = exp(-3t) sin(x) sin(y) sin(z).

**CSV source:** `results/scalar/GibouFedkiw/GibouFedkiw_Heat3D_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GibouFedkiw_Heat3D | 0.0833333 | L^2 | 4 | [4, 4, 4] | 0.0270026 | 0.0210504 | 0.0169122 | 0 | NaN | NaN | NaN |
| GibouFedkiw_Heat3D | 0.05 | L^2 | 17 | [6, 6, 6] | 0.0154584 | 0.0135272 | 0.00748165 | 0 | 1.09192 | 0.865682 | 1.5966 |
| GibouFedkiw_Heat3D | 0.0357143 | L^2 | 38 | [9, 9, 9] | 0.0106889 | 0.0092538 | 0.00534967 | 0 | 1.09652 | 1.12838 | 0.996869 |
