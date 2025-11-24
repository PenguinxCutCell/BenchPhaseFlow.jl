# GibouFedkiw Poisson1D Convergence

1D Gibou–Fedkiw Poisson problem on Ω = [-0.5, 0.5] with cut interface φ = |x| - 0.313.
Analytic solution: u(x) = 4 x^2 sin(2π x), f = u''.

**CSV source:** `results/scalar/GibouFedkiw/GibouFedkiw_Poisson1D_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GibouFedkiw_Poisson1D | 0.25 | L^2 | 1 | [3] | 0.370332 | 0.207132 | 0.306989 | 0 | NaN | NaN | NaN |
| GibouFedkiw_Poisson1D | 0.125 | L^2 | 5 | [6] | 0.352631 | 0.352202 | 0.0173794 | 0 | 0.0706609 | -0.765852 | 4.14273 |
| GibouFedkiw_Poisson1D | 0.0625 | L^2 | 9 | [11] | 0.00255331 | 0.00251254 | 0.000454476 | 0 | 7.10964 | 7.13111 | 5.25703 |
| GibouFedkiw_Poisson1D | 0.03125 | L^2 | 19 | [21] | 0.000678326 | 0.000669338 | 0.000110055 | 0 | 1.91232 | 1.90834 | 2.04597 |
| GibouFedkiw_Poisson1D | 0.015625 | L^2 | 39 | [41] | 0.000186914 | 0.000185386 | 2.38493e-05 | 0 | 1.85961 | 1.85221 | 2.20622 |
| GibouFedkiw_Poisson1D | 0.0078125 | L^2 | 79 | [81] | 4.9563e-05 | 4.92911e-05 | 5.18443e-06 | 0 | 1.91503 | 1.91113 | 2.20169 |
| GibouFedkiw_Poisson1D | 0.00390625 | L^2 | 159 | [161] | 1.27916e-05 | 1.27318e-05 | 1.2347e-06 | 0 | 1.95407 | 1.95289 | 2.07002 |
