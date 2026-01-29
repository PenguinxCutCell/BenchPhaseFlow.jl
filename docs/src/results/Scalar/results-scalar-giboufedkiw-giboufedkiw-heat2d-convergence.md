# GibouFedkiw Heat2D Convergence

2D heat equation with T_t = ΔT on Ω- (star-shaped interface). Exact solution:
T(x,y,t) = exp(-2t) sin(x) sin(y).

**CSV source:** `results/scalar/GibouFedkiw/GibouFedkiw_Heat2D_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GibouFedkiw_Heat2D | 0.0487805 | L^2 | 346 | [21, 21] | 0.0138197 | 0.0110308 | 0.00832502 | 0 | NaN | NaN | NaN |
| GibouFedkiw_Heat2D | 0.0246914 | L^2 | 1500 | [41, 41] | 0.0138112 | 0.0125121 | 0.00584776 | 0 | 0.00090245 | -0.185067 | 0.518751 |
| GibouFedkiw_Heat2D | 0.0124224 | L^2 | 6164 | [81, 81] | 0.013819 | 0.0130934 | 0.00441887 | 0 | -0.000814796 | -0.0661024 | 0.407852 |
