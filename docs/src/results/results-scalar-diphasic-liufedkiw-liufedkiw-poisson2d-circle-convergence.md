# LiuFedkiw Poisson2D Circle Convergence

Liu–Fedkiw diphasic Poisson benchmark in 2D with a circular interface.

Domain: [0, 1] x [0, 1]. Interface: circle centered at (0.5, 0.5) with radius
0.25. Interior exact solution u1(x, y) = exp(-x^2 - y^2); exterior u2(x, y) = 0.
Diffusivity: D1 = 2 (interior), D2 = 1 (exterior).

Source term:
    f1(x, y) = 8(x^2 + y^2 - 1)exp(-x^2 - y^2)  (interior)
    f2(x, y) = 0                                (exterior)

Interface jumps (phase 2 minus phase 1):
    [u] = -exp(-x^2 - y^2)
    [D∇u·n] = 8(2x^2 + 2y^2 - x - y)exp(-x^2 - y^2)
Outer boundaries: Dirichlet 0.

**CSV source:** `results/scalar/diphasic/LiuFedkiw/LiuFedkiw_Poisson2D_Circle_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_phase1 | inside_cells_phase2 | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | phase1_all_err | phase1_full_err | phase1_cut_err | phase1_empty_err | phase2_all_err | phase2_full_err | phase2_cut_err | phase2_empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LiuFedkiw_Poisson2D_Circle | 0.0625 | L^2 | 224 | 187 | 37 | [187, 37] | 0.975889 | 0.837265 | 0.501343 | 0 | 0.549405 | 0.521011 | 0.174336 | 0 | 0.975889 | 0.837265 | 0.501343 | 0 | NaN | NaN | NaN |
| LiuFedkiw_Poisson2D_Circle | 0.03125 | L^2 | 960 | 783 | 177 | [783, 177] | 0.969025 | 0.909188 | 0.335241 | 0 | 0.560522 | 0.546258 | 0.125648 | 0 | 0.969025 | 0.909188 | 0.335241 | 0 | 0.0101824 | -0.118894 | 0.5806 |
| LiuFedkiw_Poisson2D_Circle | 0.015625 | L^2 | 3968 | 3219 | 749 | [3219, 749] | 0.965507 | 0.931751 | 0.253068 | 0 | 0.565534 | 0.559224 | 0.0842473 | 0 | 0.965507 | 0.931751 | 0.253068 | 0 | 0.00524799 | -0.0353654 | 0.405677 |
| LiuFedkiw_Poisson2D_Circle | 0.0078125 | L^2 | 16128 | 13027 | 3101 | [13027, 3101] | 0.963731 | 0.946196 | 0.183005 | 0 | 0.567902 | 0.564907 | 0.0582446 | 0 | 0.963731 | 0.946196 | 0.183005 | 0 | 0.00265523 | -0.0221949 | 0.467641 |
