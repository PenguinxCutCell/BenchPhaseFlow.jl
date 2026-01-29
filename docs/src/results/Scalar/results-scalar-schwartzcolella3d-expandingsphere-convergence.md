# SchwartzColella3D ExpandingSphere Convergence

3D Schwartz–Colella heat equation on an expanding sphere Ω₃(t) = {(x,y,z): r < 0.392 + t}.
The analytic solution and source term are identical to the fixed case; the radius
grows with time while Dirichlet data track the exact field on the moving boundary
and the outer box.
NOTE : MUST USE IMPLICITINTEGRATION FOR CENTROID SURFACE

**CSV source:** `results/scalar/SchwartzColella3D_ExpandingSphere_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SchwartzColella3D_ExpandingSphere | 0.25 | L^2 | 18 | [4, 4, 4] | 0.0141104 | 0.00755132 | 0.0119198 | 0 | NaN | NaN | NaN |
| SchwartzColella3D_ExpandingSphere | 0.125 | L^2 | 147 | [8, 8, 8] | 0.00697838 | 0.00580394 | 0.00387454 | 0 | 1.0158 | 0.379695 | 1.62126 |
| SchwartzColella3D_ExpandingSphere | 0.0833333 | L^2 | 611 | [12, 12, 12] | 0.00400965 | 0.00362576 | 0.00171206 | 0 | 1.36661 | 1.16033 | 2.0143 |
| SchwartzColella3D_ExpandingSphere | 0.0625 | L^2 | 1479 | [16, 16, 16] | 0.0027139 | 0.0025211 | 0.00100466 | 0 | 1.35676 | 1.26309 | 1.8529 |
