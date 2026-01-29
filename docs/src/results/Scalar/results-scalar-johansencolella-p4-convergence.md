# JohansenColella P4 Convergence

Johansen–Colella Problem 4 (Schwartz Colella 3D Poisson):
Solve ΔΦ = -14 f inside a sphere of radius 0.392 embedded in a unit cube,
with the analytical solution Φ = f, where
    f(x,y,z) = sin(x) sin(2y) sin(3z).
Dirichlet boundary conditions enforce Φ = f on the immersed sphere surface and
Φ = 0 on the outer box. A mesh-convergence sweep records the errors.

**CSV source:** `results/scalar/JohansenColella_P4_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut | ooc_all | ooc_full | ooc_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| JohansenColella_P4 | 0.125 | L^2 | 57 | [7, 7, 7] | 0.0432373 | 0.0214265 | 0.0375549 | 0 | NaN | NaN | NaN | NaN | NaN | NaN |
| JohansenColella_P4 | 0.0625 | L^2 | 739 | [13, 13, 13] | 0.00081244 | 0.000290567 | 0.000758702 | 0 | 5.73387 | 6.20438 | 5.62932 | 5.73387 | 6.20438 | 5.62932 |
| JohansenColella_P4 | 0.03125 | L^2 | 6841 | [26, 26, 26] | 0.000310271 | 0.000103425 | 0.000292525 | 0 | 1.38873 | 1.49028 | 1.37497 | 1.38873 | 1.49028 | 1.37497 |
| JohansenColella_P4 | 0.015625 | L^2 | 60453 | [51, 51, 51] | 0.000105626 | 3.27699e-05 | 0.000100414 | 0 | 1.55456 | 1.65814 | 1.5426 | 1.55456 | 1.65814 | 1.5426 |
