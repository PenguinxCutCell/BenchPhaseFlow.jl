# JohansenColella P4 Convergence

Johansen–Colella Problem 4 (Schwartz Colella 3D Poisson):
Solve ΔΦ = -14 f inside a sphere of radius 0.392 embedded in a unit cube,
with the analytical solution Φ = f, where
    f(x,y,z) = sin(x) sin(2y) sin(3z).
Dirichlet boundary conditions enforce Φ = f on the immersed sphere surface and
Φ = 0 on the outer box. A mesh-convergence sweep records the errors.

**CSV source:** `results/scalar/JohansenColella_P4_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| JohansenColella_P4 | 0.125 | L^2 | 57 | [7, 7, 7] | 0.0429336 | 0.0210024 | 0.0374459 | 0 |
| JohansenColella_P4 | 0.0625 | L^2 | 739 | [13, 13, 13] | 0.000934657 | 0.000622896 | 0.00069684 | 0 |
| JohansenColella_P4 | 0.03125 | L^2 | 6841 | [26, 26, 26] | 0.0003118 | 0.000170552 | 0.00026102 | 0 |
| JohansenColella_P4 | 0.015625 | L^2 | 60453 | [51, 51, 51] | 0.000105052 | 4.88597e-05 | 9.2998e-05 | 0 |
