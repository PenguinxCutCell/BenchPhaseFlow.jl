# JohansenColella P5 Convergence

Johansen–Colella Problem 5: 3D heat equation inside a sphere embedded in a unit
box with the analytic solution
    Φ(x,y,z,t) = 4 / (5π (t+1)) * exp(-(x² + y² + z²)/(5 (t+1))).
The source term is
    f(x,y,z,t) = 4 * (x² + y² + z² + 5(t+1)) / (125 π (t+1)³)
                 * exp(-(x² + y² + z²)/(5 (t+1))),
and Dirichlet boundary conditions enforce Φ on the immersed sphere. A transient
simulation to `Tend` collects the final-time error across several meshes.

**CSV source:** `results/scalar/JohansenColella_P5_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| JohansenColella_P5 | 0.125 | L^2 | 57 | [7, 7, 7] | 0.0672915 | 0.0379406 | 0.0555757 | 0 | NaN | NaN | NaN |
| JohansenColella_P5 | 0.0833333 | L^2 | 257 | [10, 10, 10] | 0.0183962 | 0.0107385 | 0.0149366 | 0 | 3.19853 | 3.11294 | 3.24055 |
| JohansenColella_P5 | 0.0555556 | L^2 | 1045 | [15, 15, 15] | 0.00181244 | 0.00179404 | 0.000257565 | 0 | 5.71558 | 4.41311 | 10.0139 |
