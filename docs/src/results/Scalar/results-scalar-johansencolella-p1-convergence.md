# JohansenColella P1 Convergence

Johansen–Colella Problem 1: constant-coefficient Poisson equation inside a
star-shaped domain defined by
    Ω = {(r, θ) : r ≤ 0.30 + 0.15 cos(6θ)}.
We solve Δϕ = 7 r² cos(3θ) with Dirichlet data taken from the exact solution
ϕ(r, θ) = r⁴ cos(3θ). Convergence is measured on uniform Cartesian meshes using
Penguin's implicit-integration capacity.

Truncation errors are computed by inserting the exact solution sampled at the cut-cell centroids into the discrete operator and comparing it with the continuous operator evaluated at the same centroid locations. 
This differs slightly from the convention used by Johansen and Colella, who evaluate the exact solution at Cartesian cell centers while evaluating the right-hand side at cut-cell centroids. 
The error constants might differ.

**CSV source:** `results/scalar/JohansenColella_P1_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut | trunc_max_all | trunc_max_full | trunc_max_cut | trunc_pair_order_max_all | trunc_pair_order_max_full | trunc_pair_order_max_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| JohansenColella_P1 | 0.025 | L^1 | 411 | [36, 36] | 4.44838e-05 | 2.427e-05 | 2.02138e-05 | 0 | NaN | NaN | NaN | 0.338488 | 0.0940073 | 0.338488 | NaN | NaN | NaN |
| JohansenColella_P1 | 0.0125 | L^1 | 1825 | [72, 72] | 1.33866e-05 | 8.38378e-06 | 5.00286e-06 | 0 | 1.73249 | 1.5335 | 2.01452 | 0.609185 | 0.0394961 | 0.609185 | -0.847774 | 1.25106 | -0.847774 |
| JohansenColella_P1 | 0.00625 | L^1 | 7723 | [144, 144] | 3.51312e-06 | 2.4156e-06 | 1.09752e-06 | 0 | 1.92997 | 1.79522 | 2.18851 | 0.897194 | 0.0320099 | 0.897194 | -0.55854 | 0.303192 | -0.55854 |
| JohansenColella_P1 | 0.003125 | L^1 | 31691 | [288, 288] | 9.60925e-07 | 6.70148e-07 | 2.90776e-07 | 0 | 1.87026 | 1.84983 | 1.91626 | 0.392498 | 0.0345331 | 0.392498 | 1.19273 | -0.109461 | 1.19273 |
| JohansenColella_P1 | 0.0015625 | L^1 | 128577 | [576, 576] | 2.65409e-07 | 1.94441e-07 | 7.09671e-08 | 0 | 1.85621 | 1.78515 | 2.03469 | 1.92245 | 0.0333296 | 1.92245 | -2.29219 | 0.0511771 | -2.29219 |
