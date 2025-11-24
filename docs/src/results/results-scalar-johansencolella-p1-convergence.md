# JohansenColella P1 Convergence

Johansen–Colella Problem 1: constant-coefficient Poisson equation inside a
star-shaped domain defined by
    Ω = {(r, θ) : r ≤ 0.30 + 0.15 cos(6θ)}.
We solve Δϕ = 7 r² cos(3θ) with Dirichlet data taken from the exact solution
ϕ(r, θ) = r⁴ cos(3θ). Convergence is measured on uniform Cartesian meshes using
Penguin's implicit-integration capacity.

**CSV source:** `results/scalar/JohansenColella_P1_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| JohansenColella_P1 | 0.0625 | L^2 | 39 | [15, 15] | 0.00372574 | 0.00123643 | 0.0035146 | 0 | NaN | NaN | NaN |
| JohansenColella_P1 | 0.03125 | L^2 | 228 | [29, 29] | 0.00149813 | 0.00106683 | 0.0010518 | 0 | 1.31436 | 0.212856 | 1.7405 |
| JohansenColella_P1 | 0.015625 | L^2 | 1138 | [58, 58] | 3.09766e-05 | 1.50406e-05 | 2.70801e-05 | 0 | 5.59584 | 6.14833 | 5.27948 |
| JohansenColella_P1 | 0.0078125 | L^2 | 4868 | [116, 116] | 1.13903e-05 | 4.82107e-06 | 1.03196e-05 | 0 | 1.44338 | 1.64143 | 1.39184 |
| JohansenColella_P1 | 0.00390625 | L^2 | 20166 | [231, 231] | 3.9584e-06 | 1.3986e-06 | 3.70309e-06 | 0 | 1.52481 | 1.78537 | 1.47859 |
