# Scalar 1D Diffusion Heat Dirichlet Convergence

Scalar 1D Diffusion Heat Equation Convergence Test
This script performs a mesh-convergence study for the 1D heat equation
using the Penguin.jl library. It computes errors and estimated orders of convergence
for different mesh sizes and writes the results to a CSV file.
Analytical solution : u(x,t) = 4/π * Σ (1/(2n+1)) * sin((2n+1) * π * (x - (center - radius)) / (2 * radius)) * exp(-κ * ((2n+1) * π / (2 * radius))^2 * t)
where the sum is over n = 0 to ∞.

**CSV source:** `results/scalar/Scalar_1D_Diffusion_Heat_Dirichlet_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_1D_Diffusion_Heat_Dirichlet | 0.5 | L^2 | 1 | [1] | 0.0286856 | 0.0286856 | 0 | 0 | NaN | NaN | NaN |
| Scalar_1D_Diffusion_Heat_Dirichlet | 0.25 | L^2 | 1 | [2] | 0.105203 | 0.078495 | 0.0700437 | 0 | -1.87477 | -1.45227 | NaN |
| Scalar_1D_Diffusion_Heat_Dirichlet | 0.125 | L^2 | 3 | [4] | 0.00183269 | 0.00180393 | 0.000323434 | 0 | 5.84306 | 5.44339 | 7.75864 |
| Scalar_1D_Diffusion_Heat_Dirichlet | 0.0625 | L^2 | 7 | [8] | 0.00125706 | 0.00125527 | 6.71062e-05 | 0 | 0.543911 | 0.523146 | 2.26895 |
| Scalar_1D_Diffusion_Heat_Dirichlet | 0.03125 | L^2 | 15 | [16] | 0.000138178 | 0.000138149 | 2.8335e-06 | 0 | 3.18545 | 3.18369 | 4.56579 |
