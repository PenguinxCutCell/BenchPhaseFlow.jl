# Scalar 1D Diffusion Poisson Dirichlet Convergence

Scalar 1D Diffusion Poisson Problem with Dirichlet Boundary Conditions
This script performs a mesh-convergence study for the 1D Poisson equation
using the Penguin.jl library. It computes errors and estimated orders of convergence
for different mesh sizes and writes the results to a CSV file.
Analytical solution is : u(x) = - (x-center)^3/6 - (center*(x-center)^2)/2 + radius^2/6*(x-center) + center*radius^2/2

**CSV source:** `results/scalar/Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 0.5 | L^2 | 0 | [1] | 0.003025 | 0 | 0.003025 | 0 | NaN | NaN | NaN |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 0.25 | L^2 | 0 | [1] | 0.003025 | 0 | 0.003025 | 0 | 3.20343e-16 | NaN | 3.20343e-16 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 0.125 | L^2 | 1 | [2] | 0.000741946 | 0.000736112 | 9.28664e-05 | 0 | 2.02755 | NaN | 5.02564 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 0.0625 | L^2 | 3 | [4] | 0.000226443 | 0.000226353 | 6.36181e-06 | 0 | 1.71217 | 1.70135 | 3.86765 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 0.03125 | L^2 | 7 | [8] | 6.13316e-05 | 6.13316e-05 | 2.37516e-09 | 0 | 1.88444 | 1.88387 | 11.3872 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 0.015625 | L^2 | 13 | [15] | 1.48081e-05 | 1.47547e-05 | 1.25646e-06 | 0 | 2.05024 | 2.05545 | -9.04712 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 0.0078125 | L^2 | 27 | [29] | 3.77127e-06 | 3.76187e-06 | 2.66121e-07 | 0 | 1.97326 | 1.97165 | 2.23921 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 0.00390625 | L^2 | 55 | [57] | 9.51752e-07 | 9.49526e-07 | 6.50546e-08 | 0 | 1.98639 | 1.98617 | 2.03236 |
