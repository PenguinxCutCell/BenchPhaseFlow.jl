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
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 0.00195312 | L^2 | 111 | [113] | 2.39332e-07 | 2.38512e-07 | 1.97977e-08 | 0 | 1.99157 | 1.99315 | 1.71632 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 0.000976562 | L^2 | 225 | [226] | 6.00442e-08 | 6.00442e-08 | 4.21303e-11 | 0 | 1.99491 | 1.98996 | 8.87626 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 0.000488281 | L^2 | 449 | [451] | 1.50038e-08 | 1.49938e-08 | 5.46147e-10 | 0 | 2.0007 | 2.00166 | -3.69636 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 0.000244141 | L^2 | 901 | [902] | 3.75496e-09 | 3.75496e-09 | 1.58372e-13 | 0 | 1.99846 | 1.9975 | 11.7518 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 0.00012207 | L^2 | 1801 | [1803] | 9.38518e-10 | 9.38469e-10 | 9.61473e-12 | 0 | 2.00034 | 2.00042 | -5.92386 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 6.10352e-05 | L^2 | 3603 | [3605] | 2.34663e-10 | 2.34648e-10 | 2.64527e-12 | 0 | 1.99979 | 1.99981 | 1.86183 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 3.05176e-05 | L^2 | 7207 | [7209] | 5.86693e-11 | 5.86617e-11 | 9.4381e-13 | 0 | 1.99991 | 2.00001 | 1.48685 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 1.52588e-05 | L^2 | 14417 | [14418] | 1.46973e-11 | 1.46973e-11 | 2.51849e-14 | 0 | 1.99706 | 1.99687 | 5.22787 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 7.62939e-06 | L^2 | 28835 | [28836] | 3.5264e-12 | 3.5264e-12 | 3.54643e-15 | 0 | 2.05928 | 2.05928 | 2.82812 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 3.8147e-06 | L^2 | 57671 | [57672] | 8.69407e-13 | 8.69407e-13 | 3.69646e-16 | 0 | 2.02009 | 2.02009 | 3.26215 |
| Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence | 1.90735e-06 | L^2 | 115343 | [115344] | 5.73444e-13 | 5.73444e-13 | 1.33204e-17 | 0 | 0.60038 | 0.600379 | 4.79443 |
