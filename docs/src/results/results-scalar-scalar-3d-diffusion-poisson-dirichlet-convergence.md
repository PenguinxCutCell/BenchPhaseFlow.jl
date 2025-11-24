# Scalar 3D Diffusion Poisson Dirichlet Convergence

Scalar 3D Diffusion Poisson Problem with Dirichlet Boundary Conditions
This script performs a mesh-convergence study for the 3D Poisson equation
using the Penguin.jl library. It computes errors and estimated orders of convergence
for different mesh sizes and writes the results to a CSV file.
Analytical solution is : u(x,y,z) = 1/6 * (radius^2 - ((x - center_x)^2 + (y - center_y)^2 + (z - center_z)^2))

**CSV source:** `results/scalar/Scalar_3D_Diffusion_Poisson_Dirichlet_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_3D_Diffusion_Poisson_Dirichlet_Convergence | 0.5 | L^2 | 1 | [2, 2, 2] | 0.0117991 | 0.00590999 | 0.0102123 | 0 |
| Scalar_3D_Diffusion_Poisson_Dirichlet_Convergence | 0.25 | L^2 | 7 | [4, 4, 4] | 0.000713823 | 0.000457723 | 0.000547753 | 0 |
| Scalar_3D_Diffusion_Poisson_Dirichlet_Convergence | 0.125 | L^2 | 147 | [8, 8, 8] | 0.000272478 | 0.000157168 | 0.000222582 | 0 |
| Scalar_3D_Diffusion_Poisson_Dirichlet_Convergence | 0.0625 | L^2 | 1599 | [16, 16, 16] | 6.90082e-05 | 4.38654e-05 | 5.32725e-05 | 0 |
