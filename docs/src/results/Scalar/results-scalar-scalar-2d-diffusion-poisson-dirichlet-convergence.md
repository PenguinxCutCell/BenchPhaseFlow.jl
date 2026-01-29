# Scalar 2D Diffusion Poisson Dirichlet Convergence

Scalar 2D Diffusion Poisson Problem with Dirichlet Boundary Conditions
This script performs a mesh-convergence study for the 2D Poisson equation
using the Penguin.jl library. It computes errors and estimated orders of convergence
for different mesh sizes and writes the results to a CSV file.
Analytical solution is : u(x,y) = 1 - (x-center_x)^2 - (y-center_y)^2

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Poisson_Dirichlet_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| FrontTracking | 1 | L^2 | 1 | [2, 2] | 0.261399 | 0.0342789 | 0.259142 | 0 |
| FrontTracking | 0.5 | L^2 | 5 | [4, 4] | 0.035029 | 0.0282016 | 0.0207774 | 0 |
| FrontTracking | 0.25 | L^2 | 37 | [8, 8] | 0.00843131 | 0.00790037 | 0.00294467 | 0 |
| FrontTracking | 0.125 | L^2 | 177 | [16, 16] | 0.00237533 | 0.00227285 | 0.000690185 | 0 |
| FrontTracking | 0.0625 | L^2 | 749 | [32, 32] | 0.000705592 | 0.000520274 | 0.000476629 | 0 |
| FrontTracking | 0.03125 | L^2 | 3101 | [64, 64] | 0.000186959 | 0.000141025 | 0.000122742 | 0 |
| ImplicitIntegration | 1 | L^2 | 1 | [2, 2] | 0.261391 | 0.0342685 | 0.259135 | 0 |
| ImplicitIntegration | 0.5 | L^2 | 5 | [4, 4] | 0.0350525 | 0.0282195 | 0.0207928 | 0 |
| ImplicitIntegration | 0.25 | L^2 | 37 | [8, 8] | 0.00844843 | 0.00792336 | 0.00293195 | 0 |
| ImplicitIntegration | 0.125 | L^2 | 177 | [16, 16] | 0.00239657 | 0.00229729 | 0.000682623 | 0 |
| ImplicitIntegration | 0.0625 | L^2 | 749 | [32, 32] | 0.000714679 | 0.000535452 | 0.000473347 | 0 |
| ImplicitIntegration | 0.03125 | L^2 | 3101 | [64, 64] | 0.000189584 | 0.000144919 | 0.000122231 | 0 |
