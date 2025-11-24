# Scalar 2D Diffusion Heat Dirichlet Convergence

Scalar 2D Diffusion Heat Equation Convergence Test
This script performs a mesh-convergence study for the 2D heat equation with Dirichlet boundary conditions
using the Penguin.jl library. It computes errors and estimated orders of convergence
for different mesh sizes and writes the results to a CSV file.
Analytical solution : u(r,t) = (1 - Σ A_n * exp(-a * (α_n^2) * t / R^2) * J0(α_n * (r / R))) * (wr - w0) + w0
where A_n are coefficients based on Robin boundary conditions and α_n are the roots of the Bessel function equation.

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Dirichlet_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_2D_Diffusion_Heat_Dirichlet | 1 | L^2 | 1 | [2, 2] | 0.46687 | 0.150081 | 0.44209 | 0 | NaN | NaN | NaN |
| Scalar_2D_Diffusion_Heat_Dirichlet | 0.5 | L^2 | 5 | [4, 4] | 0.233409 | 0.18988 | 0.13574 | 0 | 1.00016 | -0.33935 | 1.7035 |
| Scalar_2D_Diffusion_Heat_Dirichlet | 0.25 | L^2 | 37 | [8, 8] | 0.12269 | 0.121004 | 0.020275 | 0 | 0.927837 | 0.650039 | 2.74307 |
| Scalar_2D_Diffusion_Heat_Dirichlet | 0.125 | L^2 | 173 | [16, 16] | 0.0233681 | 0.0232398 | 0.00244536 | 0 | 2.39241 | 2.38038 | 3.05158 |
| Scalar_2D_Diffusion_Heat_Dirichlet | 0.0625 | L^2 | 741 | [32, 32] | 0.00923052 | 0.00921257 | 0.000575305 | 0 | 1.34006 | 1.33492 | 2.08765 |
