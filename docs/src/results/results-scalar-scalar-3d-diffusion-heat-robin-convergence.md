# Scalar 3D Diffusion Heat Robin Convergence

Scalar 3D Diffusion Heat Equation with Robin Boundary Conditions
This benchmark mirrors the 2D Robin problem but inside a sphere. The initial
temperature is uniform (w0) and the interface satisfies ∂ₙw + k w = 0. The
analytical solution uses the eigenvalues μₙ obtained from μ cot(μ) + kR - 1 = 0:
w(r,t) = (2kR² w0)/r * Σ Cₙ sin(μₙ r / R) exp(-a μₙ² t / R²),
where Cₙ = sin(μₙ)[μₙ²+(kR-1)²]/(μₙ²[μₙ²+kR(kR-1)]). Th

**CSV source:** `results/scalar/Scalar_3D_Diffusion_Heat_Robin_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_3D_Diffusion_Heat_Robin | 1 | L^2 | 1 | [2, 2, 2] | 0.372861 | 0.314916 | 0.199634 | 0 | NaN | NaN | NaN |
| Scalar_3D_Diffusion_Heat_Robin | 0.5 | L^2 | 7 | [4, 4, 4] | 0.153828 | 0.0714126 | 0.136247 | 0 | 1.27732 | 2.14071 | 0.551134 |
| Scalar_3D_Diffusion_Heat_Robin | 0.25 | L^2 | 147 | [8, 8, 8] | 0.048812 | 0.0358195 | 0.0331599 | 0 | 1.65601 | 0.995435 | 2.03872 |
| Scalar_3D_Diffusion_Heat_Robin | 0.125 | L^2 | 1599 | [16, 16, 16] | 0.011079 | 0.00945192 | 0.00577969 | 0 | 2.13941 | 1.92207 | 2.52038 |
| Scalar_3D_Diffusion_Heat_Robin | 0.0625 | L^2 | 14915 | [32, 32, 32] | 0.00317916 | 0.00289967 | 0.00130345 | 0 | 1.8011 | 1.70472 | 2.14866 |
