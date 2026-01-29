# Scalar 3D Diffusion Heat Dirichlet Convergence

Scalar 3D Diffusion Heat Equation Convergence Test
This script performs a mesh-convergence study for the 3D heat equation
using the Penguin.jl library. It mirrors the 2D benchmark but for a sphere,
reporting errors and estimated convergence orders across meshes and exporting
results to a CSV file (no plotting or ad-hoc outputs).
The analytical solution used is the classical series for a cooling sphere:
u(r,t) = Tb + (T0 - Tb) * (2R)/(πr) * Σ_{n=1}^∞ (-1)^{n+1} (1/n) sin(nπr/R) exp(-κ (nπ/R)^2 t),
with the r → 0 limit handled analytically.

**CSV source:** `results/scalar/Scalar_3D_Diffusion_Heat_Dirichlet_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_3D_Diffusion_Heat_Dirichlet | 0.5 | L^2 | 7 | [4, 4, 4] | 0.103795 | 0.0825353 | 0.0629384 | 0 | NaN | NaN | NaN |
| Scalar_3D_Diffusion_Heat_Dirichlet | 0.333333 | L^2 | 57 | [6, 6, 6] | 0.0558776 | 0.0527843 | 0.0183335 | 0 | 1.52726 | 1.10247 | 3.042 |
| Scalar_3D_Diffusion_Heat_Dirichlet | 0.25 | L^2 | 147 | [8, 8, 8] | 0.0420347 | 0.0399685 | 0.0130166 | 0 | 0.989525 | 0.966769 | 1.19059 |
| Scalar_3D_Diffusion_Heat_Dirichlet | 0.2 | L^2 | 341 | [10, 10, 10] | 0.0132367 | 0.0125222 | 0.00429002 | 0 | 5.17828 | 5.20108 | 4.97406 |
