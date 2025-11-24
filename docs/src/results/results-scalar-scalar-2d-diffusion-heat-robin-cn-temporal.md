# Scalar 2D Diffusion Heat Robin CN Temporal

Temporal convergence benchmark for the scalar 2D heat equation with Robin
interface conditions. The geometry, physics, and analytical reference are the
same as `Scalar_2D_Diffusion_Heat_Robin.jl`, but the mesh is fixed while the
time step varies to compare Backward Euler (BE) and Crank–Nicolson (CN).

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Robin_CN_Temporal.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_2D_Diffusion_Heat_Robin_CN | 0.02 | L^2 | 12637 | [128, 128] | 0.027207 | 0.0271227 | 0.00214126 | 0 |
| Scalar_2D_Diffusion_Heat_Robin_CN | 0.01 | L^2 | 12637 | [128, 128] | 0.0270747 | 0.0268967 | 0.00309941 | 0 |
| Scalar_2D_Diffusion_Heat_Robin_CN | 0.005 | L^2 | 12637 | [128, 128] | 0.00689929 | 0.00680531 | 0.00113489 | 0 |
| Scalar_2D_Diffusion_Heat_Robin_CN | 0.0025 | L^2 | 12637 | [128, 128] | 0.00345928 | 0.00342039 | 0.000517234 | 0 |
