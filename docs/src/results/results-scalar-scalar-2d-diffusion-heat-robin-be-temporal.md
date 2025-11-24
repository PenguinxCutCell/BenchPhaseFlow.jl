# Scalar 2D Diffusion Heat Robin BE Temporal

Temporal convergence benchmark for the scalar 2D heat equation with Robin
interface conditions. The geometry, physics, and analytical reference are the
same as `Scalar_2D_Diffusion_Heat_Robin.jl`, but the mesh is fixed while the
time step varies to compare Backward Euler (BE) and Crank–Nicolson (CN).

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Robin_BE_Temporal.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_2D_Diffusion_Heat_Robin_BE | 0.02 | L^2 | 12637 | [128, 128] | 0.023946 | 0.0237975 | 0.0026634 | 0 |
| Scalar_2D_Diffusion_Heat_Robin_BE | 0.01 | L^2 | 12637 | [128, 128] | 0.0253701 | 0.0251742 | 0.00314626 | 0 |
| Scalar_2D_Diffusion_Heat_Robin_BE | 0.005 | L^2 | 12637 | [128, 128] | 0.00612242 | 0.00608257 | 0.000697416 | 0 |
| Scalar_2D_Diffusion_Heat_Robin_BE | 0.0025 | L^2 | 12637 | [128, 128] | 0.00307502 | 0.0030535 | 0.000363223 | 0 |
