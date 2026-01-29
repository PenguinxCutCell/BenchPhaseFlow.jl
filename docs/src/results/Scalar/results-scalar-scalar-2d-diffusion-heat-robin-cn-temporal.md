# Scalar 2D Diffusion Heat Robin CN Temporal

Temporal convergence benchmark for the scalar 2D heat equation with Robin
interface conditions. The geometry, physics, and analytical reference are the
same as `Scalar_2D_Diffusion_Heat_Robin.jl`, but the mesh is fixed while the
time step varies to compare Backward Euler (BE) and Crank–Nicolson (CN).

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Robin_CN_Temporal.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut | order_all | order_full | order_cut | order_all_fit | order_full_fit | order_cut_fit |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_2D_Diffusion_Heat_Robin_CN | 0.16 | L^2 | 50957 | [256, 256] | 0.251798 | 0.251381 | 0.0144785 | 0 | NaN | NaN | NaN | 1.5 | 1.5 | 1.2 | 1 | 1 | 0.9 |
| Scalar_2D_Diffusion_Heat_Robin_CN | 0.08 | L^2 | 50957 | [256, 256] | 0.169847 | 0.168718 | 0.0195496 | 0 | 0.568033 | 0.575264 | -0.433226 | 1.5 | 1.5 | 1.2 | 1 | 1 | 0.9 |
| Scalar_2D_Diffusion_Heat_Robin_CN | 0.04 | L^2 | 50957 | [256, 256] | 0.0781337 | 0.0779758 | 0.00496515 | 0 | 1.12021 | 1.11351 | 1.97723 | 1.5 | 1.5 | 1.2 | 1 | 1 | 0.9 |
| Scalar_2D_Diffusion_Heat_Robin_CN | 0.02 | L^2 | 50957 | [256, 256] | 0.0272105 | 0.0271749 | 0.00139113 | 0 | 1.52178 | 1.52075 | 1.83559 | 1.5 | 1.5 | 1.2 | 1 | 1 | 0.9 |
| Scalar_2D_Diffusion_Heat_Robin_CN | 0.01 | L^2 | 50957 | [256, 256] | 0.0270736 | 0.0269875 | 0.00215701 | 0 | 0.00727781 | 0.0099834 | -0.632782 | 1.5 | 1.5 | 1.2 | 1 | 1 | 0.9 |
| Scalar_2D_Diffusion_Heat_Robin_CN | 0.005 | L^2 | 50957 | [256, 256] | 0.00689741 | 0.00683556 | 0.000921606 | 0 | 1.97276 | 1.98116 | 1.22681 | 1.5 | 1.5 | 1.2 | 1 | 1 | 0.9 |
| Scalar_2D_Diffusion_Heat_Robin_CN | 0.0025 | L^2 | 50957 | [256, 256] | 0.00345711 | 0.00343209 | 0.000415167 | 0 | 0.996486 | 0.99397 | 1.15046 | 1.5 | 1.5 | 1.2 | 1 | 1 | 0.9 |
