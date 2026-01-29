# Scalar 2D Diffusion Heat Robin BE Temporal

Temporal convergence benchmark for the scalar 2D heat equation with Robin
interface conditions. The geometry, physics, and analytical reference are the
same as `Scalar_2D_Diffusion_Heat_Robin.jl`, but the mesh is fixed while the
time step varies to compare Backward Euler (BE) and Crank–Nicolson (CN).

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Robin_BE_Temporal.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut | order_all | order_full | order_cut | order_all_fit | order_full_fit | order_cut_fit |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_2D_Diffusion_Heat_Robin_BE | 0.16 | L^2 | 50957 | [256, 256] | 0.215184 | 0.214419 | 0.0181343 | 0 | NaN | NaN | NaN | 1.5 | 1.5 | 1.6 | 1 | 1 | 1 |
| Scalar_2D_Diffusion_Heat_Robin_BE | 0.08 | L^2 | 50957 | [256, 256] | 0.15327 | 0.152697 | 0.0132381 | 0 | 0.489496 | 0.489756 | 0.454032 | 1.5 | 1.5 | 1.6 | 1 | 1 | 1 |
| Scalar_2D_Diffusion_Heat_Robin_BE | 0.04 | L^2 | 50957 | [256, 256] | 0.0709314 | 0.0706601 | 0.00619751 | 0 | 1.11158 | 1.11171 | 1.09493 | 1.5 | 1.5 | 1.6 | 1 | 1 | 1 |
| Scalar_2D_Diffusion_Heat_Robin_BE | 0.02 | L^2 | 50957 | [256, 256] | 0.0239428 | 0.0238614 | 0.00197237 | 0 | 1.56683 | 1.56621 | 1.65176 | 1.5 | 1.5 | 1.6 | 1 | 1 | 1 |
| Scalar_2D_Diffusion_Heat_Robin_BE | 0.01 | L^2 | 50957 | [256, 256] | 0.0253672 | 0.0252598 | 0.00233207 | 0 | -0.0833702 | -0.0821597 | -0.241681 | 1.5 | 1.5 | 1.6 | 1 | 1 | 1 |
| Scalar_2D_Diffusion_Heat_Robin_BE | 0.005 | L^2 | 50957 | [256, 256] | 0.00611887 | 0.00609777 | 0.000507652 | 0 | 2.05163 | 2.05049 | 2.1997 | 1.5 | 1.5 | 1.6 | 1 | 1 | 1 |
| Scalar_2D_Diffusion_Heat_Robin_BE | 0.0025 | L^2 | 50957 | [256, 256] | 0.00307128 | 0.00306041 | 0.000258179 | 0 | 0.994423 | 0.994556 | 0.975471 | 1.5 | 1.5 | 1.6 | 1 | 1 | 1 |
