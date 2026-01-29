# Scalar 2D Diffusion Heat Robin Convergence

Benchmark: Scalar 2D Heat Equation with Robin Interface (legacy Heat.jl port)
This script reproduces the benchmark/Heat.jl workflow but in the BenchPhaseFlow
structure: it performs a mesh-convergence study for the 2D transient diffusion
problem around a heated circle with Robin interface conditions, logs errors, and
writes a single CSV summary (no plots or timestamped folders).

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Robin_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_2D_Diffusion_Heat_Robin | 1 | L^2 | 1 | [2, 2] | 0.152997 | 0.0435451 | 0.146669 | 0 | NaN | NaN | NaN |
| Scalar_2D_Diffusion_Heat_Robin | 0.5 | L^2 | 5 | [4, 4] | 0.0830571 | 0.0477483 | 0.0679602 | 0 | 0.881327 | -0.132939 | 1.10981 |
| Scalar_2D_Diffusion_Heat_Robin | 0.25 | L^2 | 37 | [8, 8] | 0.0229443 | 0.0177509 | 0.0145378 | 0 | 1.85597 | 1.42756 | 2.22488 |
| Scalar_2D_Diffusion_Heat_Robin | 0.125 | L^2 | 173 | [16, 16] | 0.00576651 | 0.00502576 | 0.00282743 | 0 | 1.99237 | 1.82048 | 2.36225 |
| Scalar_2D_Diffusion_Heat_Robin | 0.0625 | L^2 | 741 | [32, 32] | 0.00148637 | 0.00129613 | 0.000727564 | 0 | 1.9559 | 1.95513 | 1.95835 |
| Scalar_2D_Diffusion_Heat_Robin | 0.03125 | L^2 | 3090 | [64, 64] | 0.000376361 | 0.000331417 | 0.000178355 | 0 | 1.9816 | 1.96749 | 2.02832 |
