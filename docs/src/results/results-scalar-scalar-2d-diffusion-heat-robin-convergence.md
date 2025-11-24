# Scalar 2D Diffusion Heat Robin Convergence

Benchmark: Scalar 2D Heat Equation with Robin Interface (legacy Heat.jl port)
This script reproduces the benchmark/Heat.jl workflow but in the BenchPhaseFlow
structure: it performs a mesh-convergence study for the 2D transient diffusion
problem around a heated circle with Robin interface conditions, logs errors, and
writes a single CSV summary (no plots or timestamped folders).

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Robin_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_2D_Diffusion_Heat_Robin | 1 | L^2 | 1 | [2, 2] | 0.192595 | 0.0678528 | 0.180246 | 0 | NaN | NaN | NaN |
| Scalar_2D_Diffusion_Heat_Robin | 0.5 | L^2 | 5 | [4, 4] | 0.112221 | 0.065804 | 0.0909035 | 0 | 0.779222 | 0.044233 | 0.987562 |
| Scalar_2D_Diffusion_Heat_Robin | 0.25 | L^2 | 37 | [8, 8] | 0.0350166 | 0.0282011 | 0.0207572 | 0 | 1.68024 | 1.22242 | 2.13073 |
| Scalar_2D_Diffusion_Heat_Robin | 0.125 | L^2 | 173 | [16, 16] | 0.00786425 | 0.00698765 | 0.00360823 | 0 | 2.15466 | 2.01287 | 2.52425 |
| Scalar_2D_Diffusion_Heat_Robin | 0.0625 | L^2 | 741 | [32, 32] | 0.00226675 | 0.00206698 | 0.00093046 | 0 | 1.79469 | 1.75728 | 1.95527 |
| Scalar_2D_Diffusion_Heat_Robin | 0.03125 | L^2 | 3090 | [64, 64] | 0.000504154 | 0.000463054 | 0.00019938 | 0 | 2.16869 | 2.15827 | 2.22242 |
