# Heat 2ph 1D Convergence

Diphasic 1D heat diffusion benchmark reproduced from `benchmark/Heat_2ph_1D.jl`.
Performs a mesh-convergence study for the unsteady diphasic heat equation with
a planar interface and writes convergence data (no plotting).

**CSV source:** `results/scalar/diphasic/Heat_2ph_1D_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_phase1 | inside_cells_phase2 | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | phase1_all_err | phase1_full_err | phase1_cut_err | phase1_empty_err | phase2_all_err | phase2_full_err | phase2_cut_err | phase2_empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Heat_2ph_1D | 2 | L^2 | 3 | 1 | 2 | [1, 2] | 0.197475 | 0.114341 | 0.161005 | 0 | 0.149154 | 3.16156e-06 | 0.149154 | 0 | 0.197475 | 0.114341 | 0.161005 | 0 | NaN | NaN | NaN |
| Heat_2ph_1D | 1 | L^2 | 7 | 3 | 4 | [3, 4] | 0.113717 | 0.0931905 | 0.0651702 | 0 | 0.113717 | 0.0931905 | 0.0651702 | 0 | 0.100819 | 0.0827844 | 0.0575427 | 0 | 0.796218 | 0.295088 | 1.30482 |
| Heat_2ph_1D | 0.5 | L^2 | 15 | 7 | 8 | [7, 8] | 0.028828 | 0.0223082 | 0.0182592 | 0 | 0.028828 | 0.0223082 | 0.0182592 | 0 | 0.0270792 | 0.0209549 | 0.0171515 | 0 | 1.97991 | 2.06261 | 1.83559 |
| Heat_2ph_1D | 0.25 | L^2 | 31 | 15 | 16 | [15, 16] | 0.0138454 | 0.0138334 | 0.000576233 | 0 | 0.0138454 | 0.0138334 | 0.000576233 | 0 | 0.0134193 | 0.0134077 | 0.000558498 | 0 | 1.05806 | 0.689416 | 4.98583 |
| Heat_2ph_1D | 0.125 | L^2 | 63 | 31 | 32 | [31, 32] | 0.00174211 | 0.00174104 | 6.11541e-05 | 0 | 0.00174211 | 0.00174104 | 6.11541e-05 | 0 | 0.0017151 | 0.00171404 | 6.02059e-05 | 0 | 2.9905 | 2.99014 | 3.23613 |
| Heat_2ph_1D | 0.0625 | L^2 | 127 | 63 | 64 | [63, 64] | 0.000947431 | 0.000947351 | 1.23102e-05 | 0 | 0.000947431 | 0.000947351 | 1.23102e-05 | 0 | 0.000940058 | 0.000939979 | 1.22144e-05 | 0 | 0.878745 | 0.877977 | 2.31259 |
| Heat_2ph_1D | 0.03125 | L^2 | 255 | 127 | 128 | [127, 128] | 0.000109152 | 0.000109151 | 4.58115e-07 | 0 | 0.000109152 | 0.000109151 | 4.58115e-07 | 0 | 0.000108726 | 0.000108725 | 4.56329e-07 | 0 | 3.11768 | 3.11757 | 4.748 |
