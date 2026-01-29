# JohansenColella P2 Convergence

Johansen–Colella Problem 2: variable-coefficient Poisson equation on the same
star-shaped domain as Problem 1. The diffusion coefficient is β(r)=1−r² and the
equation ∇·(β ∇ϕ) = (7 r² − 15 r⁴) cos(3θ) retains the exact solution
ϕ(r, θ) = r⁴ cos(3θ). Interface Dirichlet data enforce the exact boundary values.

**CSV source:** `results/scalar/JohansenColella_P2_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut | trunc_max_all | trunc_max_full | trunc_max_cut | trunc_pair_order_max_all | trunc_pair_order_max_full | trunc_pair_order_max_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| JohansenColella_P2 | 0.025 | L^1 | 401 | [36, 36] | 3.41573e-05 | 1.77916e-05 | 1.63657e-05 | 0 | NaN | NaN | NaN | 0.252556 | 0.00179055 | 0.252556 | NaN | NaN | NaN |
| JohansenColella_P2 | 0.0125 | L^1 | 1813 | [72, 72] | 1.05202e-05 | 6.35284e-06 | 4.16735e-06 | 0 | 1.69903 | 1.48572 | 1.97348 | 0.672989 | 0.0166378 | 0.672989 | -1.41398 | -3.21599 | -1.41398 |
| JohansenColella_P2 | 0.00625 | L^1 | 7713 | [144, 144] | 2.79038e-06 | 1.76719e-06 | 1.02318e-06 | 0 | 1.91463 | 1.84594 | 2.02607 | 0.536003 | 0.0197929 | 0.536003 | 0.328342 | -0.250517 | 0.328342 |
| JohansenColella_P2 | 0.003125 | L^1 | 31717 | [288, 288] | 7.51003e-07 | 5.01963e-07 | 2.4904e-07 | 0 | 1.89357 | 1.81581 | 2.03862 | 0.270085 | 0.0211271 | 0.270085 | 0.988826 | -0.0941129 | 0.988826 |
| JohansenColella_P2 | 0.0015625 | L^1 | 128581 | [576, 576] | 2.08586e-07 | 1.45631e-07 | 6.29552e-08 | 0 | 1.84818 | 1.78527 | 1.98398 | 1.34173 | 0.0192736 | 1.34173 | -2.31261 | 0.132472 | -2.31261 |
| JohansenColella_P2 | 0.025 | L^Inf | 401 | [36, 36] | 0.000454436 | 0.000160453 | 0.000454436 | 0.176777 | NaN | NaN | NaN | 0.252556 | 0.00179055 | 0.252556 | NaN | NaN | NaN |
| JohansenColella_P2 | 0.0125 | L^Inf | 1813 | [72, 72] | 0.000304023 | 8.29117e-05 | 0.000304023 | 0.176777 | 0.579899 | 0.952507 | 0.579899 | 0.672989 | 0.0166378 | 0.672989 | -1.41398 | -3.21599 | -1.41398 |
| JohansenColella_P2 | 0.00625 | L^Inf | 7713 | [144, 144] | 0.000225975 | 2.44929e-05 | 0.000225975 | 0.176777 | 0.428015 | 1.75921 | 0.428015 | 0.536003 | 0.0197929 | 0.536003 | 0.328342 | -0.250517 | 0.328342 |
| JohansenColella_P2 | 0.003125 | L^Inf | 31717 | [288, 288] | 0.000108896 | 1.00478e-05 | 0.000108896 | 0.176777 | 1.05322 | 1.28548 | 1.05322 | 0.270085 | 0.0211271 | 0.270085 | 0.988826 | -0.0941129 | 0.988826 |
| JohansenColella_P2 | 0.0015625 | L^Inf | 128581 | [576, 576] | 5.62735e-05 | 3.86249e-06 | 5.62735e-05 | 0.176777 | 0.952419 | 1.37928 | 0.952419 | 1.34173 | 0.0192736 | 1.34173 | -2.31261 | 0.132472 | -2.31261 |
