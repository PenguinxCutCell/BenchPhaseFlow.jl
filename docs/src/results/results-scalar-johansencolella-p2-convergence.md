# JohansenColella P2 Convergence

Johansen–Colella Problem 2: variable-coefficient Poisson equation on the same
star-shaped domain as Problem 1. The diffusion coefficient is β(r)=1−r² and the
equation ∇·(β ∇ϕ) = (7 r² − 15 r⁴) cos(3θ) retains the exact solution
ϕ(r, θ) = r⁴ cos(3θ). Interface Dirichlet data enforce the exact boundary values.
NOTE: THIS SCRIPT DOES NOT PRODUCE ACCURATE SOLUTIONS DUE TO A BUG IN THE
VARIABLE-COEFFICIENT DIFFUSION OPERATOR IMPLEMENTATION.

**CSV source:** `results/scalar/JohansenColella_P2_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| JohansenColella_P2 | 0.0625 | L^2 | 35 | [15, 15] | 0.00345643 | 0.00129834 | 0.00320331 | 0 |
| JohansenColella_P2 | 0.03125 | L^2 | 235 | [29, 29] | 0.000153343 | 0.000124648 | 8.93132e-05 | 0 |
| JohansenColella_P2 | 0.015625 | L^2 | 1143 | [58, 58] | 0.000115947 | 0.000111447 | 3.19884e-05 | 0 |
| JohansenColella_P2 | 0.0078125 | L^2 | 4867 | [116, 116] | 0.000108222 | 0.000107622 | 1.13723e-05 | 0 |
