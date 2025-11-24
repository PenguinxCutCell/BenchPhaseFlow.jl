# Scalar 1D Diffusion Heat Robin Convergence

Scalar 1D Diffusion Heat Equation with Robin Boundary Conditions
This benchmark mirrors `examples/1D/Diffusion/Heat_robin.jl`: a 1D transient
diffusion problem on `x ∈ [0, 10]` with interface position `center = 0.25`,
diffusivity `a = 5`, and homogeneous Robin condition `uₓ + u = 0` at the cut.
The exact solution used for convergence is the complementary-error-function
profile
```
u(x,t) = erf(η) + exp(k(x-center) + a k² t) ⋅ erfc(η + k √(a t)),
η = (x - center)/(2 √(a t)), k = 1,
```
which satisfies the same boundary data (Dirichlet `u=1` on the left, `u=0`
on the right) and source-free diffusion equation `u_t = a u_{xx}`.

**CSV source:** `results/scalar/Scalar_1D_Diffusion_Heat_Robin_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_1D_Diffusion_Heat_Robin | 2.5 | L^2 | 4 | [1] | 0.518943 | 0.518943 | 0 | 0 | NaN | NaN | NaN |
| Scalar_1D_Diffusion_Heat_Robin | 1.25 | L^2 | 8 | [2] | 0.192859 | 0.192859 | 0 | 0 | 1.42803 | 1.42803 | NaN |
| Scalar_1D_Diffusion_Heat_Robin | 0.625 | L^2 | 15 | [4] | 0.0288905 | 0.0224832 | 0.018143 | 0 | 2.73888 | 3.10063 | NaN |
| Scalar_1D_Diffusion_Heat_Robin | 0.3125 | L^2 | 30 | [7] | 0.0100259 | 0.0099445 | 0.00127525 | 0 | 1.52686 | 1.17688 | 3.83056 |
| Scalar_1D_Diffusion_Heat_Robin | 0.15625 | L^2 | 61 | [13] | 0.00184667 | 0.00184495 | 7.97242e-05 | 0 | 2.44074 | 2.43032 | 3.99962 |
| Scalar_1D_Diffusion_Heat_Robin | 0.078125 | L^2 | 121 | [26] | 0.000637008 | 0.000636046 | 3.49921e-05 | 0 | 1.53554 | 1.53637 | 1.18799 |
