# Scalar 1D Diffusion Heat Robin Linf Convergence

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

**CSV source:** `results/scalar/Scalar_1D_Diffusion_Heat_Robin_Linf_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_1D_Diffusion_Heat_Robin_Linf | 2.5 | L^Inf | 4 | [1] | 0.597631 | 0.597631 | 0 | 0.111696 | NaN | NaN | NaN |
| Scalar_1D_Diffusion_Heat_Robin_Linf | 1.25 | L^Inf | 8 | [2] | 0.271818 | 0.271818 | 0 | 0.111696 | 1.13661 | 1.13661 | NaN |
| Scalar_1D_Diffusion_Heat_Robin_Linf | 0.625 | L^Inf | 15 | [4] | 0.0868778 | 0.0623712 | 0.0868778 | 0.111696 | 1.64558 | 2.12369 | NaN |
| Scalar_1D_Diffusion_Heat_Robin_Linf | 0.3125 | L^Inf | 30 | [7] | 0.0155127 | 0.0155127 | 0.00760484 | 0.186078 | 2.48553 | 2.00742 | 3.514 |
| Scalar_1D_Diffusion_Heat_Robin_Linf | 0.15625 | L^Inf | 61 | [13] | 0.00282169 | 0.00282169 | 0.00128422 | 0.186078 | 2.45882 | 2.45882 | 2.56603 |
| Scalar_1D_Diffusion_Heat_Robin_Linf | 0.078125 | L^Inf | 121 | [26] | 0.00124695 | 0.00124695 | 0.000391983 | 0.222726 | 1.17816 | 1.17816 | 1.71202 |
