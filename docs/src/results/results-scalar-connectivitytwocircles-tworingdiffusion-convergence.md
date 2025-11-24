# TwoRingDiffusion Convergence

Transient diffusion benchmark on two disconnected circular domains (inner disk
and exterior of an outer disk inside a square box), with a void annulus in
between. Uses the manufactured harmonic solution φ⋆(x,y,t) = exp(-2π²κ t)
sin(πx) sin(πy) restricted to the active regions.

**CSV source:** `results/scalar/ConnectivityTwoCircles/TwoRingDiffusion_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TwoRingDiffusion | 0.03125 | L^2 | 3104 | [3104] | 0.360636 | 0.3482 | 0.0938887 | 0 | NaN | NaN | NaN |
