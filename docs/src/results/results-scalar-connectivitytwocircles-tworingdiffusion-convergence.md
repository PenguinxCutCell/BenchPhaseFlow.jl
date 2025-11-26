# TwoRingDiffusion Convergence

Transient diffusion benchmark on two disconnected circular domains (inner disk
and exterior of an outer disk inside a square box), with a void annulus in
between. Uses the manufactured harmonic solution φ⋆(x,y,t) = exp(-2π²κ t)
sin(πx) sin(πy) restricted to the active regions.

**CSV source:** `results/scalar/ConnectivityTwoCircles/TwoRingDiffusion_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TwoRingDiffusion | 0.125 | L^2 | 168 | [168] | 0.0579924 | 0.0496782 | 0.0299199 | 0 | NaN | NaN | NaN |
| TwoRingDiffusion | 0.0625 | L^2 | 760 | [760] | 0.0564319 | 0.0528128 | 0.019884 | 0 | 0.0393523 | -0.0882754 | 0.589499 |
| TwoRingDiffusion | 0.03125 | L^2 | 3128 | [3128] | 0.0561974 | 0.0539152 | 0.0158526 | 0 | 0.00600748 | -0.0298042 | 0.326887 |
