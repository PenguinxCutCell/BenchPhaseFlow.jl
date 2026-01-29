# Diph Poisson2D SineInterface Convergence

Diphasic 2D Poisson convergence test with a sine-wave interface.

Domain: [0, 1] x [0, 1]. Interface y = 0.5 + 0.1 * sin(2πx).
Phase 1: y <= interface; Phase 2: y >= interface.
Diffusivities: D1 = 1, D2 = 2.

Manufactured solution (zero on the box boundary):
    φ(x, y) = sin(πx) * sin(πy)
    u1 = a * φ, u2 = b * φ with b = a * D1 / D2
This enforces flux continuity (β2 ∂n u2 = β1 ∂n u1) while allowing a scalar jump
[u] = (b - a) φ. Set a = 1 here.

Laplacian: Δφ = -2π^2 φ, so forcing per phase:
    f1 = 2π^2 * D1 * a * φ
    f2 = 2π^2 * D2 * b * φ = 2π^2 * D1 * a * φ (same value).

Interface conditions:
    ScalarJump = (b - a) * φ
    FluxJump   = 0
Boundary: Dirichlet exact solution (φ = 0 on all borders).

**CSV source:** `results/scalar/diphasic/Diph_Poisson2D_SineInterface_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_phase1 | inside_cells_phase2 | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | phase1_all_err | phase1_full_err | phase1_cut_err | phase1_empty_err | phase2_all_err | phase2_full_err | phase2_cut_err | phase2_empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Diph_Poisson2D_SineInterface | 0.0625 | L^2 | 234 | 123 | 111 | [123, 111] | 0.120096 | 0.114655 | 0.0357395 | 0 | 0.0575143 | 0.0512952 | 0.0260135 | 0 | 0.120096 | 0.114655 | 0.0357395 | 0 | NaN | NaN | NaN |
| Diph_Poisson2D_SineInterface | 0.03125 | L^2 | 978 | 497 | 481 | [497, 481] | 0.0582162 | 0.0571768 | 0.0109515 | 0 | 0.0286384 | 0.0270319 | 0.00945699 | 0 | 0.0582162 | 0.0571768 | 0.0109515 | 0 | 1.0447 | 1.0038 | 1.70639 |
| Diph_Poisson2D_SineInterface | 0.015625 | L^2 | 4006 | 1999 | 2007 | [1999, 2007] | 0.02868 | 0.0284232 | 0.00382925 | 0 | 0.0143553 | 0.0139522 | 0.00337782 | 0 | 0.02868 | 0.0284232 | 0.00382925 | 0 | 1.02138 | 1.00836 | 1.51599 |
| Diph_Poisson2D_SineInterface | 0.0078125 | L^2 | 16204 | 8002 | 8202 | [8002, 8202] | 0.0141978 | 0.0141387 | 0.00129485 | 0 | 0.00715411 | 0.00704987 | 0.0012168 | 0 | 0.0141978 | 0.0141387 | 0.00129485 | 0 | 1.01438 | 1.00743 | 1.56427 |
