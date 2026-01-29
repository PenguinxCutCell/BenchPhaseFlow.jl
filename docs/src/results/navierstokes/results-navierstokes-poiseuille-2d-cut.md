# Poiseuille 2D Cut

2D Navier–Stokes Poiseuille benchmark (steady) with cut-cell channel (no plots).

Adapted from `examples/2D/Stokes/poiseuille_2d_cut.jl` but uses the steady
Navier–Stokes mono solver. Computes mid-column profile errors against the
analytical parabola and writes results to `results/NavierStokes`.

**CSV source:** `results/NavierStokes/Poiseuille_2D_Cut.csv`

| iterations | residual | l2_profile | linf_profile | rel_center_err | mean_ux | max_ux | mean_uy | max_uy | mean_p | min_p | max_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 5 | 6.44891e-11 | 0.00365297 | 0.00134008 | 0.00131723 | 0.39392 | 0.998683 | -1.39913e-11 | 0.0242636 | 0.676726 | -40.0996 | 26.0195 |
