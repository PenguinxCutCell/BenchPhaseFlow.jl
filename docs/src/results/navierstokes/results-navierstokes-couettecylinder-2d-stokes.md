# CouetteCylinder 2D Stokes

Steady Stokes Couette flow between two concentric cylinders (CSV only).

Inner cylinder (R₁=0.25) rotates with ω=1, outer cylinder (R₂=0.5) is fixed.
Level-set defines an annulus domain. Solves with `StokesMono` and a
component-aware cut BC, then samples tangential velocity along x=0 to compute
L2/Linf errors against the analytic Couette profile. Writes summary CSV to
`results/NavierStokes/CouetteCylinder_2D_Stokes.csv`.

**CSV source:** `results/NavierStokes/CouetteCylinder_2D_Stokes.csv`

| l2_err | linf_err | radial_drift | n_profile | nx | ny | R_inner | R_outer |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.119062 | 0.215888 | 0 | 34 | 256 | 256 | 0.25 | 0.5 |
