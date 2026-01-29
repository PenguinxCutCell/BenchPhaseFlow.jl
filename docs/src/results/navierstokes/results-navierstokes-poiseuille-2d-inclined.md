# Poiseuille 2D Inclined

2D Navier–Stokes Poiseuille flow in inclined cut-cell channels (no plotting).

Mirrors `examples/2D/Stokes/poiseuille_2d_cut_inclined.jl` but uses the steady
Navier–Stokes mono solver. For multiple inclination angles, solves the flow,
computes volume-weighted L2 errors split by cut/full cells, and writes results
to `results/NavierStokes/Poiseuille_2D_Inclined.csv`.

**CSV source:** `results/NavierStokes/Poiseuille_2D_Inclined.csv`

| l2_total_cut | l2_total_all | iterations | n_full_cells | n_cut_cells | l2_total_full | residual | ratio | angle_deg |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.00106965 | 0.0584065 | 4 | 2331 | 252 | 0.0583967 | 1.73321e-08 | 0.018317 | 0 |
| 0.0146747 | 0.0482983 | 20 | 2290 | 274 | 0.046015 | 2.69578e-07 | 0.318911 | 5 |
| 0.014488 | 0.0445405 | 7 | 2312 | 296 | 0.0421183 | 4.04591e-07 | 0.343984 | 10 |
| 0.0138275 | 0.039165 | 7 | 2346 | 320 | 0.0366428 | 2.28705e-07 | 0.37736 | 15 |
| 0.0142688 | 0.0368522 | 100 | 2399 | 343 | 0.0339778 | 46.0217 | 0.419945 | 20 |
| 0.0145377 | 0.0377826 | 6 | 2487 | 369 | 0.0348737 | 4.20095e-07 | 0.416868 | 25 |
| 0.016784 | 0.0423109 | 6 | 2594 | 397 | 0.0388395 | 6.08314e-07 | 0.432137 | 30 |
| 0.0185595 | 0.0496719 | 100 | 2734 | 407 | 0.0460743 | 141.619 | 0.402816 | 35 |
| 0.0207379 | 0.0587875 | 100 | 2823 | 397 | 0.0550082 | 421.974 | 0.376997 | 40 |
| 0.0151164 | 0.0684598 | 100 | 2938 | 392 | 0.06677 | 35.4309 | 0.226395 | 45 |
