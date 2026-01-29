# Stokes 3D Sphere

3D steady Stokes flow around a spherical obstacle (no plots).

The flow is driven by a uniform free-stream velocity in the x-direction and a
no-slip sphere embedded via a level-set capacity. Boundary velocities are
imposed from the analytical creeping-flow solution around a sphere so the
numerical solution should converge to the exact field. Results are written to
`results/NavierStokes/Stokes_3D_Sphere.csv`.

**CSV source:** `results/NavierStokes/Stokes_3D_Sphere.csv`

| N | h | err_ux | err_uy | err_uz | rate_ux | rate_uy | rate_uz |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 8 | 1 | 0.0719353 | 0.0280612 | 0.0281089 | NaN | NaN | NaN |
| 12 | 0.666667 | 0.0450305 | 0.0180912 | 0.018077 | 1.15528 | 1.08261 | 1.08875 |
| 16 | 0.5 | 0.0326388 | 0.0134212 | 0.013423 | 1.11873 | 1.03792 | 1.03472 |
| 24 | 0.333333 | 0.021304 | 0.00862704 | 0.00862542 | 1.05214 | 1.08994 | 1.09073 |
| 32 | 0.25 | 0.0157672 | 0.00638838 | 0.00639054 | 1.04617 | 1.04428 | 1.04245 |
