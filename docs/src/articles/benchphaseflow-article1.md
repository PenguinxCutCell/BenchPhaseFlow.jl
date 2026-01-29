# A Space-Time Cartesian Cut-Cell Method for Two-Phase Diffusion Problems using a Two-Fluid Approach

Refer to the original publication for details:
[https://hal.science/hal-05457114](https://hal.science/hal-05457114)


## Benchmarks

Static problems :
Monophasic : 
- Johansen & Colella Problem 1 : 2D Star shaped poisson problem with Dirichlet BCs
- Johansen & Colella Problem 2 : 2D Flower-shaped laplace problem
- 2D Poisson equation in a disk with Robin BCs
- 3D Heat equation in a sphere with Robin BCs

Diphasic :
- Diphasic 2D unsteady diffusion with continuity/jump conditions at the interface
- Diphasic 2D unsteady diffusion with continuity/jump conditions at the interface with various Henry/diffusivity ratios
- Diphasic 3D unsteady diffusion with continuity conditions at the interface

Moving problems :
Monophasic :
- 2D Diffusion with an oscillating circular boundary
- 2D Diffusion with multiples ellipses (McCorquodale-Colella problem)
- 3D Schwartz-Colella diffusion on a expanding sphere

Diphasic :
- Manufactured solution with an oscillating two-phase interface

## Revision History

- 20/01/2026 : Paper send for review. Data and benchmarks using Penguin.jl v20/01/2026.
- 29/01/2026 : Updated benchmarks using Penguin.jl v29/01/2026 : Corrections in small cut-cell handling/time integration => Less oscillations in convergence results. Errors remains sensibly unchanged. Overall fit remain sensibly unchanged.