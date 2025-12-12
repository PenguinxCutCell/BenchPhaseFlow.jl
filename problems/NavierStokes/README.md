# Navier-Stokes Problems

## Equations 

Stokes Equations:
$$
\partial_t \mathbf{u} - \nu \nabla^2 \mathbf{u} + \nabla p = \mathbf{f} \\
\nabla \cdot \mathbf{u} = 0
$$

Navier-Stokes Equations:
$$
\partial_t \mathbf{u} + (\mathbf{u} \cdot \nabla) \mathbf{u} - \nu \nabla^2 \mathbf{u} + \nabla p = \mathbf{f} \\
\nabla \cdot \mathbf{u} = 0
$$
where $\mathbf{u}$ is the velocity field, $p$ is the pressure field, $\nu$ is the kinematic viscosity, and $\mathbf{f}$ represents external forces.

Accompagnied by various boundary conditions such as Dirichlet, Neumann, Traction, Outflow ...


## List of Problems
- 1D Stokes Flow with Dirichlet BCs : `Stokes_1D_Mono.jl`
- 2D Taylor-Green Vortex : `Taylor_Green_2D_Mono.jl`
- 2D Poiseuille Flow : `Poiseuille_2D_Cut.jl`, `Poiseuille_2D_Inclined.jl`
- Oscillatory Channel Flow : `OscillatoryChannel_2D_Mono.jl`
- Lid-Driven Cavity Flow : `LidDrivenCavity_2D_Ghia.jl`
- Kovasznay Flow : `Kovasznay_2D_Mono.jl`
- 3D Stokes Flow around a Sphere : `Stokes_3D_Sphere.jl`
- 2D Flow Past a Cylinder : `CylinderWake_Steady.jl`, `CylinderWake_ReSweep.jl`, `CylinderRe100_2D_Unsteady.jl`
- Couette Cylindrical Flow : `CouetteCylinder_2D_Stokes.jl`
- Backward Facing Step Flow : `BackwardFacingStep_Steady.jl`
