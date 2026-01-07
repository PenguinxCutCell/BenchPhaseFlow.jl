# Diphasic Diffusion Problems

This directory contains implementations of diphasic diffusion problems.

## Mathematical Formulation
Given $\Omega \subset \mathbb{R}^d$ with $d=1,2,3$, divided into two subdomains $\Omega_1$ and $\Omega_2$ by an interface $\Gamma$, we solve the diphasic transient diffusion equation:
$$
\frac{\partial u}{\partial t} = \nabla \cdot (D_i \nabla u) + f_i, \quad \text{in } \Omega_i, \quad i=1,2,
$$

with interface conditions:
$$
[[\lambda u]]_\Gamma = g_s, \quad [[D \nabla u \cdot \mathbf{n}]]_\Gamma = g_f,
$$
where $D_i$ are the diffusion coefficients, $f_i$ are source terms, $g_s$ and $g_f$ are the scalar and flux jumps across the interface, respectively.

## Implemented Problems
- **LiuFedkiw/Case1.jl**: Implements the first case from Liu and Fedkiw's benchmark with specified scalar and flux jumps.
- **LiuFedkiw/Case2.jl**: Implements the second case from Liu and Fedkiw's benchmark with different interface conditions.
- **LiuFedkiw/Case3_Poisson2D_Circle.jl**: Steady 2D Poisson benchmark with interior exponential solution, nonzero source/jumps, and a circular interface.
- **Manufactured_Continuous2D.jl**: Steady 2D manufactured Poisson test with continuous value and flux across a circular interface.
- **Heat_2ph_3D.jl**: Transient 3D bubble benchmark with continuous value/flux, using an integral analytical solution for convergence.
- **Poisson2D_SineInterface.jl**: Steady 2D Poisson test with a sine-wave interface, manufactured solution `sin(πx)sin(πy)`, and unequal diffusivities.
- **UnsteadyStarManufactured.jl**: Transient 2D diphasic diffusion with a star-shaped interface and manufactured exponential solution.
- **Heat_2ph_1D.jl**: A 1D diphasic heat diffusion problem with homothetic interface jump conditions.
- **Heat_2ph_2D.jl**: A 2D diphasic heat diffusion problem with a circular interface and Henry-law jump conditions.
