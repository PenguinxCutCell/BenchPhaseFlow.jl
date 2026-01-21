# Gibou-Fedkiw Problems

This directory contains implementations of the Gibou-Fedkiw benchmark problems
for diffusion equations with interfaces, as described in:
- Gibou, F. and Fedkiw, R. (2005). A second-order accurate symmetric discretization
  of the Poisson equation on irregular domains. Journal of Computational Physics,
  200(2), 492-520. https://doi.org/10.1016/j.jcp.2004.08.015

## Problem 1: Poisson Equation in 1D with interface
Implemented in `Poisson1D.jl`, this problem involves solving the Poisson equation
in a one-dimensional domain with an interface. The exact solution is known, allowing for
error analysis and convergence studies.

## Problem 2: Poisson Equation in 2D with star-shaped interface
Implemented in `Poisson2D.jl`, this problem involves solving the Poisson equation
in a two-dimensional domain with a star-shaped interface. The exact solution is known,
allowing for error analysis and convergence studies.

## Problem 3: Poisson Equation in 3D with star-shaped interface
Implemented in `Poisson3D.jl`, this problem involves solving the Poisson equation
in a three-dimensional domain with a star-shaped interface. The exact solution is known,
allowing for error analysis and convergence studies.

## Problem 4: Heat Equation in 1D with interface
Implemented in `Heat1D.jl`, this problem involves solving the heat equation
in a one-dimensional domain with an interface. The exact solution is known,
allowing for error analysis and convergence studies.

## Problem 5: Heat Equation in 2D with star-shaped interface
Implemented in `Heat2D.jl`, this problem involves solving the heat equation
in a two-dimensional domain with a star-shaped interface. The exact solution is known,
allowing for error analysis and convergence studies.

## Problem 6: Heat Equation in 3D with star-shaped interface
Implemented in `Heat3D.jl`, this problem involves solving the heat equation
in a three-dimensional domain with a star-shaped interface. The exact solution is known,
allowing for error analysis and convergence studies.


# Results

## Problem 4

| Case | N / Points | Error | Order |
|------|------------|-------|-------|
| This work | 41 | 1.7988e-2 | – |
|          | 81 | 9.4519e-3 | 0.93 |
|          | 161 | 4.8470e-3 | 0.97 |
| Gibou et al. | 41  | 1.443e-2 | – |
|              | 81  | 7.240e-3 | 0.99 |
|              | 161 | 3.634e-3 | 0.99 |

| Case | N / Points | Error | Order |
|------|------------|-------|-------|
| This work | 41 | 4.5533e-4 | – |
|          | 81 | 1.1607e-4 | 1.97 |
|          | 161 | 2.8908e-5 | 2.01 |
| Gibou et al. | 41  | 4.084e-4 | – |
|              | 81  | 9.907e-5 | 2.01 |
|              | 161 | 2.424e-5 | 2.03 |
