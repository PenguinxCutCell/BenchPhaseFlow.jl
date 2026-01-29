# Scalar 2D Diffusion Heat Dirichlet Convergence

Scalar 2D diffusion (radial) heat equation convergence test with Dirichlet boundary conditions. The benchmark uses a disk of radius $R$ centered at $c$ and compares against the radial series solution

```math
u(r,t) = \Big(1 - \sum_{n=1}^{\infty} A_n\, e^{-a\,\alpha_n^2 t / R^2}\, J_0\!\big(\alpha_n\, r / R\big)\Big)\,(w_r - w_0) + w_0,
```

where $r=\lVert x-c\rVert$, $J_0$ is the Bessel function of order zero, $\alpha_n$ are roots of the radial eigenvalue equation, and $A_n$ are the corresponding coefficients. Dirichlet boundary condition $w_r$ is applied at the boundary. Errors and observed orders are reported for a mesh refinement sweep.

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Dirichlet_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Dirichlet_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```