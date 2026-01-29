# Scalar 2D Diffusion Heat Robin Convergence

Temporal and spatial convergence benchmark for the scalar 2D heat equation with Robin boundary conditions. The geometry, physics, and analytical reference are the same as `Scalar_2D_Diffusion_Heat_Robin.jl`, with errors and observed orders reported for a mesh refinement sweep.

Analytical solution:

```math
u(r,t) = \Big(1 - \sum_{n=1}^{\infty} A_n\, e^{-a\,\alpha_n^2 t / R^2}\, J_0\!\big(\alpha_n\, r / R\big)\Big)\,(w_r - w_0) + w_0,
```

where $r=\lVert x-c\rVert$, $J_0$ is the Bessel function of order zero, $\alpha_n$ are roots of the radial eigenvalue equation, and $A_n$ are the corresponding coefficients. Robin boundary condition $w_r$ is applied at the boundary.



**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Robin_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Robin_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Robin_Linf_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Robin_Linf_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```