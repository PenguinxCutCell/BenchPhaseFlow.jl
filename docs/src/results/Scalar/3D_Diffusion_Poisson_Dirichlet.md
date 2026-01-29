# Scalar 3D Diffusion Poisson Dirichlet Convergence

Scalar 3D Diffusion Poisson Problem with Dirichlet Boundary Conditions
This script performs a mesh-convergence study for the 3D Poisson equation
using the Penguin.jl library. It computes errors and estimated orders of convergence
for different mesh sizes and writes the results to a CSV file.

The analytical solution is given by:

```math
u(x, y, z) = \frac{1}{6} \left( \text{radius}^2 - \left( (x - \text{center}_x)^2 + (y - \text{center}_y)^2 + (z - \text{center}_z)^2 \right) \right)
```

**CSV source:** `results/scalar/Scalar_3D_Diffusion_Poisson_Dirichlet_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_3D_Diffusion_Poisson_Dirichlet_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```