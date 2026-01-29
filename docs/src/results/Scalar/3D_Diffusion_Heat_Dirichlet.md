# Scalar 3D Diffusion Heat Dirichlet Convergence

Scalar 3D Diffusion Heat Equation Convergence Test
This script performs a mesh-convergence study for the 3D heat equation
using the Penguin.jl library. It mirrors the 2D benchmark but for a sphere,
reporting errors and estimated convergence orders across meshes and exporting
results to a CSV file (no plotting or ad-hoc outputs).
The analytical solution used is the classical series for a cooling sphere:

```math
T(r,t) = T_b + (T_0 - T_b) \frac{2R}{\pi r} \sum_{n=1}^{\infty} (-1)^{n+1} \frac{1}{n} \sin\!\left(\frac{n \pi r}{R}\right) e^{-\kappa \left(\frac{n \pi}{R}\right)^2 t},
```

with boundary temperature \(T_b\), initial temperature \(T_0\), radius \(R\), thermal diffusivity \(\kappa\), and radial coordinate \(r\).

**CSV source:** `results/scalar/Scalar_3D_Diffusion_Heat_Dirichlet_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_3D_Diffusion_Heat_Dirichlet_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```