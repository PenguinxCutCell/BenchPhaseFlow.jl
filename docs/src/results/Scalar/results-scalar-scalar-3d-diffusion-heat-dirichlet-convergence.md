# Scalar 3D Diffusion Heat Dirichlet Convergence

Scalar 3D Diffusion Heat Equation Convergence Test
This script performs a mesh-convergence study for the 3D heat equation
using the Penguin.jl library. It mirrors the 2D benchmark but for a sphere,
reporting errors and estimated convergence orders across meshes and exporting
results to a CSV file (no plotting or ad-hoc outputs).
The analytical solution used is the classical series for a cooling sphere:
u(r,t) = Tb + (T0 - Tb) * (2R)/(πr) * Σ_{n=1}^∞ (-1)^{n+1} (1/n) sin(nπr/R) exp(-κ (nπ/R)^2 t),
with the r → 0 limit handled analytically.

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