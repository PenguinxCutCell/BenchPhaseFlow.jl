# Scalar 3D Diffusion Heat Moving Convergence

Scalar 3D heat equation with an oscillating spherical interface. The radius
oscillates periodically, and the exact solution is imposed on the moving
boundary. This script performs a mesh convergence study and writes a CSV
summary without producing plots or timestamped folders.

**CSV source:** `results/scalar/Scalar_3D_Diffusion_Heat_Moving_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_3D_Diffusion_Heat_Moving_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```