# Scalar 2D Diffusion Heat Moving Convergence

Scalar 2D heat equation with an oscillating circular interface. The radius
oscillates periodically, and the exact solution is imposed on the moving
boundary. This script performs a mesh convergence study and writes a CSV
summary without producing plots or timestamped folders.
# Might need to adjust interface centroid computation for moving bodies : bary_interface vs compute_interface_centroid()

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Moving_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Moving_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```