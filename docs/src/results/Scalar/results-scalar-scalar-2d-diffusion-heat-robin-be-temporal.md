# Scalar 2D Diffusion Heat Robin BE Temporal

Temporal convergence benchmark for the scalar 2D heat equation with Robin
interface conditions. The geometry, physics, and analytical reference are the
same as `Scalar_2D_Diffusion_Heat_Robin.jl`, but the mesh is fixed while the
time step varies to compare Backward Euler (BE) and Crank–Nicolson (CN).

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Robin_BE_Temporal.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Robin_BE_Temporal.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```