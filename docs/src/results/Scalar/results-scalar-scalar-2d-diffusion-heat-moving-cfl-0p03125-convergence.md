# Scalar 2D Diffusion Heat Moving CFL 0p03125 Convergence

No description available; please add a docstring to the originating problem script.

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Moving_CFL_0p03125_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Moving_CFL_0p03125_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```