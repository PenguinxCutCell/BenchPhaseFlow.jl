# Scalar 1D Diffusion Heat Moving ConstantBC Convergence

1D heat equation with a moving interval and constant Dirichlet boundary data.
The solution is manufactured to remain identically 1 when both the immersed
boundary and the domain edges enforce `u = 1`. This script verifies that the
numerical solution preserves this constant state while the cut location moves.

**CSV source:** `results/scalar/Scalar_1D_Diffusion_Heat_Moving_ConstantBC_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_1D_Diffusion_Heat_Moving_ConstantBC_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```