# Scalar 2D Diffusion Poisson DiskSteepK Convergence

2D variable-coefficient Poisson problem in a disk with a smooth steep k(r).
Solve -div(k grad u) = f on r <= R, u = 0 on r = R.
Exact solution (non-radial, m=2 by default):
    u(x,y) = (1 - r^2/R^2) * P_m(x,y), with P_m = Re((x + i y)^m).

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Poisson_DiskSteepK_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Poisson_DiskSteepK_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```