# Scalar 2D Diffusion Poisson DiskSteepK Parametric

Parametric study for the steep k-profile disk Poisson problem.
Sweeps eps (steepness), k_max/k_min (jump), and mean type (arithmetic/harmonic).

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Poisson_DiskSteepK_Parametric.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Poisson_DiskSteepK_Parametric.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```