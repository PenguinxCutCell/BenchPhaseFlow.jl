# Scalar 2D Diffusion Heat NeumannMass Mass

2D transient diffusion inside a circle with homogeneous Neumann conditions on
the immersed boundary (cut cells) and on the outer box. The initial field is
identically 1 and the source term is zero, so the total integral of the
solution should remain constant in time. This script marches the solution with
Backward Euler, records the volume integral, and reports the drift together
with the unexpectedly large absolute values produced by the null-space of the
pure Neumann operator.

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_NeumannMass_Mass.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_NeumannMass_Mass.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```