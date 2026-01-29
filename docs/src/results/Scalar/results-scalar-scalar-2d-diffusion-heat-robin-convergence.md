# Scalar 2D Diffusion Heat Robin Convergence

Benchmark: Scalar 2D Heat Equation with Robin Interface (legacy Heat.jl port)
This script reproduces the benchmark/Heat.jl workflow but in the BenchPhaseFlow
structure: it performs a mesh-convergence study for the 2D transient diffusion
problem around a heated circle with Robin interface conditions, logs errors, and
writes a single CSV summary (no plots or timestamped folders).

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Robin_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Robin_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```