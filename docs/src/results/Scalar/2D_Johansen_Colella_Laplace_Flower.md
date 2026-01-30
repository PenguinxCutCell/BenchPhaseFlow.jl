# JohansenColella P3 Diagnostics

Johansen–Colella Problem 3: Laplace equation inside a flower-shaped interface
($r = 0.25 + 0.05 \cos 6\theta$) embedded in $[0,1]^2$. Dirichlet data are $\phi = 1$ on the
immersed boundary and $\phi = 0$ on the outer box. This script mirrors the problem 3 from Johansen and Colella setup and focuses  on overshoot diagnostics, logging
statistics versus mesh size to a CSV file.


**CSV source:** `results/scalar/JohansenColella_P3_Diagnostics.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "JohansenColella_P3_Diagnostics.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```
