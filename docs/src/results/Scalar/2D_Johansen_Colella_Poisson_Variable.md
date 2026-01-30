# JohansenColella P2 Convergence

Johansen–Colella Problem 2: variable-coefficient Poisson equation on the same
star-shaped domain as Problem 1. The diffusion coefficient is $\beta(r)=1−r^2$ and the
equation is

```math
    \nabla \cdot \big(\beta(r) \nabla \phi\big) = (7 r^2 - 15 r^4) \cos(3\theta)
```

with Dirichlet boundary condition taken from the exact solution

```math
    \phi(r, \theta) = r^4 \cos(3\theta)
```
on the boundary $\partial\Omega$.

**CSV source:** `results/scalar/JohansenColella_P2_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "JohansenColella_P2_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

**CSV source:** `results/scalar/JohansenColella_P2_Truncation_Means.csv`


```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "JohansenColella_P2_Truncation_Means.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```