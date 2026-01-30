# JohansenColella P1 Convergence

Johansen–Colella Problem 1: constant-coefficient Poisson equation inside a
star-shaped domain defined by

```math
    \Omega = \{(r, \theta) : r \leq 0.30 + 0.15 \cos(6\theta)\}.
```
We solve 

```math
    \Delta \phi = 7 r^2 \cos(3\theta)
```

with Dirichlet boundary condition taken from the exact solution

```math
    \phi(r, \theta) = r^4 \cos(3\theta)
```
on the boundary $\partial\Omega$. Convergence is measured on uniform Cartesian meshes using
Penguin's implicit-integration capacity.

Truncation errors are computed by inserting the exact solution sampled at the cut-cell centroids into the discrete operator and comparing it with the continuous operator evaluated at the same centroid locations. 
This differs slightly from the convention used by Johansen and Colella, who evaluate the exact solution at Cartesian cell centers while evaluating the right-hand side at cut-cell centroids. 
The error constants might differ.

**CSV source:** `results/scalar/JohansenColella_P1_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "JohansenColella_P1_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

**CSV source:** `results/scalar/JohansenColella_P1_1997_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "JohansenColella_P1_1997_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

**CSV source:** `results/scalar/JohansenColella_P1_InfNorm_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "JohansenColella_P1_InfNorm_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```
