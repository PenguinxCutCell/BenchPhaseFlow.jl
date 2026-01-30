# JohansenColella P4 Convergence

Schwartz–Colella Problem 4 (Schwartz Colella 3D Poisson):
Solve 

```math
    - \Delta \Phi = -14 f
```

inside a sphere of radius 0.392 embedded in a unit cube, with the analytical solution $\Pfi = f$, where
    
```math
    f(x,y,z) = sin(x) sin(2y) sin(3z)
```

Dirichlet boundary conditions enforce Φ = f on the immersed sphere surface and
Φ = 0 on the outer box. A mesh-convergence sweep records the errors.

**CSV source:** `results/scalar/JohansenColella_P4_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "JohansenColella_P4_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

**CSV source:** `results/scalar/JohansenColella_P4_Neumann_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "JohansenColella_P4_Neumann_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```