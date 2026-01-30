# JohansenColella P5 Convergence

Johansen–Colella Problem 5: 3D heat equation inside a sphere embedded in a unit
box with the analytic solution

```math
    \Phi(x,y,z,t) = \frac{4}{5\pi (t+1)} \exp\left(-\frac{x^2 + y^2 + z^2}{5 (t+1)}\right).
``

The source term is
```math
    f(x,y,z,t) = \frac{4 (x^2 + y^2 + z^2 + 5(t+1))}{125 \pi (t+1)^3}
                 \exp\left(-\frac{x^2 + y^2 + z^2}{5 (t+1)}\right),
```

and Dirichlet boundary conditions enforce Φ on the immersed sphere. A transient
simulation to `Tend` collects the final-time error across several meshes.

**CSV source:** `results/scalar/JohansenColella_P5_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "JohansenColella_P5_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```