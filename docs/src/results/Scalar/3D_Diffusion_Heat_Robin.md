# Scalar 3D Diffusion Heat Robin Convergence 

Scalar 3D Diffusion Heat Equation with Robin Boundary Conditions
This benchmark mirrors the 2D Robin problem but inside a sphere. The initial
temperature is uniform $w_0$ and the interface satisfies :

```math
\frac{\partial w}{\partial n} + k w = 0
```

The analytical solution uses the eigenvalues μₙ obtained from 

```math
μ \cot(μ) + kR - 1 = 0:
```
```math
w(r,t) = \frac{2kR^2 w_0}{r} \sum C_n \sin\left(\frac{μ_n r}{R}\right) \exp\left(-a \frac{μ_n^2 t}{R^2}\right),
\text{where } C_n = \frac{\sin(μ_n)[μ_n^2+(kR-1)^2]}{μ_n^2[μ_n^2+kR(kR-1)]}.
```

**CSV source:** `results/scalar/Scalar_3D_Diffusion_Heat_Robin_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_3D_Diffusion_Heat_Robin_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

**CSV source:** `results/scalar/Scalar_3D_Diffusion_Heat_Robin_Convergence_Linf.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_3D_Diffusion_Heat_Robin_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```