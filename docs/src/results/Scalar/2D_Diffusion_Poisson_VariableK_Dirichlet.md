# Scalar 2D Diffusion Poisson DiskSteepK Convergence

2D variable-coefficient Poisson problem in a disk with a smooth steep k(r).
Solve :

```math
- \nabla \cdot (k(r) \nabla u) = f, \quad r = \lVert x - c \rVert,
```
with Dirichlet boundary conditions on the disk boundary.

Exact solution (non-radial, m=2 by default):

```math
    u(x,y) = (1 - r^2/R^2) * P_m(x,y), \quad P_m = \mathrm{Re}((x + i y)^m).
```
The conductivity profile is

```math
    k(r) = k_{\text{min}} + \frac{1}{2}(k_{\text{max}} - k_{\text{min}}) \left(1 + \tanh\left(\gamma (r - r_s)\right)\right),
```

where $r_s$ is the steep transition radius, and $\gamma$ controls the steepness.

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

## Scalar 2D Diffusion Poisson DiskSteepK Parametric Study


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