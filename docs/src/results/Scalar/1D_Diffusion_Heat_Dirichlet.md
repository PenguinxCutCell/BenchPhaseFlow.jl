# Scalar 1D Diffusion (Dirichlet)

This test exercises the 1D heat equation on the interval $[0.25, 0.75]$ with homogeneous Dirichlet boundaries, evolved to $t=0.1$ using Penguin.jl’s cut-cell stack. The reference field is the classical odd-sine series:

```math
u(x,t) = \frac{4}{\pi} \sum_{n=0}^{\infty} \frac{1}{2n+1}\;
\sin\!\left(\frac{(2n+1)\pi\,(x - (c - R))}{2R}\right)
\exp\!\left(-\kappa \Big(\frac{(2n+1)\pi}{2R}\Big)^2 t\right),
```

with center $c=0.5$, half-length $R=0.25$, and diffusivity $\kappa=1$. A CN stepper with $\Delta t = 0.5\,\Delta x^2$ advances the solution; errors are measured in $L^2$ against the truncated series (400 terms).

**CSV source:** `results/scalar/Scalar_1D_Diffusion_Heat_Dirichlet_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_1D_Diffusion_Heat_Dirichlet_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```
