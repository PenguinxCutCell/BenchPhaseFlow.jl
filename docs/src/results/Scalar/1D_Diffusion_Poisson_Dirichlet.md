# Scalar 1D Diffusion Poisson Dirichlet Convergence

Scalar 1D diffusion–Poisson on $[0,1]$ with homogeneous Dirichlet boundaries. The cut cell is centered at $c=0.5$ with radius $r=0.11$, and the right-hand side is chosen so that

```math
u(x) = -\tfrac{1}{6}(x-c)^3 - \tfrac{c}{2}(x-c)^2 + \tfrac{r^2}{6}(x-c) + \tfrac{c r^2}{2}
```

is the exact solution. Each mesh size $h=1/n_x$ is solved with Penguin’s cut-cell Poisson operator (implicit integration capacity, Dirichlet on both ends); errors are reported in $L^2$ over the whole domain and cut/uncut partitions, with empirical orders via successive mesh-pair ratios.

**CSV source:** `results/scalar/Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence.csv`
```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_1D_Diffusion_Poisson_Dirichlet_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```
