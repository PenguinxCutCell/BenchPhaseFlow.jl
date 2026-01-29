# Scalar 2D Diffusion Poisson Robin Convergence

Steady 2D Poisson problem with an embedded circular Robin boundary.

The domain is [0, 4]^2 with a circular obstacle of radius `ly / 4` centered at
`(lx/2, ly/2) + (0.01, 0.01)`. Dirichlet zero data is imposed on the outer box
and Robin(1, 1, 1) data on the embedded boundary. The manufactured solution
`u(x, y) = 7/4 - ((x - cx)^2 + (y - cy)^2)/4` is used to measure both L2 and H1
(gradient) convergence.

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Poisson_Robin_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Poisson_Robin_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```