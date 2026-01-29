# Scalar 1D Diffusion Poisson MMS NonHomogeneous Convergence

1D MMS Poisson with non-homogeneous Dirichlet boundary conditions.
Analytical solution:

```math
u(x) = u_{\text{left}} + (u_{\text{right}} - u_{\text{left}}) \; \frac{x-dx}{lx-dx} + A \; \sin\!\big(\pi (x-dx)/(lx-dx)\big)
```

**CSV source:** `results/scalar/Scalar_1D_Diffusion_Poisson_MMS_NonHomogeneous_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_1D_Diffusion_Poisson_MMS_NonHomogeneous_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

