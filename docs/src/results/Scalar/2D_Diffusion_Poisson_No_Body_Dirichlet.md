# Poisson2D MMS Convergence

Manufactured Poisson solution on the unit square $[0,1]^2$ with homogeneous Dirichlet boundaries on all sides. The exact field and forcing are

```math
u_{\text{exact}} = \sin\!\big(\pi (x-dx)/(Lx-dx)\big) \; \sin\!\big(\pi (y-dy)/(Ly-dy)\big)
```

```math
-\Delta u = f = ((\pi/(Lx-dx))^2 + (\pi/(Ly-dy))^2)\, u_{\text{exact}}
```

No Cut cells, so the table reports whole-domain errors and pairwise convergence orders for each mesh refinement $h = 1/n_x$.

**CSV source:** `results/scalar/Poisson2D_MMS_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Poisson2D_MMS_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

```
