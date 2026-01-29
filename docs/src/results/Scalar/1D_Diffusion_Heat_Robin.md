# Scalar 1D Diffusion (Robin)

This benchmark mirrors `examples/1D/Diffusion/Heat_robin.jl`: a 1D transient
diffusion problem on $x \in [0, 10]$ with interface position $\text{center}=0.25$,
diffusivity $a=5$, and homogeneous Robin condition $u_x + u = 0$ at the cut.
The exact solution used for convergence is a complementary error-function
profile:

```math
u(x,t) = \operatorname{erf}(\eta)
	+ \exp\big(k(x-\text{center}) + a k^2 t\big)\; \operatorname{erfc}(\eta + k \sqrt{a t}),
\qquad \eta = \frac{x - \text{center}}{2 \sqrt{a t}},\; k = 1.
```

This satisfies the same boundary data (Dirichlet $u=1$ on the left, $u=0$ on
the right) and the source-free diffusion equation $u_t = a\,u_{xx}$.

**CSV source:** `results/scalar/Scalar_1D_Diffusion_Heat_Robin_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_1D_Diffusion_Heat_Robin_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

**CSV source:** `results/scalar/Scalar_1D_Diffusion_Heat_Robin_Linf_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_1D_Diffusion_Heat_Robin_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```