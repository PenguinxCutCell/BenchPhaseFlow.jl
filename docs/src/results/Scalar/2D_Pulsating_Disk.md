# Scalar 2D Diffusion Heat Moving Convergence

**Equations**

Governing equation for the scalar field $\Phi(x,y,t)$ in the time-dependent domain $\Omega(t)$:

```math
\frac{\partial \Phi}{\partial t} = D\,\nabla^2 \Phi + S(x,y,t), \qquad (x,y) \in \Omega(t),
```

Dirichlet boundary condition prescribed on the moving interface $\Gamma(t)$ (the exact solution is imposed on the interface):

```math
\Phi(x,y,t) = \Phi_{\Gamma}(x,y,t), \qquad (x,y) \in \Gamma(t).
```

The oscillating circular interface is defined with center $c=(c_x,c_y)$ and radius

```math
R(t) = r_{\mathrm{mean}} + r_{\mathrm{amp}}\,\sin\!\left(\frac{2\pi t}{T}\right),
```

so that $\Gamma(t)=\{(x,y): r(x,y)=R(t)\}$ with $r(x,y)=\sqrt{(x-c_x)^2+(y-c_y)^2}$. The normal velocity of the interface is

```math
V_n(t)=\frac{dR}{dt}=r_{\mathrm{amp}}\frac{2\pi}{T}\,\cos\!\left(\frac{2\pi t}{T}\right).
```

Initial and exterior conditions used in the tests:

```math
\Phi(x,y,0)=\Phi_0(x,y), \qquad \Phi(x,y,t)=0 \quad\text{for } r(x,y)>R(t).
```

Parameters used in the implementation:

- diffusion coefficient: $D$,
- forcing/source term: $S(x,y,t)$,
- mean radius $r_{\mathrm{mean}}$, amplitude $r_{\mathrm{amp}}$, period $T$.

Scalar 2D heat equation with an oscillating circular interface. The radius
oscillates periodically, and the exact solution is imposed on the moving
boundary. This script performs a mesh convergence study and writes a CSV
summary without producing plots or timestamped folders.
**Might need to adjust interface centroid computation for moving bodies : bary_interface vs compute_interface_centroid()**

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Moving_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Moving_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Moving_Convergence_D1.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Moving_Convergence_D1.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

# Scalar 2D Diffusion Heat Moving CFL 2p0 Convergence

No description available; please add a docstring to the originating problem script.

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Moving_CFL_2p0_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Moving_CFL_2p0_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Moving_CFL_1p0_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Moving_CFL_1p0_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Moving_CFL_0p5_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Moving_CFL_0p5_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Moving_CFL_0p25_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Moving_CFL_0p25_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Moving_CFL_0p125_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Moving_CFL_0p125_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Moving_CFL_0p0625_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Moving_CFL_0p0625_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Moving_CFL_0p03125_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "Scalar_2D_Diffusion_Heat_Moving_CFL_0p03125_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```