# Liu Fedkiw Test cases

Liu-Fedkiw diphasic diffusion benchmark (Case 1).

Problem: u_xx = 0 on [0, 1] with u(0) = 0, u(1) = 2, interface at x = 0.5, jump
conditions [u] = 1 and [u_x] = 0. For convergence, the analytical solution is
evaluated on the effective domain using L_eff = L - dx:
u_left(x) = (x - dx) / L_eff and u_right(x) = u_left(x) + 1.

**CSV source:** `results/scalar/diphasic/LiuFedkiw/LiuFedkiw_Case1_Convergence.csv`

```@eval
using CSV, DataFrames, Markdown, Printf

csv_path = joinpath(@__DIR__, "..", "..", "..", "..", "results", "scalar", "diphasic", "LiuFedkiw", "LiuFedkiw_Case1_Convergence.csv")
df = CSV.read(csv_path, DataFrame)

fmt(x) = x isa Float64 ? @sprintf("%.6g", x) : string(x)

header = join(string.(names(df)), " | ")
divider = join(fill("---", ncol(df)), " | ")
rows = [join([fmt(row[c]) for c in names(df)], " | ") for row in eachrow(df)]
Markdown.parse("| $header |\n| $divider |\n" * join("| " .* rows .* " |\n"))
```
