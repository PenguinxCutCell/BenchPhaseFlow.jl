# Stokes 1D Mono

1D monophasic Stokes benchmark (no plotting).

This is a CSV-only version of `examples/1D/Stokes/stokes_mono.jl`. It solves a
staggered 1D Stokes system with homogeneous velocity Dirichlet boundaries and a
pinned pressure gauge, then writes key diagnostics to
`results/NavierStokes/Stokes_1D_Mono.csv`.

**CSV source:** `results/NavierStokes/Stokes_1D_Mono.csv`

| mean_uω | max_uω | mean_uγ | max_uγ | min_p | max_p | mean_p | adjoint_rel_error | viscous_sym_error | rayleigh_quotient |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 1 | 1 | -9.44269e-14 | 0 | -6.59747e-14 | 1.00058 | 0 | 6765.22 |
