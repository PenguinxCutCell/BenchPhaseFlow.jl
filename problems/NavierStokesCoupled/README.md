# Navier-Stokes + Scalar Coupled Benchmarks

Coupled Navier–Stokes/heat regression tests with embedded `@testset`s.

- `PureConduction_PassiveCoupling.jl`: buoyancy disabled (β=0) with isothermal top/bottom; checks zero velocity and small temperature error.
- `Hydrostatic_Stratification.jl`: linear temperature stratification with buoyancy; ensures hydrostatic balance yields zero velocity and small residual.
- `DifferentiallyHeatedCavity_Ra1e3.jl`: classic Ra=1e3, Pr=0.71 cavity; asserts Nusselt and peak velocities against de Vahl Davis reference tolerances.

Run with `julia --project problems/NavierStokesCoupled/<script>.jl` or via the top-level `scripts/run_all_problems.jl`.
