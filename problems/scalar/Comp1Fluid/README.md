# Comp1Fluid Diphasic Heat Cases

Three Henry-law jump variants of the 2D diphasic heat benchmark sharing a disk
of radius `R0 = 1.0` centered at `(0.125/π, 0.125*sqrt(3/2))` inside the domain
`[-20, 20] × [-20, 20]`. The final time is tied to the liquid diffusivity via
`Tend = 3.2 * Dl / R0^2 = 0.557312`. Boundary conditions are left empty to
mirror the base benchmark setup.

`He` enters the interface condition as `c_g = He * c_l` (via `ScalarJump(1, He, 0)`).

## Cases
- **air_to_water**: He = 30, `Dg/Dl = 10` with `Dl = 0.17416`, `Dg = 1.7416`.
- **water_to_air**: inverse jump and ratio with `He = 1/30`, `Dg/Dl = 0.1`
  (`Dl = 0.17416`, `Dg = 0.017416`).
- **near_unity**: He = 2.0, `Dg/Dl = 10` (`Dl = 0.17416`, `Dg = 1.7416`).

## Running
`Heat_2ph_2D_Comp1Fluid.jl` exposes helpers:
- `main(case=:air_to_water)` to run a single case (default) and write CSV under
  `results/scalar/Comp1Fluid/`.
- `run_all_cases()` to execute all three case definitions with the same mesh
  sequence.

Adjust `nx_list`, `norm`, or `relative` via keyword arguments when calling
`main` or `run_all_cases`. Initial concentrations mirror the base benchmark
(`cg0 = 1`, `cl0 = 0`).
