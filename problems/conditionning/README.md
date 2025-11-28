# Conditioning benchmarks

CSV-only versions of the original conditioning scripts. Each script writes its output next to the file.

- `print_poisson_spectrum.jl`: 1D cut-cell Poisson spectrum sweep over mesh size and diffusivity ratio.
- `conditioning_poisson.jl`: compares eigenvalue extremes and condition numbers for cut-cell vs. standard FV Poisson in 1D.
- `conditioning_sweep.jl`: 1D two-phase diffusion conditioning sweep across mesh sizes and diffusivity ratios.
- `conditioning_3D_poisson_2ph.jl`: 3D two-phase Poisson conditioning sweep (sphere in a box) for small meshes.

Run with `julia --project problems/conditionning/<script>.jl`; outputs default to `results/conditionning/*.csv` under the repo root. Override with `csv_path=<path>` when calling `main` programmatically.
