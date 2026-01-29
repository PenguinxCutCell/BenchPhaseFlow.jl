# BenchPhaseFlow: Penguin.jl Benchmarks

BenchPhaseFlow hosts benchmark problems and convergence studies for the cut-cell CFD stack built on [Penguin.jl](https://github.com/PenguinxCutCell/Penguin.jl).

## What’s here
- A curated set of benchmark definitions under `problems/` (scalar diffusion, Navier–Stokes, multiphase, and conditioning studies).
- CSV outputs in `results/` with mesh sizes, error norms, and convergence rates that accompany each benchmark.
- Utility routines in `utils/convergence.jl` to post-process solver output and compute norms and rates.

## Rebuild locally
1. Run the benchmarks (produces fresh CSVs in `results/`):
   ```bash
   julia --project scripts/run_all_problems.jl
   ```
2. Build this landing page:
   ```bash
   julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
   ```
3. Open `docs/build/index.html` in your browser.

## Where to look next
- Benchmark scripts: `problems/`
- Generated CSVs: `results/`
- Post-processing helpers: `utils/convergence.jl`
