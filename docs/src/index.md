# BenchPhaseFlow Results

This micro-site is generated with [Documenter.jl](https://documenter.julialang.org/) and summarizes the CSV convergence studies stored under `results/`. Each benchmark script inside `problems/` writes a CSV and ships with a short description at the top of the file. The build step reuses that text so you can quickly remember the mathematical setup behind each dataset.

## Updating the site

1. Make sure the latest CSVs are available locally (for CI the `benchmark-results` branch can be fetched and copied into `results/`).
2. Run the benchmark suite if needed:
   ```bash
   julia --project scripts/run_all_problems.jl
   ```
3. Build the documentation:
   ```bash
   julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
   ```
4. The HTML output lives in `docs/build`. Serve it locally with any static file server or publish it via GitHub Pages/Documenter deploys.

The following sections are generated automatically—each CSV becomes its own page containing the textual problem statement and a Markdown preview of the data table. No manual editing is required.
