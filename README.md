# BenchPhaseFlow.jl

BenchPhaseFlow collects the verification and regression problems that exercise [Penguin.jl](https://github.com/PenguinxCutCell/Penguin.jl). Every file inside `problems/` is a standalone Julia script that computes convergence tables and writes the resulting CSVs under `results/`.

## Getting started

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project scripts/run_all_problems.jl                 # run every benchmark
julia --project scripts/run_all_problems.jl scalar 2D        # run scripts whose path matches either filter
julia --project scripts/run_all_problems.jl --list           # print the manifest without running anything
```

The runner automatically discovers every `problems/**/*.jl` file, spawns a fresh Julia process for each script, and propagates the first non-zero exit code, which makes it safe to use in CI environments. Each script already contains lightweight `@testset`s that validate the reported convergence orders and ensure the CSV files are written.

## Directory structure

- `problems/`: benchmark definitions grouped by physics category (e.g. `scalar`, `NavierStokes`, `NavierStokesCoupled`, `conditionning`)
- `results/`: CSV convergence tables (overwritten on each run)
- `src/`, `utils/`, `test/`: helper code that backs the benchmark utilities
- `scripts/run_all_problems.jl`: convenience driver for CI and local automation

## Continuous integration

The repository ships with a GitHub Actions workflow (`.github/workflows/ci.yml`) that runs on pushes to `main` and all pull requests. The job performs:

1. `Pkg.instantiate` + `Pkg.test()` via the Julia GitHub Actions helpers
2. `julia --project scripts/run_all_problems.jl`

The workflow fails if any benchmark script raises an exception, any embedded `@testset` fails, or any problem exits with a non-zero status.

## Documentation site

The `docs/` folder contains a lightweight Documenter.jl site that scans `results/` and generates one page per CSV, reusing the descriptive prolog that appears at the top of every benchmark script. To rebuild it locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
```

By default the site reads the CSVs present in `results/`. If you want to use the read-only snapshots from the `benchmark-results` branch, check out that branch in a temporary directory and copy its `results/` folder into your working tree before running the command above.
