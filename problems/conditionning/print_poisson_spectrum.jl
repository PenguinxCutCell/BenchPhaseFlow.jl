using Penguin
using LinearAlgebra
using SparseArrays

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# Meshes and diffusivity ratios to sweep
const MESHES = [8, 16, 32]
const RATIOS = [1.0, 10.0, 100.0]

# Geometry
const LX = 1.0
const X0 = 0.0
const CENTER = 0.5
const RADIUS = 0.21
body(x, _=0) = sqrt((x - CENTER)^2) - RADIUS

# Source and boundary conditions
g(x, y, _=0) = x
const BC_DIRICHLET = Dirichlet(0.0)
const BC_B = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => BC_DIRICHLET, :bottom => BC_DIRICHLET))
const BC_I = Dirichlet(0.0)

"""
Assemble the cut-cell monophasic Poisson operator for a given resolution and diffusivity,
trim zero rows/columns, and return spectral metrics.
"""
function poisson_stats(nx::Int, diffusivity::Float64)
    mesh = Penguin.Mesh((nx,), (LX,), (X0,))
    capacity = Capacity(body, mesh)
    a = (x, y, _=0) -> diffusivity
    operator = DiffusionOps(capacity)
    phase = Phase(capacity, operator, g, a)
    solver = DiffusionSteadyMono(phase, BC_B, BC_I)

    A = solver.A
    b = solver.b
    Atrim, _, _, _ = remove_zero_rows_cols!(copy(A), copy(b))

    evals = eigvals(Matrix(Atrim))
    λmin = minimum(evals)
    λmax = maximum(evals)
    cond2 = λmax / λmin
    return (; λmin, λmax, cond2, rows=size(Atrim, 1), cols=size(Atrim, 2), nnz=nnz(Atrim))
end

function main(; csv_path=nothing)
    results_dir = isnothing(csv_path) ? joinpath(BENCH_ROOT, "results", "conditionning") : dirname(csv_path)
    mkpath(results_dir)
    outfile = isnothing(csv_path) ? joinpath(results_dir, "poisson_spectrum.csv") : csv_path
    open(outfile, "w") do io
        println(io, "scheme,nx,D,lambda_min,lambda_max,cond2,rows,cols,nnz")
        for nx in MESHES
            for D in RATIOS
                stats = poisson_stats(nx, D)
                println(io, join([
                    "poisson",
                    string(nx),
                    string(D),
                    string(stats.λmin),
                    string(stats.λmax),
                    string(stats.cond2),
                    string(stats.rows),
                    string(stats.cols),
                    string(stats.nnz)
                ], ","))
            end
        end
    end
    println("Wrote Poisson spectrum sweep to ", outfile)
end

main()
