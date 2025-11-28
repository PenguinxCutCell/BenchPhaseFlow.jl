using Penguin
using LinearAlgebra
using SparseArrays

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# Mesh sizes to sweep
const MESHES = [8, 16, 32]

# Domain and geometry
const LX = 1.0
const X0 = 0.0
const CENTER = 0.5
const RADIUS = 0.21
body(x, _=0) = sqrt((x - CENTER)^2) - RADIUS

# Boundary conditions: Dirichlet 0 on both ends and on the embedded boundary
const BC_LEFT = Dirichlet(0.0)
const BC_RIGHT = Dirichlet(0.0)
const BC_B = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => Dirichlet(0.0), :bottom => Dirichlet(0.0)))
const BC_I = Dirichlet(0.0)

# Source and diffusion coefficient (constant)
g(x, y, _=0) = x
a(x, y, _=0) = 1.0

"""
Build the cut-cell monophasic steady diffusion (Poisson) matrix for nx,
trim zero rows/cols, and return λmin, λmax, and cond2.
"""
function cutcell_stats_poisson(nx::Int)
    mesh = Penguin.Mesh((nx,), (LX,), (X0,))
    capacity = Capacity(body, mesh)
    operator = DiffusionOps(capacity)
    fluide = Phase(capacity, operator, g, a)

    solver = DiffusionSteadyMono(fluide, BC_B, BC_I)

    A = solver.A
    b = solver.b
    Atrim, _, _, _ = remove_zero_rows_cols!(copy(A), copy(b))

    evals = eigvals(Matrix(Atrim))
    λmin = minimum(evals)
    λmax = maximum(evals)
    cond2 = λmax / λmin
    return (; λmin, λmax, cond2, rows=size(Atrim, 1), cols=size(Atrim, 2), nnz=nnz(Atrim))
end

"""
Standard 1D finite-volume Poisson matrix on [0,L] with Dirichlet 0 at both ends,
constant diffusion a=1. Returns eigen extrema and cond2.
"""
function fv_stats_poisson(nx::Int)
    L = LX
    dx = L / nx
    N = nx
    T = 1.0 / dx  # face transmissibility with D=1

    A = zeros(Float64, N, N)
    for i in 1:N
        A[i, i] = 2T
        if i > 1
            A[i, i-1] = -T
            A[i-1, i] = -T
        end
        if i < N
            A[i, i+1] = -T
            A[i+1, i] = -T
        end
    end

    evals = eigvals(Symmetric(A))
    λmin = minimum(evals)
    λmax = maximum(evals)
    cond2 = λmax / λmin
    return (; λmin, λmax, cond2, rows=size(A, 1), cols=size(A, 2), nnz=nnz(sparse(A)))
end

function main(; csv_path=nothing)
    results_dir = isnothing(csv_path) ? joinpath(BENCH_ROOT, "results", "conditionning") : dirname(csv_path)
    mkpath(results_dir)
    outpath = isnothing(csv_path) ? joinpath(results_dir, "conditioning_poisson.csv") : csv_path
    open(outpath, "w") do io
        println(io, "scheme,nx,lambda_min,lambda_max,cond2,rows,cols,nnz")
        for nx in MESHES
            stats_cut = cutcell_stats_poisson(nx)
            println(io, join([
                "cutcell",
                string(nx),
                string(stats_cut.λmin),
                string(stats_cut.λmax),
                string(stats_cut.cond2),
                string(stats_cut.rows),
                string(stats_cut.cols),
                string(stats_cut.nnz)
            ], ","))

            stats_fv = fv_stats_poisson(nx)
            println(io, join([
                "fv",
                string(nx),
                string(stats_fv.λmin),
                string(stats_fv.λmax),
                string(stats_fv.cond2),
                string(stats_fv.rows),
                string(stats_fv.cols),
                string(stats_fv.nnz)
            ], ","))
        end
    end

    println("CSV written to: ", outpath)
end

main()
