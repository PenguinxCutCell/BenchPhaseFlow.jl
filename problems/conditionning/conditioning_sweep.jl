using Penguin
using LinearAlgebra
using SparseArrays

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# Sweep parameters
const RATIOS = [1.0, 10.0, 100.0]   # D2/D1
const MESHES = [8, 16, 32]          # number of cells (nx)

# Problem setup
const LX = 8.0
const X0 = 0.0
const XINT = 4.05

# Interface body definitions
body(x, _=0) = (x - XINT)
body_c(x, _=0) = -(x - XINT)

# Boundary and interface conditions
const BC1 = Dirichlet(0.0)
const BC0 = Dirichlet(1.0)
const BC_B = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => BC0, :bottom => BC1))
const IC = InterfaceConditions(ScalarJump(1.0, 1.0, 0.0), FluxJump(1.0, 1.0, 0.0))

# Zero source terms
f1(x, y, z, t) = 0.0
f2(x, y, z, t) = 0.0

"""
Build the cut-cell (diphasic) system matrix A for given nx and diffusivity ratio r = D2/D1.
Returns a trimmed sparse A (zero rows/cols removed), along with λmin, λmax, cond2.
"""
function cutcell_stats(nx::Int, ratio::Float64)
    mesh = Penguin.Mesh((nx,), (LX,), (X0,))
    capacity = Capacity(body, mesh)
    capacity_c = Capacity(body_c, mesh)

    operator = DiffusionOps(capacity)
    operator_c = DiffusionOps(capacity_c)

    D1val = 1.0
    D2val = ratio
    D1 = (x, y, z) -> D1val
    D2 = (x, y, z) -> D2val

    fluide_1 = Phase(capacity, operator, f1, D1)
    fluide_2 = Phase(capacity_c, operator_c, f2, D2)

    # Initial condition vectors (as in example)
    u0ₒ1 = zeros(nx + 1)
    u0ᵧ1 = zeros(nx + 1)
    u0ₒ2 = ones(nx + 1)
    u0ᵧ2 = ones(nx + 1)
    u0 = vcat(u0ₒ1, u0ᵧ1, u0ₒ2, u0ᵧ2)

    # Time stepping (small dt, one step is enough to assemble A)
    Δt = 0.5 * (LX / nx)^2
    Tend = Δt
    solver = DiffusionUnsteadyDiph(fluide_1, fluide_2, BC_B, IC, Δt, u0, "CN")

    # Run a single step to ensure A is built
    solve_DiffusionUnsteadyDiph!(
        solver, fluide_1, fluide_2,
        Δt, Tend, BC_B, IC, "CN";
        method = Base.:\
    )

    # Extract and trim A
    A = solver.A
    b = solver.b
    Atrim, _, _, _ = remove_zero_rows_cols!(copy(A), copy(b))

    # Convert to dense symmetric for eigen analysis
    evals = eigvals(Matrix(Atrim))
    evals = real.(evals)
    λmin = minimum(evals)
    λmax = maximum(evals)
    cond2 = λmax / λmin
    return (; λmin, λmax, cond2, rows=size(Atrim, 1), cols=size(Atrim, 2), nnz=nnz(Atrim))
end

"""
Build a standard 1D finite-volume diffusion matrix (Dirichlet at x=0 and x=L),
with piecewise constant D: D=D1 for x<xint, D=D2 for x>xint, uniform grid nx.
Returns dense symmetric A, along with λmin, λmax, cond2.
"""
function fv_stats(nx::Int, ratio::Float64)
    D1 = 1.0
    D2 = ratio
    L = LX
    dx = L / nx
    N = nx

    # Face transmissibilities T_{i+1/2} = D_face/dx
    faces = collect(0:N) # index i for face i+1/2
    T = zeros(Float64, N + 1)
    x_face(i) = (i + 0.5) * dx
    for i in faces
        xF = x_face(i)
        if isapprox(xF, XINT; atol=1e-12)
            Df = 2.0 / (1.0 / D1 + 1.0 / D2)
            T[i + 1] = Df / dx
        elseif xF < XINT
            T[i + 1] = D1 / dx
        else
            T[i + 1] = D2 / dx
        end
    end

    # Assemble symmetric tridiagonal matrix A (Dirichlet BCs incorporated on diagonal)
    A = zeros(Float64, N, N)
    for i in 1:N
        Tm = T[i]
        Tp = T[i + 1]
        A[i, i] = Tm + Tp
        if i > 1
            A[i, i-1] = -Tm
            A[i-1, i] = -Tm
        end
        if i < N
            A[i, i+1] = -Tp
            A[i+1, i] = -Tp
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
    outpath = isnothing(csv_path) ? joinpath(results_dir, "conditioning_sweep.csv") : csv_path
    open(outpath, "w") do io
        println(io, "scheme,nx,ratio,lambda_min,lambda_max,cond2,rows,cols,nnz")
        for nx in MESHES
            for ratio in RATIOS
                stats_cut = cutcell_stats(nx, ratio)
                println(io, join([
                    "cutcell",
                    string(nx),
                    string(ratio),
                    string(stats_cut.λmin),
                    string(stats_cut.λmax),
                    string(stats_cut.cond2),
                    string(stats_cut.rows),
                    string(stats_cut.cols),
                    string(stats_cut.nnz)
                ], ","))

                stats_fv = fv_stats(nx, ratio)
                println(io, join([
                    "fv",
                    string(nx),
                    string(ratio),
                    string(stats_fv.λmin),
                    string(stats_fv.λmax),
                    string(stats_fv.cond2),
                    string(stats_fv.rows),
                    string(stats_fv.cols),
                    string(stats_fv.nnz)
                ], ","))
            end
        end
    end

    println("CSV written to: ", outpath)
end

main()
