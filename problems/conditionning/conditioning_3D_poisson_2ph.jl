using Penguin
using LinearAlgebra
using SparseArrays

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# 3D diphasic Poisson (steady diffusion) conditioning sweep

# Parameters to sweep
const MESHES = [(4, 4, 4), (8, 8, 8), (16, 16, 16)]  # (nx,ny,nz)
const RATIOS = [1.0, 10.0, 100.0]   # D2/D1

# Domain and body (sphere in box)
const LX, LY, LZ = 4.0, 4.0, 4.0
const X0, Y0, Z0 = 0.0, 0.0, 0.0
const RADIUS = LY / 4
const CENTER = (LX / 2, LY / 2, LZ / 2)
sphere(x, y, z) = sqrt((x - CENTER[1])^2 + (y - CENTER[2])^2 + (z - CENTER[3])^2) - RADIUS
sphere_c(x, y, z) = -sphere(x, y, z)

# Boundary and interface conditions
const BC0 = Dirichlet(0.0)
const BC_B = BorderConditions(Dict{Symbol, AbstractBoundary}(
    :left => BC0, :right => BC0, :top => BC0, :bottom => BC0, :front => BC0, :back => BC0
))
const IC = InterfaceConditions(ScalarJump(1.0, 1.0, 0.0), FluxJump(1.0, 1.0, 0.0))

# Sources
f1(x, y, z) = 0.0
f2(x, y, z) = 0.0

"""
Build cut-cell (diphasic) matrix for 3D Poisson at given resolution and ratio.
Returns trimmed A, λmin, λmax, cond2.
"""
function cutcell_stats_3d(nx::Int, ny::Int, nz::Int, ratio::Float64)
    mesh = Penguin.Mesh((nx, ny, nz), (LX, LY, LZ), (X0, Y0, Z0))
    capacity = Capacity(sphere, mesh)
    capacity_c = Capacity(sphere_c, mesh)

    operator = DiffusionOps(capacity)
    operator_c = DiffusionOps(capacity_c)

    D1val = 1.0
    D2val = ratio
    D1 = (x, y, z) -> D1val
    D2 = (x, y, z) -> D2val

    phase1 = Phase(capacity, operator, f1, D1)
    phase2 = Phase(capacity_c, operator_c, f2, D2)

    solver = DiffusionSteadyDiph(phase1, phase2, BC_B, IC)

    A = solver.A
    b = solver.b
    Atrim, _, _, _ = remove_zero_rows_cols!(copy(A), copy(b))

    evals = eigvals(Matrix(Atrim))
    evals = real.(evals)
    λmin = minimum(evals)
    λmax = maximum(evals)
    cond2 = λmax / λmin
    return (; λmin, λmax, cond2, rows=size(Atrim, 1), cols=size(Atrim, 2), nnz=nnz(Atrim))
end

"""
Assemble a standard 7-point FV Poisson matrix on a uniform grid with Dirichlet=0 on all faces.
Piecewise constant D: D=D1 inside sphere, D=D2 outside. Interface handled via harmonic averaging
between neighboring cell-centered diffusivities.
Returns A, λmin, λmax, cond2.
"""
function fv3d_stats(nx::Int, ny::Int, nz::Int, ratio::Float64)
    Lx, Ly, Lz = LX, LY, LZ
    dx, dy, dz = Lx / nx, Ly / ny, Lz / nz

    lin(i, j, k) = (k - 1) * nx * ny + (j - 1) * nx + i
    N = nx * ny * nz

    xc(i) = X0 + (i - 0.5) * dx
    yc(j) = Y0 + (j - 0.5) * dy
    zc(k) = Z0 + (k - 0.5) * dz

    D1, D2 = 1.0, ratio
    Dcell = Array{Float64}(undef, nx, ny, nz)
    for k in 1:nz, j in 1:ny, i in 1:nx
        x, y, z = xc(i), yc(j), zc(k)
        Dcell[i, j, k] = sphere(x, y, z) <= 0 ? D1 : D2
    end

    harm(a, b) = 2.0 / (1.0 / a + 1.0 / b)

    A = spzeros(Float64, N, N)
    for k in 1:nz, j in 1:ny, i in 1:nx
        p = lin(i, j, k)
        diag_val = 0.0

        # x- neighbor
        if i > 1
            dface = harm(Dcell[i-1, j, k], Dcell[i, j, k])
            T = dface / dx^2
            q = lin(i - 1, j, k)
            A[p, q] -= T
            A[q, p] -= T
            diag_val += T
        else
            dface = Dcell[i, j, k]
            T = dface / dx^2
            diag_val += T
        end

        # x+ neighbor
        if i < nx
            dface = harm(Dcell[i, j, k], Dcell[i+1, j, k])
            T = dface / dx^2
            q = lin(i + 1, j, k)
            A[p, q] -= T
            A[q, p] -= T
            diag_val += T
        else
            dface = Dcell[i, j, k]
            T = dface / dx^2
            diag_val += T
        end

        # y- neighbor
        if j > 1
            dface = harm(Dcell[i, j-1, k], Dcell[i, j, k])
            T = dface / dy^2
            q = lin(i, j - 1, k)
            A[p, q] -= T
            A[q, p] -= T
            diag_val += T
        else
            dface = Dcell[i, j, k]
            T = dface / dy^2
            diag_val += T
        end

        # y+ neighbor
        if j < ny
            dface = harm(Dcell[i, j, k], Dcell[i, j+1, k])
            T = dface / dy^2
            q = lin(i, j + 1, k)
            A[p, q] -= T
            A[q, p] -= T
            diag_val += T
        else
            dface = Dcell[i, j, k]
            T = dface / dy^2
            diag_val += T
        end

        # z- neighbor
        if k > 1
            dface = harm(Dcell[i, j, k-1], Dcell[i, j, k])
            T = dface / dz^2
            q = lin(i, j, k - 1)
            A[p, q] -= T
            A[q, p] -= T
            diag_val += T
        else
            dface = Dcell[i, j, k]
            T = dface / dz^2
            diag_val += T
        end

        # z+ neighbor
        if k < nz
            dface = harm(Dcell[i, j, k], Dcell[i, j, k+1])
            T = dface / dz^2
            q = lin(i, j, k + 1)
            A[p, q] -= T
            A[q, p] -= T
            diag_val += T
        else
            dface = Dcell[i, j, k]
            T = dface / dz^2
            diag_val += T
        end

        A[p, p] += diag_val
    end

    evals = eigvals(Matrix(A))
    evals = real.(evals)
    λmin = minimum(evals)
    λmax = maximum(evals)
    cond2 = λmax / λmin
    return (; λmin, λmax, cond2, rows=size(A, 1), cols=size(A, 2), nnz=nnz(A))
end

function main(; csv_path=nothing)
    results_dir = isnothing(csv_path) ? joinpath(BENCH_ROOT, "results", "conditionning") : dirname(csv_path)
    mkpath(results_dir)
    outpath = isnothing(csv_path) ? joinpath(results_dir, "conditioning_3D_poisson_2ph.csv") : csv_path
    open(outpath, "w") do io
        println(io, "scheme,nx,ny,nz,ratio,lambda_min,lambda_max,cond2,rows,cols,nnz")
        for (nx, ny, nz) in MESHES
            for ratio in RATIOS
                stats_cut = cutcell_stats_3d(nx, ny, nz, ratio)
                println(io, join([
                    "cutcell",
                    string(nx), string(ny), string(nz), string(ratio),
                    string(stats_cut.λmin), string(stats_cut.λmax), string(stats_cut.cond2),
                    string(stats_cut.rows), string(stats_cut.cols), string(stats_cut.nnz)
                ], ","))

                stats_fv = fv3d_stats(nx, ny, nz, ratio)
                println(io, join([
                    "fv",
                    string(nx), string(ny), string(nz), string(ratio),
                    string(stats_fv.λmin), string(stats_fv.λmax), string(stats_fv.cond2),
                    string(stats_fv.rows), string(stats_fv.cols), string(stats_fv.nnz)
                ], ","))
            end
        end
    end
    println("CSV written to: ", outpath)
end

main()
