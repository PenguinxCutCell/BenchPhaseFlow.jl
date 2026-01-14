using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using CSV
using Test

"""
1D MMS Poisson with non-homogeneous Dirichlet boundary conditions.
Analytical solution:
    u(x) = u_left + (u_right - u_left) * (x-dx)/(lx-dx) + A * sin(pi*(x-dx)/(lx-dx))
"""

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

function get_analytical_solution(dx, lx, u_left, u_right, A)
    u_linear = (x) -> u_left + (u_right - u_left) * (x - dx) / (lx - dx)
    u_sin = (x) -> sin(pi * (x - dx) / (lx - dx))
    return (x) -> u_linear(x) + A * u_sin(x)
end

function get_source_term(dx, lx, A)
    lambda = (pi / (lx - dx))^2
    u_sin = (x) -> sin(pi * (x - dx) / (lx - dx))
    return (x, _y=0, _=0) -> A * lambda * u_sin(x)
end

function run_mesh_convergence_mms_1d(
    nx_list::Vector{Int};
    lx::Float64 = 1.0,
    x0::Float64 = 0.0,
    u_left::Float64 = 1.0,
    u_right::Float64 = 2.0,
    amplitude::Float64 = 0.5,
    norm = 2,
    relative::Bool = false
)
    h_vals = Float64[]
    err_vals = Float64[]
    err_full_vals = Float64[]
    err_cut_vals = Float64[]
    err_empty_vals = Float64[]
    inside_cells = Int[]
    inside_cells_by_dim = Vector{Vector{Int}}()

    for nx in nx_list
        dx = lx / nx
        mesh = Penguin.Mesh((nx,), (lx,), (x0,))

        body = (x, _=0) -> -1.0
        capacity = Capacity(body, mesh; method="ImplicitIntegration", compute_centroids=false)
        operator = DiffusionOps(capacity)

        u_exact = get_analytical_solution(dx, lx, u_left, u_right, amplitude)
        f = get_source_term(dx, lx, amplitude)
        D = (x, _y=0, _=0) -> 1.0

        bc = Dirichlet(u_left)
        bc_b = BorderConditions(Dict(
            :bottom => Dirichlet(u_left),
            :top => Dirichlet(u_right)
        ))

        phase = Phase(capacity, operator, f, D)
        solver = DiffusionSteadyMono(phase, bc_b, bc)
        solve_DiffusionSteadyMono!(solver; method=Base.:\)

        _, _, global_err, full_err, cut_err, empty_err =
            check_convergence(u_exact, solver, capacity, norm, relative)

        push!(h_vals, dx)
        push!(err_vals,       global_err)
        push!(err_full_vals,  full_err)
        push!(err_cut_vals,   cut_err)
        push!(err_empty_vals, empty_err)
        push!(inside_cells, count_inside_cells(capacity))
        push!(inside_cells_by_dim, [nx])
    end

    return (
        h_vals = h_vals,
        err_vals = err_vals,
        err_full_vals = err_full_vals,
        err_cut_vals = err_cut_vals,
        err_empty_vals = err_empty_vals,
        inside_cells = inside_cells,
        inside_cells_by_dim = inside_cells_by_dim,
        orders = compute_orders(h_vals, err_vals, err_full_vals, err_cut_vals),
        norm = norm
    )
end

function write_convergence_csv(method_name, data; csv_path=nothing)
    df = make_convergence_dataframe(method_name, data)

    results_dir = csv_path === nothing ?
        joinpath(BENCH_ROOT, "results", "scalar") :
        dirname(csv_path)

    mkpath(results_dir)
    csv_out = csv_path === nothing ?
        joinpath(results_dir, "$(method_name).csv") :
        csv_path

    CSV.write(csv_out, df)
    return (csv_path = csv_out, table = df)
end

function main(; csv_path=nothing, nx_list=nothing)
    nx_vals = isnothing(nx_list) ? [5, 10, 20, 40, 80, 160] : nx_list

    data = run_mesh_convergence_mms_1d(
        nx_vals;
        lx=1.0,
        x0=0.0,
        u_left=1.0,
        u_right=2.0,
        amplitude=0.5,
        norm=2,
        relative=false
    )

    csv_info = write_convergence_csv("Scalar_1D_Diffusion_Poisson_MMS_NonHomogeneous_Convergence", data; csv_path=csv_path)

    return (data=data, csv_path=csv_info.csv_path, table=csv_info.table)
end

results = main()

@testset "Poisson 1D MMS non-homogeneous convergence" begin
    orders = results.data.orders
    @test orders.all > 1.0
    @test length(results.data.h_vals) == length(results.data.err_vals)
    @test results.data.h_vals[1] > results.data.h_vals[end]
    @test minimum(results.data.err_vals) < maximum(results.data.err_vals)
    @test isfile(results.csv_path)
end
