using Penguin
using LinearSolve
using CSV
using Test

"""
2D MMS Poisson, homogeneous Dirichlet on all borders.
u_exact = sin(pi*(x-dx)/(Lx-dx)) * sin(pi*(y-dy)/(Ly-dy))
-Laplace(u) = f = ( (pi/(Lx-dx))^2 + (pi/(Ly-dy))^2 ) * u_exact
"""

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

function run_mms_convergence_2d(nx_list; lx=1.0, ly=1.0, x0=0.0, y0=0.0,
                                method_capacity="ImplicitIntegration",
                                solver_alg=KrylovJL_GMRES(),
                                norm=2,
                                relative=false)
    h_vals = Float64[]
    err_vals = Float64[]
    err_full_vals = Float64[]
    err_cut_vals = Float64[]
    err_empty_vals = Float64[]
    inside_cells = Int[]
    inside_cells_by_dim = Vector{Vector{Int}}()

    for nx in nx_list
        ny = nx
        dx = lx / nx
        dy = ly / ny
        lx_eff = lx - dx
        ly_eff = ly - dy

        mesh = Penguin.Mesh((nx, ny), (lx, ly), (x0, y0))

        # Full domain
        body(x, y, _=0) = -1.0

        capacity = Capacity(body, mesh; method=method_capacity, compute_centroids=false)
        operator = DiffusionOps(capacity)

        # MMS exact (aligned with the 1D convention)
        u_exact(x, y) = sin(pi * (x - dx) / lx_eff) * sin(pi * (y - dy) / ly_eff)

        lambda_x = (pi / lx_eff)^2
        lambda_y = (pi / ly_eff)^2
        f(x, y, _=0) = (lambda_x + lambda_y) * u_exact(x, y)
        D(x, y, _=0) = 1.0

        # Homogeneous Dirichlet everywhere
        bc0 = Dirichlet(0.0)
        bc_b = BorderConditions(Dict(:left => bc0, :right => bc0, :bottom => bc0, :top => bc0))

        phase = Phase(capacity, operator, f, D)
        solver = DiffusionSteadyMono(phase, bc_b, bc0)

        solve_DiffusionSteadyMono!(solver; algorithm=solver_alg, log=false)

        _, _, global_err, full_err, cut_err, empty_err =
            check_convergence(u_exact, solver, capacity, norm, relative)

        push!(h_vals, 1.0 / nx)
        push!(err_vals,       global_err)
        push!(err_full_vals,  full_err)
        push!(err_cut_vals,   cut_err)
        push!(err_empty_vals, empty_err)
        push!(inside_cells, count_inside_cells(capacity))
        push!(inside_cells_by_dim, [nx, ny])
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
    nx_vals = isnothing(nx_list) ? [8, 16, 32, 64, 128] : nx_list

    data = run_mms_convergence_2d(
        nx_vals;
        lx=1.0,
        ly=1.0,
        norm=2,
        relative=false
    )

    csv_info = write_convergence_csv("Poisson2D_MMS_Convergence", data; csv_path=csv_path)

    return (data=data, csv_path=csv_info.csv_path, table=csv_info.table)
end

results = main()

@testset "Poisson 2D MMS convergence" begin
    orders = results.data.orders
    @test orders.all > 1.0
    @test length(results.data.h_vals) == length(results.data.err_vals)
    @test results.data.h_vals[1] > results.data.h_vals[end]
    @test minimum(results.data.err_vals) < maximum(results.data.err_vals)
    @test isfile(results.csv_path)
end
