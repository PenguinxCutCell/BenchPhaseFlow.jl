using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using CSV
using Test

"""
2D variable-coefficient Poisson problem in a disk with a smooth steep k(r).
Solve -div(k grad u) = f on r <= R, u = 0 on r = R.
Exact solution (non-radial, m=2 by default):
    u(x,y) = (1 - r^2/R^2) * P_m(x,y), with P_m = Re((x + i y)^m).
TODO : Do a parametric study on k-profile steepness (eps), jump (k_max/k_min), and arithmetic vs harmonic mean for diffusivity in cut cells.
"""

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

function k_profile(r; k_min=1.0, k_max=10.0, r0=0.25, eps=0.025)
    delta_k = k_max - k_min
    return k_min + 0.5 * delta_k * (1.0 + tanh((r - r0) / eps))
end

function k_profile_prime(r; k_min=1.0, k_max=10.0, r0=0.25, eps=0.025)
    delta_k = k_max - k_min
    z = (r - r0) / eps
    sech2 = 1.0 / cosh(z)^2
    return 0.5 * delta_k / eps * sech2
end

function poly_real_part(x, y, m)
    z = complex(x, y)
    return real(z^m)
end

function exact_solution(x, y, r, R, m, center)
    xr = x - center[1]
    yr = y - center[2]
    return (1.0 - (r / R)^2) * poly_real_part(xr, yr, m)
end

function forcing_term(x, y, r, R, m, center; k_min=1.0, k_max=10.0, r0=0.25, eps=0.025)
    k_val = k_profile(r; k_min=k_min, k_max=k_max, r0=r0, eps=eps)
    k_prime = k_profile_prime(r; k_min=k_min, k_max=k_max, r0=r0, eps=eps)
    if r == 0.0
        return 0.0
    end

    xr = x - center[1]
    yr = y - center[2]
    p_m = poly_real_part(xr, yr, m)
    a_prime = m * r^(m - 1) - (m + 2) * r^(m + 1) / (R^2)
    return (4.0 * (m + 1) / (R^2)) * k_val * p_m - (k_prime * a_prime) * (p_m / r^m)
end

function run_disk_convergence(
    nx_list::Vector{Int};
    lx::Float64 = 1.0,
    ly::Float64 = 1.0,
    center::Tuple{Float64,Float64} = (0.5, 0.5),
    R::Float64 = 0.5,
    k_min::Float64 = 1.0,
    k_max::Float64 = 10.0,
    r0::Float64 = 0.25,
    eps::Float64 = 0.025,
    m::Int = 2,
    norm::Real = 2,
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
        ny = nx
        mesh = Penguin.Mesh((nx, ny), (lx, ly), (0.0, 0.0))
        r_of_xy = (x, y) -> sqrt((x - center[1])^2 + (y - center[2])^2)
        body = (x, y, _=0) -> r_of_xy(x, y)- R
        capacity = Capacity(body, mesh; method="ImplicitIntegration")
        operator = DiffusionOps(capacity)

        u_exact = (x, y) -> exact_solution(x, y, r_of_xy(x, y), R, m, center)
        f = (x, y, _=0) -> forcing_term(x, y, r_of_xy(x, y), R, m, center; k_min=k_min, k_max=k_max, r0=r0, eps=eps)
        D = (x, y, _=0) -> k_profile(r_of_xy(x, y); k_min=k_min, k_max=k_max, r0=r0, eps=eps)

        bc_outer = Dirichlet(0.0)
        bc_b = BorderConditions(Dict(
            :left   => bc_outer,
            :right  => bc_outer,
            :top    => bc_outer,
            :bottom => bc_outer
        ))

        phase = Phase(capacity, operator, f, D)
        solver = DiffusionSteadyMonoVariable(phase, bc_b, Dirichlet((x,y,_=0)->u_exact(x,y)))
        solve_DiffusionSteadyMono!(solver; method=Base.:\)

        _, _, global_err, full_err, cut_err, empty_err =
            check_convergence(u_exact, solver, capacity, norm, relative)

        push!(h_vals, min(lx / nx, ly / ny))
        push!(err_vals, global_err)
        push!(err_full_vals, full_err)
        push!(err_cut_vals, cut_err)
        push!(err_empty_vals, empty_err)
        push!(inside_cells, count_inside_cells(capacity))
        dx = lx / nx
        dy = ly / ny
        coverage_x = ceil(Int, 2 * R / dx)
        coverage_y = ceil(Int, 2 * R / dy)
        push!(inside_cells_by_dim, [coverage_x, coverage_y])
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
    results_dir = isnothing(csv_path) ? joinpath(BENCH_ROOT, "results", "scalar") : dirname(csv_path)
    mkpath(results_dir)
    csv_out = isnothing(csv_path) ? joinpath(results_dir, "$(method_name)_Convergence.csv") : csv_path
    CSV.write(csv_out, df)
    return (csv_path = csv_out, table = df)
end

function main(; csv_path=nothing, nx_list=nothing)
    nx_vals = isnothing(nx_list) ? [40, 80, 160, 320] : nx_list
    data = run_disk_convergence(
        nx_vals;
        lx=1.0,
        ly=1.0,
        center=(0.5, 0.5),
        R=0.5,
        k_min=1.0,
        k_max=10.0,
        r0=0.25,
        eps=0.025,
        m=2,
        norm=2,
        relative=false
    )

    csv_info = write_convergence_csv("Scalar_2D_Diffusion_Poisson_DiskSteepK", data; csv_path=csv_path)
    return (data = data, csv_path = csv_info.csv_path, table = csv_info.table)
end

results = main()

@testset "Poisson disk steep-k convergence" begin
    orders = results.data.orders
    @test !isnan(orders.all)
    @test length(results.data.h_vals) == length(results.data.err_vals)
    @test results.data.h_vals[1] > results.data.h_vals[end]
    @test minimum(results.data.err_vals) < maximum(results.data.err_vals)
    @test isfile(results.csv_path)
end

using CairoMakie
function plot_solutionn(solver, mesh, capacity, u_exact)
    u_ana, u_num, _, _, _, _ = check_convergence(u_exact, solver, capacity)

    cell_types = capacity.cell_types
    u_num_promasked = copy(u_num)
    u_ana_promasked = copy(u_ana)
    u_num_promasked[cell_types .== 0] .= NaN
    u_ana_promasked[cell_types .== 0] .= NaN

    xg = mesh.nodes[1]
    yg = mesh.nodes[2]
    nx_plot = length(xg)
    ny_plot = length(yg)

    sol_grid = reshape(u_num_promasked, (nx_plot, ny_plot))'
    ana_grid = reshape(u_ana_promasked, (nx_plot, ny_plot))'
    err_grid = sol_grid .- ana_grid

    fig_sol = Figure(resolution=(500,400))
    ax1 = Axis(fig_sol[1,1], title="Solution (numerical)", xlabel="x", ylabel="y", aspect=DataAspect())
    hm1 = heatmap!(ax1, xg, yg, sol_grid, colormap = :viridis)
    Colorbar(fig_sol[1,2], hm1, label = "Value")
    display(fig_sol)
    save("disk_steepk_solution.png", fig_sol)

    fig_ana = Figure(resolution=(500,400))
    ax2 = Axis(fig_ana[1,1], title="Solution (analytical)", xlabel="x", ylabel="y", aspect=DataAspect())
    hm2 = heatmap!(ax2, xg, yg, ana_grid, colormap = :viridis)
    Colorbar(fig_ana[1,2], hm2, label = "Value")
    display(fig_ana)
    save("disk_steepk_solution_analytical.png", fig_ana)

    fig_logerr = Figure(resolution=(500,400))
    ax3 = Axis(fig_logerr[1,1], title="Log10 Error Magnitude", xlabel="x", ylabel="y", aspect=DataAspect())
    log_err_grid = log10.(abs.(err_grid))
    hm3 = heatmap!(ax3, xg, yg, log_err_grid, colormap = :viridis)
    Colorbar(fig_logerr[1,2], hm3, label = "Log10 |Error|")
    display(fig_logerr)
    save("disk_steepk_log_error.png", fig_logerr)
end

function plot_solution_and_error(; nx=320, ny=320)
    lx = 1.0
    ly = 1.0
    center = (0.5, 0.5)
    R = 0.5
    k_min = 1.0
    k_max = 10.0
    r0 = 0.25
    eps = 0.025
    m = 2

    mesh = Penguin.Mesh((nx, ny), (lx, ly), (0.0, 0.0))
    r_of_xy = (x, y) -> sqrt((x - center[1])^2 + (y - center[2])^2)
    body = (x, y, _=0) -> r_of_xy(x, y)- R
    capacity = Capacity(body, mesh; method="ImplicitIntegration")
    operator = DiffusionOps(capacity)

    u_exact = (x, y) -> exact_solution(x, y, r_of_xy(x, y), R, m, center)
    f = (x, y, _=0) -> forcing_term(x, y, r_of_xy(x, y), R, m, center; k_min=k_min, k_max=k_max, r0=r0, eps=eps)
    D = (x, y, _=0) -> k_profile(r_of_xy(x, y); k_min=k_min, k_max=k_max, r0=r0, eps=eps)

    bc_outer = Dirichlet(0.0)
    bc_b = BorderConditions(Dict(
        :left   => bc_outer,
        :right  => bc_outer,
        :top    => bc_outer,
        :bottom => bc_outer
    ))

    phase = Phase(capacity, operator, f, D)
    solver = DiffusionSteadyMonoVariable(phase, bc_b, Dirichlet((x,y,_=0)->u_exact(x,y)))
    solve_DiffusionSteadyMono!(solver; method=Base.:\)

    plot_solutionn(solver, mesh, capacity, u_exact)
end

plot_solution_and_error()
