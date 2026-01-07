using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using CSV
using Test

"""
Johansen–Colella Problem 1: constant-coefficient Poisson equation inside a
star-shaped domain defined by
    Ω = {(r, θ) : r ≤ 0.30 + 0.15 cos(6θ)}.
We solve Δϕ = 7 r² cos(3θ) with Dirichlet data taken from the exact solution
ϕ(r, θ) = r⁴ cos(3θ). Convergence is measured on uniform Cartesian meshes using
Penguin's implicit-integration capacity.
"""

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

center_default() = (0.49, 0.5)
base_radius() = 0.30
perturb_radius() = 0.15
star_max_radius() = base_radius() + perturb_radius()

function star_level_set(x, y, _=0)
    dx = x - center_default()[1]
    dy = y - center_default()[2]
    θ = atan(dy, dx)
    r = sqrt(dx^2 + dy^2)
    return r - (base_radius() + perturb_radius() * cos(6θ))
end

function polar_features(x::Float64, y::Float64, center::Tuple{Float64,Float64})
    dx = x - center[1]
    dy = y - center[2]
    r2 = dx^2 + dy^2
    θ = atan(dy, dx)
    return r2, θ
end

function exact_phi(x::Float64, y::Float64, center=center_default())
    r2, θ = polar_features(x, y, center)
    return (r2^2) * cos(3θ)
end

function forcing_problem1(x::Float64, y::Float64, center=center_default())
    r2, θ = polar_features(x, y, center)
    return -7.0 * r2 * cos(3θ)
end

function run_star_poisson_convergence(
    nx_list::Vector{Int},
    ny_list::Vector{Int},
    center::Tuple{Float64,Float64},
    u_exact::Function;
    source::Function,
    diffusivity::Function = (x,y,_)->1.0,
    lx::Float64 = 1.0,
    ly::Float64 = 1.0,
    norm::Int = 2,
    relative::Bool = false
)
    h_vals = Float64[]
    err_vals = Float64[]
    err_full = Float64[]
    err_cut = Float64[]
    err_empty = Float64[]
    inside_cells = Int[]
    inside_cells_by_dim = Vector{Vector{Int}}()

    body = star_level_set

    for (nx, ny) in zip(nx_list, ny_list)
        mesh = Penguin.Mesh((nx, ny), (lx, ly), (0.0, 0.0))
        capacity = Capacity(body, mesh; method="ImplicitIntegration")
        operator = DiffusionOps(capacity)

        bc_outer = Dirichlet(0.0)
        bc_b = BorderConditions(Dict(
            :left   => bc_outer,
            :right  => bc_outer,
            :top    => bc_outer,
            :bottom => bc_outer
        ))
        phase = Phase(capacity, operator, source, diffusivity)
        solver = DiffusionSteadyMono(phase, bc_b, Dirichlet((x,y,_=0)->u_exact(x,y)))
        solve_DiffusionSteadyMono!(solver; method=Base.:\)

        _, _, global_err, full_err, cut_err, empty_err =
            check_convergence(u_exact, solver, capacity, norm, relative)

        push!(h_vals, min(lx / nx, ly / ny))
        push!(err_vals, global_err)
        push!(err_full, full_err)
        push!(err_cut, cut_err)
        push!(err_empty, empty_err)
        push!(inside_cells, count_inside_cells(capacity))
        Δx = lx / nx
        Δy = ly / ny
        coverage_x = ceil(Int, 2 * star_max_radius() / Δx)
        coverage_y = ceil(Int, 2 * star_max_radius() / Δy)
        push!(inside_cells_by_dim, [coverage_x, coverage_y])
    end

    return (
        h_vals = h_vals,
        err_vals = err_vals,
        err_full_vals = err_full,
        err_cut_vals = err_cut,
        err_empty_vals = err_empty,
        inside_cells = inside_cells,
        inside_cells_by_dim = inside_cells_by_dim,
        orders = compute_orders(h_vals, err_vals, err_full, err_cut),
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

function main(; csv_path=nothing, nx_list=nothing, ny_list=nothing)
    center = center_default()
    nx_vals = isnothing(nx_list) ? [16, 32, 64, 128, 256, 512] : nx_list
    ny_vals = isnothing(ny_list) ? nx_vals : ny_list
    u_exact = (x,y) -> exact_phi(x,y,center)

    data = run_star_poisson_convergence(
        nx_vals, ny_vals, center, u_exact;
        source = (x,y,_)->forcing_problem1(x,y,center),
        diffusivity = (x,y,_)->1.0,
        lx = 1.0,
        ly = 1.0,
        norm = 2
    )

    csv_info = write_convergence_csv("JohansenColella_P1", data; csv_path=csv_path)
    return (data = data, csv_path = csv_info.csv_path, table = csv_info.table)
end

results = main()

@testset "Johansen-Colella Problem 1" begin
    orders = results.data.orders
    @test !isnan(orders.all)
    @test orders.all > 1.0
    @test length(results.data.h_vals) == length(results.data.err_vals)
    @test results.data.h_vals[1] > results.data.h_vals[end]
    @test minimum(results.data.err_vals) < maximum(results.data.err_vals)
    @test isfile(results.csv_path)
    println(orders)
end


# Using cairomakie, plot the solution and the error on the finest mesh
using CairoMakie
function plot_solutionn(solver, mesh, body::Function, capacity; state_i=1)
   
    # Use check_convergence to obtain analytical and numerical cell-centered fields
    u_ana, u_num, _, _, _, _ = check_convergence((x,y)->exact_phi(x,y), solver, capacity)

    # Mask inactive cells (empty cells) using capacity.cell_types
    cell_types = capacity.cell_types
    u_num_promasked = copy(u_num)
    u_ana_promasked = copy(u_ana)
    u_num_promasked[cell_types .== 0] .= NaN
    u_ana_promasked[cell_types .== 0] .= NaN

    # Prepare plotting grid and reshape arrays to (ny, nx) as used by heatmap
    xg = mesh.nodes[1]
    yg = mesh.nodes[2]
    nx_plot = length(xg)
    ny_plot = length(yg)


    sol_grid = reshape(u_num_promasked, (nx_plot, ny_plot))'
    ana_grid = reshape(u_ana_promasked, (nx_plot, ny_plot))'
    err_grid = sol_grid .- ana_grid

    # Solution figure
    fig_sol = Figure(resolution=(500,400))
    ax1 = Axis(fig_sol[1,1], title="Solution (numerical)", xlabel="x", ylabel="y", aspect=DataAspect())
    hm1 = heatmap!(ax1, xg, yg, sol_grid, colormap = :viridis)
    Colorbar(fig_sol[1,2], hm1, label = "Value")
    display(fig_sol)
    save("johansen_colella_problem1_solution.png", fig_sol)

    # Error figure (separate)
    fig_err = Figure(resolution=(500,400))
    ax3 = Axis(fig_err[1,1], title="Error (numerical - analytical)", xlabel="x", ylabel="y", aspect=DataAspect())
    hm3 = heatmap!(ax3, xg, yg, err_grid, colormap = :viridis)
    Colorbar(fig_err[1,2], hm3, label = "Error")
    display(fig_err)

    # Log error figure (separate)
    fig_logerr = Figure(resolution=(500,400))
    ax4 = Axis(fig_logerr[1,1], title="Log10 Error Magnitude", xlabel="x", ylabel="y", aspect=DataAspect())
    log_err_grid = log10.(abs.(err_grid))
    hm4 = heatmap!(ax4, xg, yg, log_err_grid, colormap = :viridis)
    Colorbar(fig_logerr[1,2], hm4, label = "Log10 |Error|")
    display(fig_logerr)
    save("johansen_colella_problem1_log_error.png", fig_logerr)
end
function plot_solution_and_error(; nx=512, ny=512, center=center_default())
    lx = 1.0
    ly = 1.0
    mesh = Penguin.Mesh((nx, ny), (lx, ly), (0.0, 0.0))
    capacity = Capacity(star_level_set, mesh; method="ImplicitIntegration")
    operator = DiffusionOps(capacity) 
    bc_outer = Dirichlet(0.0)
    bc_b = BorderConditions(Dict(
        :left   => bc_outer,
        :right  => bc_outer,
        :top    => bc_outer,
        :bottom => bc_outer
    ))
    phase = Phase(capacity, operator,
        (x,y,_)->forcing_problem1(x,y,center),
        (x,y,_)->1.0
    )
    solver = DiffusionSteadyMono(phase, bc_b, Dirichlet((x,y,_=0)->exact_phi(x,y,center)))
    solve_DiffusionSteadyMono!(solver; method=Base.:\)
    
    plot_solutionn(solver, mesh, star_level_set, capacity)
end

plot_solution_and_error()